`timescale 1ns/1ps

// ============================================================================
// systolic_matmul : 4 x 4 output-stationary systolic array, C = A x B
//
//   A rows enter from the left (row i delayed by i cycles) and move right.
//   B cols enter from the top  (col j delayed by j cycles) and move down.
//   PE(i,j) keeps C[i][j] locally; at compute cycle t it uses k = t - i - j.
//
// Timing contract (matches the supplied testbench):
//   edge E0      : start && !busy  -> latch m/n/k, clear all PEs, busy = 1, t = 0
//   cycles 0..T-1: PEs enabled; E(t+1) registers the cycle-t MAC in every PE
//   edge E_T     : last MAC registered, busy -> 0, done -> 1 (one-cycle pulse)
//   where T = K + M + N - 2 compute cycles.
// ============================================================================
module systolic_matmul #(
    parameter int ARRAY_DIM = 4,
    parameter int DATA_W    = 16,
    parameter int ACC_W     = 32,
    parameter int MAX_K     = 8,
    parameter int DIM_W     = $clog2(ARRAY_DIM + 1),
    parameter int K_W       = $clog2(MAX_K + 1)
)(
    input  logic                     clk,
    input  logic                     rst_n,

    input  logic                     start,
    input  logic [DIM_W-1:0]         m,
    input  logic [DIM_W-1:0]         n,
    input  logic [K_W-1:0]           k,

    input  logic signed [DATA_W-1:0]
                                     mat_a [0:ARRAY_DIM-1][0:MAX_K-1],
    input  logic signed [DATA_W-1:0]
                                     mat_b [0:MAX_K-1][0:ARRAY_DIM-1],

    output logic                     busy,
    output logic                     done,

    output logic signed [ACC_W-1:0]
                                     mat_c [0:ARRAY_DIM-1][0:ARRAY_DIM-1]
);

    // Valid commands satisfy:
    //   1 <= m <= ARRAY_DIM,  1 <= n <= ARRAY_DIM,  1 <= k <= MAX_K
    // The testbench keeps command inputs and matrices stable until done.

    // Longest possible run: K + M + N - 2 with every dimension at its maximum.
    localparam int T_MAX = MAX_K + 2*ARRAY_DIM - 2;
    localparam int T_W   = $clog2(T_MAX + 1);

    // ------------------------------------------------------------------
    // 1. Controller
    // ------------------------------------------------------------------
    logic [DIM_W-1:0] m_r, n_r;     // latched command
    logic [K_W-1:0]   k_r;
    logic [T_W-1:0]   t;            // zero-based compute-cycle index
    logic             accept;       // command-acceptance edge
    int               t_last;       // index of the final compute cycle

    assign accept = start && !busy;

    // Final compute cycle index = T - 1 = K + M + N - 3.
    // Arithmetic is done in int so that nothing wraps as an unsigned value.
    assign t_last = int'(k_r) + int'(m_r) + int'(n_r) - 3;

    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            busy <= 1'b0;
            done <= 1'b0;
            t    <= '0;
            m_r  <= '0;
            n_r  <= '0;
            k_r  <= '0;
        end else begin
            done <= 1'b0;                       // default: done is a 1-cycle pulse

            if (accept) begin                   // setup edge (not a compute cycle)
                m_r  <= m;
                n_r  <= n;
                k_r  <= k;
                t    <= '0;
                busy <= 1'b1;
            end else if (busy) begin
                if (int'(t) == t_last) begin    // this edge registers the last MAC
                    busy <= 1'b0;
                    done <= 1'b1;
                    t    <= '0;
                end else begin
                    t <= t + 1'b1;
                end
            end
        end
    end

    // All PEs are cleared on the acceptance edge and step while busy.
    logic pe_clear, pe_enable;
    assign pe_clear  = accept;
    assign pe_enable = busy;

    // ------------------------------------------------------------------
    // 2. Skewed boundary feeders
    //    A_left[i] = A[i][t-i] when 0 <= t-i < K (and row i is active), else 0
    //    B_top [j] = B[t-j][j] when 0 <= t-j < K (and col j is active), else 0
    //    Rows i >= M and columns j >= N are fed zeros: those PEs produce no
    //    useful output, so they are kept quiet.
    // ------------------------------------------------------------------
    logic signed [DATA_W-1:0] a_left [0:ARRAY_DIM-1];
    logic signed [DATA_W-1:0] b_top  [0:ARRAY_DIM-1];

    always_comb begin
        for (int i = 0; i < ARRAY_DIM; i++) begin
            if (busy && (i < int'(m_r)) &&
                (int'(t) - i >= 0) && (int'(t) - i < int'(k_r)))
                a_left[i] = mat_a[i][int'(t) - i];
            else
                a_left[i] = '0;
        end
    end

    always_comb begin
        for (int j = 0; j < ARRAY_DIM; j++) begin
            if (busy && (j < int'(n_r)) &&
                (int'(t) - j >= 0) && (int'(t) - j < int'(k_r)))
                b_top[j] = mat_b[int'(t) - j][j];
            else
                b_top[j] = '0;
        end
    end

    // ------------------------------------------------------------------
    // 3. 4 x 4 PE array with neighbour-to-neighbour links
    //    a_bus[i][j] = A input of PE(i,j); a_bus[i][j+1] = its A output
    //    b_bus[i][j] = B input of PE(i,j); b_bus[i+1][j] = its B output
    //    (column ARRAY_DIM of a_bus and row ARRAY_DIM of b_bus leave the array)
    // ------------------------------------------------------------------
    wire signed [DATA_W-1:0] a_bus [0:ARRAY_DIM-1][0:ARRAY_DIM];
    wire signed [DATA_W-1:0] b_bus [0:ARRAY_DIM][0:ARRAY_DIM-1];

    genvar gi, gj;
    generate
        for (gi = 0; gi < ARRAY_DIM; gi++) begin : g_left_edge
            assign a_bus[gi][0] = a_left[gi];
        end

        for (gj = 0; gj < ARRAY_DIM; gj++) begin : g_top_edge
            assign b_bus[0][gj] = b_top[gj];
        end

        for (gi = 0; gi < ARRAY_DIM; gi++) begin : g_row
            for (gj = 0; gj < ARRAY_DIM; gj++) begin : g_col
                systolic_pe #(
                    .DATA_W (DATA_W),
                    .ACC_W  (ACC_W)
                ) u_pe (
                    .clk    (clk),
                    .rst_n  (rst_n),
                    .clear  (pe_clear),
                    .enable (pe_enable),
                    .a_in   (a_bus[gi][gj]),
                    .b_in   (b_bus[gi][gj]),
                    .a_out  (a_bus[gi][gj+1]),     // to right neighbour
                    .b_out  (b_bus[gi+1][gj]),     // to lower neighbour
                    .psum   (mat_c[gi][gj])        // C[i][j] stays here
                );
            end
        end
    endgenerate

endmodule
