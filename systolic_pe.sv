`timescale 1ns/1ps

// ============================================================================
// systolic_pe : one processing element of the output-stationary systolic array
//
//   a_in  ---> [ a_out register ] ---> a_out   (A moves one PE to the right)
//   b_in  ---> [ b_out register ] ---> b_out   (B moves one PE downward)
//   psum  <= psum + a_in * b_in                (C[i][j] stays in this PE)
//
// Per clock edge:
//   rst_n = 0 : asynchronous reset of all state
//   clear = 1 : zero psum AND both forwarding registers (start of a new
//               transaction). Zeroing the forwarding registers matters: at the
//               end of a run the bottom-right PEs still hold the last A/B
//               values, which would otherwise be multiplied into the next run.
//   enable = 1: one compute step - forward operands and accumulate
//   otherwise : hold everything (results stay visible after done)
// clear has priority over enable.
// ============================================================================
module systolic_pe #(
    parameter int DATA_W = 16,
    parameter int ACC_W  = 32
)(
    input  logic                     clk,
    input  logic                     rst_n,

    // Clear the local accumulator before a new matrix operation.
    input  logic                     clear,

    // One compute step. When enable = 1, operands advance by one PE
    // and the local product is accumulated.
    input  logic                     enable,

    input  logic signed [DATA_W-1:0] a_in,
    input  logic signed [DATA_W-1:0] b_in,

    output logic signed [DATA_W-1:0] a_out,
    output logic signed [DATA_W-1:0] b_out,
    output logic signed [ACC_W-1:0]  psum
);

    // Full-precision signed product (16b x 16b -> 32b).
    // Both operands are declared signed, so SystemVerilog performs a signed
    // multiply and sign-extends the operands to the 2*DATA_W result width.
    // (Mixing in any unsigned operand would silently make it unsigned.)
    logic signed [2*DATA_W-1:0] product;
    assign product = a_in * b_in;

    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            a_out <= '0;
            b_out <= '0;
            psum  <= '0;
        end else if (clear) begin
            a_out <= '0;
            b_out <= '0;
            psum  <= '0;
        end else if (enable) begin
            a_out <= a_in;            // pass A to the right neighbour
            b_out <= b_in;            // pass B to the lower neighbour
            psum  <= psum + product;  // signed accumulate (wraps mod 2^ACC_W)
        end
    end

endmodule
