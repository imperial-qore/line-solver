/*
 * Compute T-matrix using either Sylvester or NARE method
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
package jline.lib.fjcodes;

import java.util.Objects;

import jline.util.matrix.Matrix;

public final class ComputeT {
    private ComputeT() {}

    /**
     * Result of computeT.
     */
    public static final class ComputeTResult {
        public final Matrix T;
        public final Matrix S;
        public final Matrix A_jump;
        public final Matrix S_Arr;
        public final Matrix sum_Ajump;

        public ComputeTResult(Matrix T, Matrix S, Matrix A_jump, Matrix S_Arr, Matrix sum_Ajump) {
            this.T = T;
            this.S = S;
            this.A_jump = A_jump;
            this.S_Arr = S_Arr;
            this.sum_Ajump = sum_Ajump;
        }

        public Matrix component1() { return T; }
        public Matrix component2() { return S; }
        public Matrix component3() { return A_jump; }
        public Matrix component4() { return S_Arr; }
        public Matrix component5() { return sum_Ajump; }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (!(o instanceof ComputeTResult)) return false;
            ComputeTResult that = (ComputeTResult) o;
            return Objects.equals(T, that.T) && Objects.equals(S, that.S)
                    && Objects.equals(A_jump, that.A_jump) && Objects.equals(S_Arr, that.S_Arr)
                    && Objects.equals(sum_Ajump, that.sum_Ajump);
        }

        @Override
        public int hashCode() {
            return Objects.hash(T, S, A_jump, S_Arr, sum_Ajump);
        }
    }

    /**
     * Compute T-matrix using specified method.
     *
     * The T-matrix is computed using either the Sylvester equation approach
     * or the NARE (Nonsymmetric Algebraic Riccati Equation) method.
     *
     * @param arrival Arrival process (lambda0, lambda1)
     * @param services Service process for single subtask
     * @param service_h Service representation for 2-node job
     * @param C Capacity parameter
     * @param tMode Method to use: "NARE" (default) or "Sylvest"
     * @return ComputeTResult with T, S, A_jump, S_Arr, sum_Ajump
     */
    public static ComputeTResult computeT(FJArrival arrival, FJService services, FJServiceH service_h, int C, String tMode) {
        // Build S and A_jump matrices
        SAResult saResult = BuildSA.build_SA(services, service_h, C);
        Matrix S = saResult.getS();
        Matrix A_jump = saResult.getA_jump();

        int d0 = arrival.lambda0.getNumRows();
        Matrix S_Arr = S.kron(Matrix.eye(d0));
        Matrix A_jump_Arr = A_jump.kron(Matrix.eye(d0));

        // Compute T-matrix using specified method
        Matrix T;
        if (tMode != null && tMode.toLowerCase().contains("sylvest")) {
            // Sylvester equation approach
            T = ComputeT_Sylvester.computeT_Sylvester(arrival.lambda0, arrival.lambda1, S_Arr, A_jump);
        } else {
            // NARE method (default)
            T = ComputeT_NARE.computeT_NARE(arrival.lambda0, arrival.lambda1, S_Arr, A_jump);
        }

        // Compute sum of A_jump_Arr rows
        Matrix sum_Ajump = new Matrix(A_jump_Arr.getNumRows(), 1);
        for (int i = 0; i < A_jump_Arr.getNumRows(); i++) {
            double sum = 0.0;
            for (int j = 0; j < A_jump_Arr.getNumCols(); j++) {
                sum += A_jump_Arr.get(i, j);
            }
            sum_Ajump.set(i, 0, sum);
        }

        return new ComputeTResult(T, S, A_jump, S_Arr, sum_Ajump);
    }

    public static ComputeTResult computeT(FJArrival arrival, FJService services, FJServiceH service_h, int C) {
        return computeT(arrival, services, service_h, C, "NARE");
    }
}
