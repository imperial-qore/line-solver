/**
 * @file Quasi-Birth-Death process R and G matrix computation
 *
 * Computes fundamental R and G matrices for QBD analysis of MAP/MAP/1 queues.
 * Essential for solving infinite-dimensional Markov chains with structured transitions.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import java.util.Map;
import java.util.Objects;

import jline.lib.smc.QBD_CR;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Qbd_rg {
    private Qbd_rg() {}

    /**
     * Result class for qbd_rg containing R, G, B, L, F, U.
     */
    public static final class QbdRgResult {
        public final Matrix R;
        public final Matrix G;
        public final Matrix B;
        public final Matrix L;
        public final Matrix F;
        public final Matrix U;

        public QbdRgResult(Matrix R, Matrix G, Matrix B, Matrix L, Matrix F, Matrix U) {
            this.R = R;
            this.G = G;
            this.B = B;
            this.L = L;
            this.F = F;
            this.U = U;
        }

        public Matrix component1() { return R; }
        public Matrix component2() { return G; }
        public Matrix component3() { return B; }
        public Matrix component4() { return L; }
        public Matrix component5() { return F; }
        public Matrix component6() { return U; }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (!(o instanceof QbdRgResult)) return false;
            QbdRgResult that = (QbdRgResult) o;
            return Objects.equals(R, that.R) && Objects.equals(G, that.G)
                    && Objects.equals(B, that.B) && Objects.equals(L, that.L)
                    && Objects.equals(F, that.F) && Objects.equals(U, that.U);
        }

        @Override
        public int hashCode() {
            return Objects.hash(R, G, B, L, F, U);
        }
    }

    /**
     * Result from simplified QBD CR solver (kept for backward compatibility).
     */
    public static final class QbdCrResult {
        public final Matrix G;
        public final Matrix R;
        public final Matrix U;

        public QbdCrResult(Matrix G, Matrix R, Matrix U) {
            this.G = G;
            this.R = R;
            this.U = U;
        }
    }

    /**
     * Compute R and G matrices for MAP/MAP/1 queue using QBD approach.
     *
     * @param MAPa Arrival MAP (MatrixCell with D0 and D1)
     * @param MAPs Service MAP (MatrixCell with D0 and D1)
     * @param util Optional utilization parameter for scaling (may be null)
     * @return QbdRgResult containing R, G, B, L, F, and U matrices
     */
    public static QbdRgResult qbd_rg(MatrixCell MAPa, MatrixCell MAPs, Double util) {
        // int na = MAPa.get(0).getNumRows(); // currently unused (kept from Kotlin)
        // int ns = MAPs.get(0).getNumRows();

        MatrixCell scaledMAPs = MAPs;
        if (util != null) {
            double lambdaA = Map_lambda.map_lambda(MAPa);
            scaledMAPs = Map_scale.map_scale(MAPs, util / lambdaA);
        }

        // Construct QBD matrices using Kronecker sum for L (matching MATLAB krons)
        int na = MAPa.get(0).getNumRows();
        int ns = MAPs.get(0).getNumRows();
        Matrix IA = Matrix.eye(na);
        Matrix IS = Matrix.eye(ns);
        Matrix F = MAPa.get(1).kron(IS);
        Matrix L = MAPa.get(0).kron(IS).add(IA.kron(scaledMAPs.get(0)));  // Kronecker sum
        Matrix B = IA.kron(scaledMAPs.get(1));

        // Solve for G, R, U using SMC library Cyclic Reduction
        Map<String, Matrix> qbdResult = QBD_CR.QBD_CR(B, L, F, null, null, null, null);
        Matrix G = qbdResult.get("G");
        Matrix R = qbdResult.get("R");
        Matrix U = qbdResult.get("U");
        if (G == null) throw new RuntimeException("QBD_CR failed to compute G matrix");
        if (R == null) throw new RuntimeException("QBD_CR failed to compute R matrix");
        if (U == null) throw new RuntimeException("QBD_CR failed to compute U matrix");

        return new QbdRgResult(R, G, B, L, F, U);
    }

    public static QbdRgResult qbd_rg(MatrixCell MAPa, MatrixCell MAPs) {
        return qbd_rg(MAPa, MAPs, null);
    }
}
