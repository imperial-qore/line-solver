/**
 * @file Square-root Non-iterative (SQNI) approximate MVA
 *
 * Implements the Square-root Non-iterative method for analyzing multi-class closed
 * queueing networks. Provides efficient approximation by reducing multi-class networks to
 * single-queue representations with interpolation-based corrections.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.mva;

import java.util.Objects;

import jline.util.matrix.Matrix;

public final class Pfqn_sqni {
    private Pfqn_sqni() {}

    /**
     * Result of pfqn_sqni containing Q, U and X matrices.
     */
    public static final class PfqnSqniResult {
        public final Matrix Q;
        public final Matrix U;
        public final Matrix X;

        public PfqnSqniResult(Matrix Q, Matrix U, Matrix X) {
            this.Q = Q;
            this.U = U;
            this.X = X;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (!(o instanceof PfqnSqniResult)) return false;
            PfqnSqniResult that = (PfqnSqniResult) o;
            return Objects.equals(Q, that.Q) && Objects.equals(U, that.U) && Objects.equals(X, that.X);
        }

        @Override
        public int hashCode() {
            return Objects.hash(Q, U, X);
        }

        @Override
        public String toString() {
            return "PfqnSqniResult(Q=" + Q + ", U=" + U + ", X=" + X + ")";
        }
    }

    public static PfqnSqniResult pfqn_sqni(Matrix N, Matrix L, Matrix Z) {
        int queueIdx = 0;
        int C = L.length();
        double Nt = N.elementSum();
        Matrix Q = new Matrix(2, C);
        Matrix U = new Matrix(2, C);
        Matrix X = new Matrix(1, C);

        if (Nt <= 0) {
            return new PfqnSqniResult(Q, U, X);
        }

        if (Nt == 1.0) {
            for (int r = 0; r < C; r++) {
                double Xr = N.get(0, r) / (Z.get(r) + L.get(r));
                X.set(0, r, Xr);
                U.set(queueIdx, r, Xr * L.get(r));
                Q.set(queueIdx, r, Xr * L.get(r));
            }
        } else {
            // A Z=0 class (self-looping) has no delay to interpolate through: its
            // queue length is its whole population and it is solved after the loop.
            for (int r = 0; r < C; r++) {
                if (Z.get(r) == 0.0) {
                    Q.set(queueIdx, r, N.get(0, r));
                }
            }

            for (int r = 0; r < C; r++) {
                if (Z.get(r) == 0.0) {
                    continue;
                }
                double Nr = N.get(0, r);
                double Lr = L.get(r);
                double Zr = Z.get(r);

                Matrix Nvec_1r = N.copy();
                Nvec_1r.set(0, r, Nvec_1r.get(0, r) - 1);

                double sumN = N.elementSum();
                // sumBrPart runs over EVERY class, class r included, as in MATLAB
                // pfqn_sqni; skipping r shifted X by 1.6% on a 2-class model.
                double sumBrPart = 0.0;
                for (int i = 0; i < C; i++) {
                    double Zi = Z.get(i);
                    double Li = L.get(i);
                    double Ni = Nvec_1r.get(0, i);
                    sumBrPart += Zi * Ni / (Zi + Li + Li * (sumN - 2));
                }

                Matrix BrVec = new Matrix(1, C);
                for (int i = 0; i < C; i++) {
                    double Zi = Z.get(i);
                    double Li = L.get(i);
                    double Ni = N.get(0, i);
                    double denom = Zi + Li + Li * (sumN - 1 - sumBrPart);
                    BrVec.set(0, i, Ni / denom * Zi);
                }

                double BrSum = 0.0;
                for (int i = 0; i < C; i++) {
                    if (i != r) {
                        BrSum += BrVec.get(0, i);
                    }
                }

                double Br = Lr * BrSum;
                double Xr;
                if (Lr == 0.0) {
                    Xr = Nr / Zr;
                } else {
                    double sqrtTerm = Math.sqrt(Br * Br - 2 * Br * Lr * Nt - 2 * Br * Zr
                            + Lr * Lr * Nt * Nt + 2 * Lr * Nt * Zr - 4 * Nr * Lr * Zr + Zr * Zr);
                    Xr = (Zr - sqrtTerm - Br + Lr * Nt) / (2 * Lr * Zr);
                }

                X.set(0, r, Xr);
                U.set(queueIdx, r, Xr * Lr);
                Q.set(queueIdx, r, Nr - Xr * Zr);
            }
        }

        for (int r = 0; r < C; r++) {
            if (Z.get(r) == 0.0) {
                double denom = L.get(r) * (1 + Q.sumRows(queueIdx));
                double Xr = N.get(0, r) / denom;
                X.set(0, r, Xr);
                U.set(queueIdx, r, Xr * L.get(r));
                Q.set(queueIdx, r, N.get(0, r) - Xr * Z.get(r));
            }
        }

        return new PfqnSqniResult(Q, U, X);
    }
}
