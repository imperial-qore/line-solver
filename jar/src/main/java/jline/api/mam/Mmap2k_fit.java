/**
 * @file MMAP(2,K) closed-form fitting: a marked MAP of second order, K classes
 *
 * The fit factorizes into two independent inverse problems:
 *
 * 1. the underlying MAP(2) from (M1, M2, M3, GAMMA), through the acyclic
 *    canonical inverse of amap2_fit_gamma. Order two loses nothing here: every
 *    MAP(2) is equivalent to one of the two acyclic canonical forms;
 * 2. the marking, which is LINEAR in the class characteristics. Writing the
 *    marking as fractions of the aggregate D1,
 *      form 1: D1c = D1 .* [q1c 0; q2c q3c]
 *      form 2: D1c = D1 .* [0 q1c; q2c q3c]
 *    the map (q1c, q2c, q3c) -> (p_c, p_c F_c, p_c B_c) is linear with a 3x3
 *    matrix that depends only on (h1, h2, r1, r2) and is the SAME for every
 *    class. Its inverse is pre-computed below and applied per class, so the
 *    cost does not grow with K and no quadratic program is involved.
 *
 * The inverse was derived symbolically in SageMath, see
 * io/sage/proofs/mmap2k_marking_inverse.py. It is singular exactly on the six
 * degenerate loci r1 in {0,1}, r2 in {0,1}, h1-h2+h2*r1 = 0 and, per form,
 * h1*r2-h2 = 0 or h1*r1*r2-h1*r1+h1-h2 = 0, which are the branches
 * Mamap2m_fit_fb_multiclass handles; there, and whenever the closed form leaves
 * the unit box, this class falls back to that quadratic program.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import java.util.List;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Mmap2k_fit {
    private Mmap2k_fit() {}

    private static final double DEGENTOL = 1e-8;
    private static final double FEASTOL = 1e-8;

    /** Result of a fit: the MMAP and whether the closed form matched exactly. */
    public static final class Result {
        public final MatrixCell mmap;
        public final boolean exact;

        public Result(MatrixCell mmap, boolean exact) {
            this.mmap = mmap;
            this.exact = exact;
        }
    }

    /**
     * Fits an MMAP(2,K) to the inter-arrival moments, the autocorrelation decay
     * rate and the per-class probabilities, forward moments and backward moments.
     *
     * @param M1 first raw moment of the inter-arrival times
     * @param M2 second raw moment
     * @param M3 third raw moment
     * @param GAMMA autocorrelation decay rate
     * @param P class probabilities, summing to one
     * @param F first-order forward moments, with sum_c P[c] F[c] = M1
     * @param B first-order backward moments, with sum_c P[c] B[c] = M1
     */
    public static Result mmap2k_fit(double M1, double M2, double M3, double GAMMA,
                                    double[] P, double[] F, double[] B) {
        int K = P.length;
        if (F.length != K || B.length != K) {
            throw new IllegalArgumentException("mmap2k_fit: P, F and B must have the same length");
        }

        List<MatrixCell> candidates = Amap2_fit_gamma.amap2_fitall_gamma(M1, M2, M3, GAMMA);
        if (candidates.isEmpty()) {
            jline.util.Pair<MatrixCell, List<MatrixCell>> fit =
                Amap2_fit_gamma.amap2_fit_gamma(M1, M2, M3, GAMMA);
            candidates = new java.util.ArrayList<MatrixCell>();
            candidates.add(fit.getLeft());
        }

        Matrix bestD0 = null;
        Matrix bestD1 = null;
        int bestForm = 0;
        double[][] bestQ = null;
        double bestErr = Double.POSITIVE_INFINITY;

        for (MatrixCell cand : candidates) {
            Matrix D0 = cand.get(0);
            Matrix D1 = cand.get(1);
            if (D0.getNumRows() != 2 || D0.get(1, 0) != 0.0) {
                continue;
            }
            int form;
            if (D1.get(0, 1) == 0.0) {
                form = 1;
            } else if (D1.get(0, 0) == 0.0) {
                form = 2;
            } else {
                continue;
            }
            double h1 = -1.0 / D0.get(0, 0);
            double h2 = -1.0 / D0.get(1, 1);
            double r1 = D0.get(0, 1) * h1;
            double r2 = D1.get(1, 1) * h2;

            double[][] q = new double[3][K];
            boolean ok = true;
            for (int c = 0; c < K && ok; c++) {
                double[] qc = markingInverse(form, h1, h2, r1, r2, P[c], F[c], B[c]);
                if (qc == null) {
                    ok = false;
                } else {
                    q[0][c] = qc[0];
                    q[1][c] = qc[1];
                    q[2][c] = qc[2];
                }
            }
            if (!ok) {
                continue;
            }

            double viol = 0.0;
            for (int j = 0; j < 3; j++) {
                double sum = 0.0;
                for (int c = 0; c < K; c++) {
                    viol = Math.max(viol, -q[j][c]);
                    viol = Math.max(viol, q[j][c] - 1.0);
                    sum += q[j][c];
                }
                viol = Math.max(viol, Math.abs(sum - 1.0));
            }
            if (viol < bestErr) {
                bestErr = viol;
                bestD0 = D0;
                bestD1 = D1;
                bestForm = form;
                bestQ = q;
            }
        }

        if (bestQ != null && bestErr <= FEASTOL) {
            MatrixCell out = new MatrixCell(2 + K);
            out.set(0, bestD0.copy());
            out.set(1, bestD1.copy());
            for (int c = 0; c < K; c++) {
                Matrix Dc = Matrix.zeros(2, 2);
                double q1 = Math.min(Math.max(bestQ[0][c], 0.0), 1.0);
                double q2 = Math.min(Math.max(bestQ[1][c], 0.0), 1.0);
                double q3 = Math.min(Math.max(bestQ[2][c], 0.0), 1.0);
                if (bestForm == 1) {
                    Dc.set(0, 0, bestD1.get(0, 0) * q1);
                } else {
                    Dc.set(0, 1, bestD1.get(0, 1) * q1);
                }
                Dc.set(1, 0, bestD1.get(1, 0) * q2);
                Dc.set(1, 1, bestD1.get(1, 1) * q3);
                out.set(2 + c, Dc);
            }
            return new Result(out, true);
        }

        // degenerate underlying form, or characteristics outside the feasible set
        MatrixCell fallback = Mamap2m_fit.mamap2m_fit_gamma_fb(M1, M2, M3, GAMMA, P, F, B);
        return new Result(fallback, false);
    }

    /**
     * Pre-computed inverse of the per-class linear map; null on the singular loci.
     */
    private static double[] markingInverse(int form, double h1, double h2, double r1, double r2,
                                           double p, double Fc, double Bc) {
        if (form == 1) {
            double d1 = (r2 - 1) * (r1 - 1) * (h1 * r2 - h2);
            double d2 = (r2 - 1) * r1;
            double d3 = r1 * r2 * (h1 + h2 * r1 - h2);
            if (Math.abs(d1) < DEGENTOL || Math.abs(d2) < DEGENTOL || Math.abs(d3) < DEGENTOL
                    || Math.abs(h1 + h2 * r1 - h2) < DEGENTOL || Math.abs(h1 * r2 - h2) < DEGENTOL) {
                return null;
            }
            double W = r1 * r2 - r2 + 1;
            double q1 = p * W * ((h1 * r2 - h1 - h2) + Bc) / d1;
            double q2 = p * W * ((h1 * h1 * (r2 - 1) + h1 * h2 * r1 * (r2 - 1) - h2 * h2 * r1)
                        / ((h1 + h2 * r1 - h2) * (h1 * r2 - h2))
                        - Fc / (h1 + h2 * r1 - h2)
                        + Bc / (h1 * r2 - h2)) / d2;
            double q3 = p * W * ((h1 + h2 * r1) - Fc) / d3;
            return new double[]{q1, q2, q3};
        } else {
            double U = h1 * r1 * r2 - h1 * r1 + h1 - h2;
            double d1 = (r2 - 1) * (r1 - 1) * U;
            double d2 = (r2 - 1) * (h1 + h2 * r1 - h2);
            if (Math.abs(d1) < DEGENTOL || Math.abs(d2) < DEGENTOL || Math.abs(r2) < DEGENTOL
                    || Math.abs(U) < DEGENTOL || Math.abs(h1 + h2 * r1 - h2) < DEGENTOL) {
                return null;
            }
            double V = r1 * r2 - r1 - r2 + 2;
            double q1 = p * V * ((h1 * r1 * r2 - h1 * r1 - h2) + Bc) / d1;
            double q2 = p * V * (h2 - Fc) / d2;
            double q3 = p * V * ((h1 * h1 + h1 * h2 * r1 * r2 - h2 * h2)
                        / ((h1 + h2 * r1 - h2) * U)
                        - Fc / (h1 + h2 * r1 - h2)
                        - Bc / U) / r2;
            return new double[]{q1, q2, q3};
        }
    }
}
