/**
 * @file Distance measures for discrete-time MAPs (D-MAPs)
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;

public final class Dmap_dist {
    private Dmap_dist() {}

    /**
     * Computes the joint PMF inner product of two D-MAPs via recursive discrete Lyapunov equations.
     */
    public static double dmap_geo_mul_sum(Matrix D0A, Matrix D1A, Matrix D0B, Matrix D1B, int L,
                                           Matrix alA, Matrix alB) {
        Matrix Z = Matrix.dlyap(D0B.transpose(), D0A, alB.transpose().mult(alA));
        for (int i = 1; i < L; i++) {
            Z = Matrix.dlyap(D0B.transpose(), D0A, D1B.transpose().mult(Z).mult(D1A));
        }
        int NA = D0A.getNumRows();
        int NB = D0B.getNumRows();
        Matrix dA = Matrix.eye(NA).add(-1.0, D0A).sumRows();
        Matrix dB = Matrix.eye(NB).add(-1.0, D0B).sumRows();
        return dB.transpose().mult(Z).mult(dA).toDouble();
    }

    /**
     * Computes the squared L2 distance between lag-L joint PMFs of two D-MAPs.
     */
    public static double dmap_dist(Matrix D0A, Matrix D1A, Matrix D0B, Matrix D1B, int L) {
        Matrix alA = Dmap_pie.dmap_pie(D0A, D1A);
        Matrix alB = Dmap_pie.dmap_pie(D0B, D1B);
        return dmap_dist(D0A, D1A, D0B, D1B, L, alA, alB);
    }

    public static double dmap_dist(Matrix D0A, Matrix D1A, Matrix D0B, Matrix D1B, int L,
                                    Matrix alA, Matrix alB) {
        return dmap_geo_mul_sum(D0A, D1A, D0A, D1A, L + 1, alA, alA)
                - 2 * dmap_geo_mul_sum(D0A, D1A, D0B, D1B, L + 1, alA, alB)
                + dmap_geo_mul_sum(D0B, D1B, D0B, D1B, L + 1, alB, alB);
    }

    /**
     * Computes the geometric sum for the discrete autocorrelation distance of two D-MAPs.
     */
    public static double dmap_geo_mul_sum_acf(Matrix D0A, Matrix D1A, Matrix D0B, Matrix D1B,
                                               Matrix alA, Matrix alB) {
        int NA = D0A.getNumRows();
        int NB = D0B.getNumRows();
        Matrix D0Ai = Matrix.eye(NA).add(-1.0, D0A).inv();
        Matrix D0Bi = Matrix.eye(NB).add(-1.0, D0B).inv();
        Matrix PAh = D0Ai.mult(D1A).add(-1.0, Matrix.ones(NA, 1).mult(alA));
        Matrix PBh = D0Bi.mult(D1B).add(-1.0, Matrix.ones(NB, 1).mult(alB));
        Matrix M = Matrix.eye(NA * NB).add(-1.0, PBh.transpose().kron(PAh));
        if (Math.abs(M.det()) < 1e-10) {
            return Double.MAX_VALUE;
        }
        Matrix X = Matrix.dlyap(PAh, PBh, D0Ai.sumRows().mult(alB.mult(D0Bi)));
        // Total inner product = sum of all entries of the 1xNB row vector;
        // sumCols().toDouble() kept only the first column (see Map_dist_acf)
        return alA.mult(D0Ai).mult(X).mult(D0Bi).elementSum();
    }

    /**
     * Computes the squared L2 distance between autocorrelation functions of two D-MAPs.
     */
    public static double dmap_dist_acf(Matrix D0A, Matrix D1A, Matrix D0B, Matrix D1B) {
        Matrix alA = Dmap_pie.dmap_pie(D0A, D1A);
        Matrix alB = Dmap_pie.dmap_pie(D0B, D1B);
        return dmap_dist_acf(D0A, D1A, D0B, D1B, alA, alB);
    }

    public static double dmap_dist_acf(Matrix D0A, Matrix D1A, Matrix D0B, Matrix D1B,
                                        Matrix alA, Matrix alB) {
        Matrix D0Ai_A = Matrix.eye(D0A.getNumRows()).add(-1.0, D0A).inv();
        Matrix D0Ai_B = Matrix.eye(D0B.getNumRows()).add(-1.0, D0B).inv();
        Matrix eA = Matrix.ones(D0A.getNumRows(), 1);
        Matrix eB = Matrix.ones(D0B.getNumRows(), 1);
        double muA = alA.mult(D0Ai_A).mult(eA).toDouble();
        double muB = alB.mult(D0Ai_B).mult(eB).toDouble();
        double m2A = alA.mult(Matrix.eye(D0A.getNumRows()).add(1.0, D0A)).mult(D0Ai_A).mult(D0Ai_A).mult(eA).toDouble();
        double m2B = alB.mult(Matrix.eye(D0B.getNumRows()).add(1.0, D0B)).mult(D0Ai_B).mult(D0Ai_B).mult(eB).toDouble();
        double varA = m2A - muA * muA;
        double varB = m2B - muB * muB;
        double cA = (m2A + muA) / 2;
        double cB = (m2B + muB) / 2;
        return (dmap_geo_mul_sum_acf(D0A, D1A, D0A, D1A, alA, alA) - cA * cA) / (varA * varA)
                - 2 * (dmap_geo_mul_sum_acf(D0A, D1A, D0B, D1B, alA, alB) - cA * cB) / (varA * varB)
                + (dmap_geo_mul_sum_acf(D0B, D1B, D0B, D1B, alB, alB) - cB * cB) / (varB * varB);
    }
}
