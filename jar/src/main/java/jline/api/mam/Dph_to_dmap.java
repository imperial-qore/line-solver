package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * Renewal DMAP of a discrete phase-type law, and its inverse.
 *
 * <p>Maps (alpha, A) to {D0, D1} with D0 = A and D1 = a*alpha, a = e - A*e. The
 * process renews the phase at every event, so its interevent times are i.i.d.
 * copies of the DPH and D0+D1 is stochastic by construction.
 *
 * <p>MATLAB twins: dph_to_dmap.m, dmap_to_dph.m, dmap_is_renewal.m
 */
public final class Dph_to_dmap {
    private Dph_to_dmap() {}

    /** Builds the renewal DMAP {A, a*alpha} of the DPH (alpha, A). */
    public static MatrixCell dph_to_dmap(Matrix alpha, Matrix A) {
        int m = A.getNumRows();
        Matrix rowSums = A.mult(Matrix.ones(m, 1));
        Matrix D1 = new Matrix(m, m);
        for (int i = 0; i < m; i++) {
            double a = 1.0 - rowSums.get(i, 0);
            for (int j = 0; j < m; j++) {
                D1.set(i, j, a * alpha.get(0, j));
            }
        }
        return new MatrixCell(A.copy(), D1);
    }

    /**
     * True when the DMAP renews at every event, i.e. D1 has rank one. Tested by
     * rebuilding a*alpha from the row masses rather than by a rank routine, so
     * the tolerance is on the entries the caller will actually use.
     */
    public static boolean dmap_is_renewal(Matrix D0, Matrix D1) {
        int m = D1.getNumRows();
        if (m == 1) {
            return true;
        }
        Matrix rowMass = D1.mult(Matrix.ones(m, 1));
        int pivot = -1;
        double best = 0;
        for (int i = 0; i < m; i++) {
            if (rowMass.get(i, 0) > best) {
                best = rowMass.get(i, 0);
                pivot = i;
            }
        }
        if (pivot < 0 || best <= 1e-14) {
            return false;
        }
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                double expected = rowMass.get(i, 0) * D1.get(pivot, j) / best;
                if (Math.abs(D1.get(i, j) - expected) > 1e-8) {
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * Discrete phase-type law underlying a renewal DMAP: alpha is the
     * normalized pivot row of D1 and A is D0.
     */
    public static Matrix dmap_to_dph_alpha(Matrix D0, Matrix D1) {
        if (!dmap_is_renewal(D0, D1)) {
            throw new RuntimeException(
                    "The DMAP does not renew at events, so it has no discrete phase-type form.");
        }
        int m = D1.getNumRows();
        Matrix rowMass = D1.mult(Matrix.ones(m, 1));
        int pivot = 0;
        for (int i = 1; i < m; i++) {
            if (rowMass.get(i, 0) > rowMass.get(pivot, 0)) {
                pivot = i;
            }
        }
        if (rowMass.get(pivot, 0) <= 1e-14) {
            throw new RuntimeException("The DMAP has no events: D1 is the zero matrix.");
        }
        Matrix alpha = new Matrix(1, m);
        for (int j = 0; j < m; j++) {
            alpha.set(0, j, D1.get(pivot, j) / rowMass.get(pivot, 0));
        }
        return alpha;
    }
}
