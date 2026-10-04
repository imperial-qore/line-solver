/**
 * @file Limited Joint Dependence (LJD) lattice indexing helper for population-vector lookup tables.
 *
 * Maps a per-class population vector to a linearized index. This is the
 * mechanism shared by any state-dependent scaling that is supplied as a table
 * over a bounded population range (for example the throughput table produced by
 * the flow-equivalent-server aggregation, see FESAggregator).
 *
 * Mirrors the MATLAB helper ljd_linearize.m.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.ld;

import org.apache.commons.math3.util.FastMath;

import jline.util.matrix.Matrix;

public final class Ljd {
    private Ljd() {}

    /**
     * Convert per-class population vector to linearized index (0-indexed).
     *
     * idx = n1 + n2*(N1+1) + n3*(N1+1)*(N2+1) + ...
     *
     * @param nvec    per-class population vector
     * @param cutoffs per-class cutoffs [N1, N2, ..., NR]
     * @return the linearized index
     */
    public static int ljd_linearize(Matrix nvec, Matrix cutoffs) {
        int K = nvec.length();
        int idx = 0;
        int multiplier = 1;

        for (int k = 0; k < K; k++) {
            int nk = (int) FastMath.min(nvec.get(k), cutoffs.get(k));
            idx += nk * multiplier;
            multiplier *= ((int) cutoffs.get(k) + 1);
        }
        return idx;
    }

    /**
     * Inverse of {@link #ljd_linearize}: the population vector behind an index.
     *
     * <p>The forward map is a mixed-radix numeral with class k in radix
     * (Nk+1), so the inverse is the digit-by-digit division that reads it
     * back. It is what lets a caller walk a tabulated dependence in index
     * order and still know which population each entry belongs to.</p>
     *
     * @param idx 0-based linearized index
     * @param cutoffs per-class population cutoffs
     * @return the per-class population vector
     */
    public static Matrix ljd_delinearize(int idx, Matrix cutoffs) {
        int K = cutoffs.length();
        Matrix nvec = new Matrix(1, K);
        int rem = idx;
        for (int k = 0; k < K; k++) {
            int radix = (int) cutoffs.get(k) + 1;
            nvec.set(0, k, rem % radix);
            rem /= radix;
        }
        return nvec;
    }
}
