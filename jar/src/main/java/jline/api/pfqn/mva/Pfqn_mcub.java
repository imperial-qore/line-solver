/**
 * @file Multiclass Composite Upper Bound (Kerola 1986)
 *
 * Kerola's composite bound method (Perf. Eval. 6:1-9, eqs. 10-16): a multiclass BJB lower
 * bound (eq. 10) seeds a residual-utilization composite upper bound (eqs. 13-16). Named
 * mcub because pfqn_cub is the unrelated cubature NC method. Ported at parity from MATLAB
 * pfqn_mcub.m.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.mva;

import jline.util.matrix.Matrix;

public final class Pfqn_mcub {
    private Pfqn_mcub() {}

    /**
     * Multiclass composite upper bound and its BJB lower seed.
     *
     * @param L service demand matrix, station x class (M x R)
     * @param N population vector (1 x R)
     * @param Z think time vector (1 x R)
     * @return {Xub, Xlb} 1xR upper (composite) and lower (multiclass BJB) throughput bounds
     */
    public static Matrix[] pfqn_mcub(Matrix L, Matrix N, Matrix Z) {
        int M = L.getNumRows();
        int R = L.getNumCols();
        double Ntot = N.elementSum();
        Matrix Xlb = new Matrix(1, R);
        Matrix Xub = new Matrix(1, R);
        // eq (10): multiclass Balanced Job Bounds lower throughput bound.
        for (int r = 0; r < R; r++) {
            double R0 = 0.0, Lb = 0.0;
            for (int k = 0; k < M; k++) {
                R0 += L.get(k, r);
                if (L.get(k, r) > Lb) Lb = L.get(k, r);
            }
            double zr = (Z != null && Z.getNumCols() > r) ? Z.get(0, r) : 0.0;
            Xlb.set(0, r, N.get(0, r) / (R0 + zr + (Ntot - 1) * Lb));
        }
        // eqs (13)-(16): composite upper bound per class.
        for (int r = 0; r < R; r++) {
            double best = Double.POSITIVE_INFINITY;
            for (int k = 0; k < M; k++) {
                if (L.get(k, r) > 0) {
                    double uoth = 0.0;
                    for (int s = 0; s < R; s++) {
                        if (s != r) uoth += Xlb.get(0, s) * L.get(k, s);
                    }
                    double dev = (1.0 - uoth) / L.get(k, r);
                    if (dev < best) best = dev;
                }
            }
            Xub.set(0, r, best);
        }
        return new Matrix[]{Xub, Xlb};
    }
}
