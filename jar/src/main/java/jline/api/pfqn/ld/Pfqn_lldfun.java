/**
 * @file Limited load-dependent function evaluation with clamped linear interpolation
 *
 * Evaluates limited load-dependent (LLD) scaling functions at continuous queue-length values.
 * Supports multi-server stations with load-dependent service rates and interpolates linearly
 * between the discrete scaling points, clamping at the ends of the lattice.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.ld;

import org.apache.commons.math3.util.FastMath;

import jline.util.Maths;
import jline.util.Utils;
import jline.util.matrix.Matrix;

public final class Pfqn_lldfun {
    private Pfqn_lldfun() {}

    /**
     * Evaluate limited-load dependent (LLD) function.
     *
     * @param n          Queue-length values. The values can be continuous.
     * @param lldscaling If not null, then the LLD function uses lldscaling to interpolate continuous values of n
     * @param nservers   If not null, then the LLD function is set to be for a multi-server with nserver stations
     * @return Interpolated LLD function values
     */
    public static Matrix pfqn_lldfun(Matrix n, Matrix lldscaling, Matrix nservers) {
        int M = n.length();
        Matrix r = new Matrix(M, 1);
        r.fill(1.0);
        int smax = lldscaling.getNumCols();
        double alpha = 20.0;

        for (int i = 0; i < M; i++) {
            if (!(nservers == null || nservers.isEmpty())) {
                if (Utils.isInf(nservers.get(i))) {
                    r.set(i, 0, 1);
                } else {
                    double softminValue = r.get(i, 0) / Maths.softmin(n.get(i), nservers.get(i), alpha);
                    if (Double.isNaN(softminValue)) {
                        r.set(i, 0, 1.0 / FastMath.min(n.get(i), nservers.get(i)));
                    } else {
                        r.set(i, 0, softminValue);
                    }
                }
            }

            if (!lldscaling.isEmpty()) {
                Matrix lldscaling_i = new Matrix(1, smax);
                Matrix.extract(lldscaling, i, i + 1, 0, smax, lldscaling_i, 0, 0);
                if (lldscaling_i.elementMax() != lldscaling_i.elementMin()) {
                    r.set(i, 0, r.get(i, 0) / interpLatticeClamped(lldscaling, i, smax, n.get(i)));
                }
            }
        }
        return r;
    }

    /**
     * Clamped linear interpolation of row {@code i} of the lattice, sampled at n = 1..smax.
     *
     * <p>The lattice alpha(n) is defined only at integer n, and AMVA-QD evaluates it at the
     * interpolated mean queue length, so an interpolation rule is required. Linear is exact for
     * the piecewise-linear {@code min(1:N,c)} lattice that load dependence is overwhelmingly used
     * to express, whereas a cubic spline overshoots between the knots (at n=1.44 on min(n,2):
     * spline 1.6444 vs the exact 1.4400), granting more capacity than the model declares.
     *
     * <p>Clamping at n=smax mirrors the state-space convention (MATLAB State.afterEventStation.m
     * indexes {@code lldscaling(ist, min(ni, lldlimit))}): the last lattice entry is the saturated
     * rate. It also replaces the previous Apache SplineInterpolator, which THREW an
     * OutOfRangeException past the lattice ("3.018 out of [1, 3] range") for unbounded open
     * queues, where MATLAB's spline instead extrapolated to negative rates.
     */
    private static double interpLatticeClamped(Matrix lldscaling, int i, int smax, double n) {
        if (smax == 1) {
            return lldscaling.get(i, 0);
        }
        double x = n;
        if (!(x > 1.0)) { // also catches NaN -> clamp to the first knot
            x = 1.0;
        } else if (x > smax) {
            x = smax;
        }
        int lo = (int) FastMath.floor(x);
        if (lo > smax - 1) {
            lo = smax - 1; // keep lo, lo+1 inside the lattice at the right edge
        }
        double frac = x - lo;
        // lattice index j holds alpha(j+1), so knot lo is column lo-1
        return lldscaling.get(i, lo - 1) * (1.0 - frac) + lldscaling.get(i, lo) * frac;
    }
}
