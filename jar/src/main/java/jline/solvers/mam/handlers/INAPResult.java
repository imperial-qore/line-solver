/**
 * @file INAP fixed-point iteration result
 *
 * @since LINE 3.1.0
 */
package jline.solvers.mam.handlers;

import java.util.List;
import java.util.Objects;

import jline.util.matrix.Matrix;

/**
 * Result of the INAP/INAPplus fixed-point iteration.
 */
public final class INAPResult {
    public final Matrix x;
    public final List<Matrix> pi;
    public final List<Matrix> Q;
    public final int iter;
    // Per-process geometric-tail decay set by the matrix-geometric 'inapinf'
    // method; rhoProc/isGeomProc are null for the finite-state inap/inapplus.
    public final double[] rhoProc;
    public final boolean[] isGeomProc;
    public final double rcatRes;

    public INAPResult(Matrix x, List<Matrix> pi, List<Matrix> Q, int iter) {
        this(x, pi, Q, iter, null, null, 0.0);
    }

    public INAPResult(Matrix x, List<Matrix> pi, List<Matrix> Q, int iter,
                      double[] rhoProc, boolean[] isGeomProc, double rcatRes) {
        this.x = x;
        this.pi = pi;
        this.Q = Q;
        this.iter = iter;
        this.rhoProc = rhoProc;
        this.isGeomProc = isGeomProc;
        this.rcatRes = rcatRes;
    }

    public Matrix getX() { return x; }
    public List<Matrix> getPi() { return pi; }
    public List<Matrix> getQ() { return Q; }
    public int getIter() { return iter; }

    public Matrix component1() { return x; }
    public List<Matrix> component2() { return pi; }
    public List<Matrix> component3() { return Q; }
    public int component4() { return iter; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof INAPResult)) return false;
        INAPResult that = (INAPResult) o;
        return iter == that.iter
                && Objects.equals(x, that.x)
                && Objects.equals(pi, that.pi)
                && Objects.equals(Q, that.Q);
    }

    @Override
    public int hashCode() {
        return Objects.hash(x, pi, Q, iter);
    }

    @Override
    public String toString() {
        return "INAPResult(x=" + x + ", pi=" + pi + ", Q=" + Q + ", iter=" + iter + ")";
    }
}
