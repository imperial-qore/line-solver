/**
 * @file Norlund-Rice Edgeworth (NRE) method for load-dependent normalizing constants
 *
 * Port of matlab/src/api/pfqn/pfqn_nre.m. Evaluates the Norlund-Rice integral
 * form of the limited load-dependent normalizing constant by steepest descent
 * on a saddle-tilted contour, with a second-order Edgeworth correction.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.util.FastMath;

import jline.GlobalConstants;
import jline.api.pfqn.ld.Pfqn_gld;
import jline.api.pfqn.ld.Pfqn_lldsingle;
import jline.io.InputOutput;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;

/**
 * Normalizing constant via the saddle-tilted Edgeworth (NRE) approximation.
 *
 * Two corrections are applied over Pfqn_nrl and Pfqn_nrp: the integrand is
 * invariant under t -> t + c*1, so the redundant direction is quotiented out
 * and the integral is (R-1)-dimensional; and the contour radii are tilted per
 * class to the saddle point, making the origin a stationary point of the
 * phase. A second-order Edgeworth term built from the third and fourth
 * cumulants of the tilted distribution then gives a relative error of
 * O(1/sum(N)^2). Every integrand evaluation is at real positive demands, so
 * unlike Pfqn_nrl and Pfqn_nrp no complex arithmetic is involved.
 */
public final class Pfqn_nre {
    private Pfqn_nre() {}

    /** Finite-difference step; results are flat over [5e-3,5e-2]. */
    private static final double HSTEP = 2e-2;
    /** Stencil offsets stay within [-4,4] along each dimension. */
    private static final int SPAN = 9;
    /** Beyond 8 classes the fourth-cumulant tensor is no longer affordable. */
    private static final int MAX_DIM = 7;

    /**
     * MATLAB's {@code [lG,G,lGs,vsad]}: the constant, the saddlepoint term on
     * its own and the tilt the expansion was taken about.
     *
     * <p>{@code lG - lGs} is the Edgeworth correction, so a caller that wants
     * the plain saddlepoint estimate reads {@code lGs} rather than re-deriving
     * it. {@code vsad} is null on the shortcut arms that never solve a saddle
     * point, matching the reference's empty {@code []}.
     */
    public static final class Result {
        /** Logarithm of the normalizing constant, Edgeworth correction included. */
        public final double lG;
        /** The normalizing constant, {@code exp(lG)}. */
        public final double G;
        /** Logarithm of the saddlepoint term alone, i.e. lG with the Edgeworth factor omitted. */
        public final double lGs;
        /** The tilt actually used (1xd), null when no saddle point was solved. */
        public final Matrix vsad;

        Result(double lG, double lGs, Matrix vsad) {
            this.lG = lG;
            this.G = FastMath.exp(lG);
            this.lGs = lGs;
            this.vsad = vsad;
        }
    }

    /**
     * Logarithm of the normalizing constant of a limited load-dependent model.
     *
     * @param Lin     service demand matrix (MxR)
     * @param N       population vector (1xR)
     * @param Z       think time vector (1xR or DxR)
     * @param alphaIn load-dependent rate matrix (Mx sum(N))
     * @param options solver options
     * @return logarithm of the normalizing constant
     */
    public static double pfqn_nre(Matrix Lin, Matrix N, Matrix Z, Matrix alphaIn, SolverOptions options) {
        return pfqn_nre_full(Lin, N, Z, alphaIn, options, null).lG;
    }

    /**
     * The full form of the reference's four outputs, named alike in the native
     * python and C++ ports: the constant, the saddlepoint term alone and the tilt used.
     *
     * @param Lin     service demand matrix (MxR)
     * @param N       population vector (1xR)
     * @param Z       think time vector (1xR or DxR)
     * @param alphaIn load-dependent rate matrix (Mx sum(N))
     * @param options solver options
     * @param vfix    tilt to use instead of solving the saddle-point equation,
     *                null for the standard estimator. Supplying the tilt
     *                obtained at a nearby population makes numerator and
     *                denominator of a ratio share one expansion point, which is
     *                the Tierney-Kadane arrangement.
     * @return the constant, the saddlepoint term and the tilt
     */
    public static Result pfqn_nre_full(Matrix Lin, Matrix N, Matrix Z, Matrix alphaIn,
                                       SolverOptions options, Matrix vfix) {
        double Ntd = N.elementSum();
        if (Ntd < 0) {
            return new Result(GlobalConstants.NegInf, GlobalConstants.NegInf, null);
        }
        if (Ntd == 0.0) {
            return new Result(0.0, 0.0, null);
        }
        int Nt = (int) FastMath.round(Ntd);
        if (alphaIn.getNumCols() < Nt) {
            throw new RuntimeException("pfqn_nre: the load-dependent rate matrix must have at least sum(N) columns.");
        }
        // trim so that every rate used downstream is positive
        Matrix alpha = new Matrix(alphaIn.getNumRows(), Nt);
        for (int i = 0; i < alphaIn.getNumRows(); i++) {
            for (int j = 0; j < Nt; j++) {
                alpha.set(i, j, alphaIn.get(i, j));
            }
        }
        Matrix L = Lin.copy();
        if (Z != null && !Z.isEmpty() && Z.elementSum() > 0) {
            Matrix Zrow = Z.getNumRows() > 1 ? Z.sumCols() : Z;
            L = Matrix.concatRows(L, Zrow, null);
            Matrix delayRates = new Matrix(1, Nt);
            for (int j = 0; j < Nt; j++) {
                delayRates.set(0, j, j + 1.0); // the delay is an infinite server station
            }
            alpha = Matrix.concatRows(alpha, delayRates, null);
        }
        int M = L.getNumRows();
        int R = L.getNumCols();
        if (M == 1) {
            double lGone = Pfqn_gld.pfqn_gld(L, N, alpha, options).lG;
            return new Result(lGone, lGone, null);
        }

        // scale demands in [0,1] per class, the residual factor is exact by homogeneity
        double lGscale = 0.0;
        for (int r = 0; r < R; r++) {
            double colMax = 0.0;
            for (int i = 0; i < M; i++) {
                if (L.get(i, r) > colMax) {
                    colMax = L.get(i, r);
                }
            }
            if (colMax <= 0) {
                colMax = 1.0;
            }
            for (int i = 0; i < M; i++) {
                L.set(i, r, L.get(i, r) / colMax);
            }
            lGscale += N.get(r) * FastMath.log(colMax);
        }

        Matrix Ntmat = new Matrix(1, 1);
        Ntmat.set(0, 0, Nt);
        if (R == 1) {
            // coefficient extraction is the identity in a single class
            double lGone = Pfqn_lldsingle.pfqn_lldsingle(L, Ntmat, alpha, options).lG + lGscale;
            return new Result(lGone, lGone, null);
        }

        int d = R - 1; // dimension of the quotient torus
        if (d > MAX_DIM) {
            throw new RuntimeException("pfqn_nre: pfqn_nre is limited to 8 classes, use nrl or clw beyond that.");
        }
        double[] Nd = new double[d];
        for (int a = 0; a < d; a++) {
            Nd[a] = N.get(a);
        }
        Cgf cgf = new Cgf(L, Ntmat, alpha, options, d);

        // ---- saddle point: minimise the convex F(v) = K(v) - Nd*v ----
        // A tilt supplied by the caller is used as given, so that a ratio of two
        // constants can be expanded about one common point rather than two.
        double[] vbase = new double[d];
        boolean converged = false;
        boolean tiltGiven = vfix != null && !vfix.isEmpty();
        if (tiltGiven) {
            if (vfix.length() < d) {
                throw new RuntimeException("pfqn_nre: the supplied tilt must have one entry per "
                        + "quotient dimension (R-1).");
            }
            for (int a = 0; a < d; a++) {
                vbase[a] = vfix.get(a);
            }
            converged = true;
        }
        for (int it = 0; !tiltGiven && it < 100; it++) {
            cgf.reset(vbase);
            Matrix grad = new Matrix(d, 1);
            Matrix hess = new Matrix(d, d);
            for (int a = 0; a < d; a++) {
                grad.set(a, 0, (cgf.at(unitoff(d, a, 1)) - cgf.at(unitoff(d, a, -1))) / (2 * HSTEP) - Nd[a]);
            }
            for (int a = 0; a < d; a++) {
                for (int b = 0; b < d; b++) {
                    hess.set(a, b, secondDiff(cgf, d, a, b));
                }
            }
            Matrix step = hess.inv().mult(grad);
            double F0 = cgf.at(new int[d]) - dot(Nd, vbase);
            double tau = 1;
            double[] vtry = new double[d];
            while (tau > 1e-10) {
                for (int a = 0; a < d; a++) {
                    vtry[a] = vbase[a] - tau * step.get(a, 0);
                }
                if (cgf.atPoint(vtry) - dot(Nd, vtry) <= F0) {
                    break;
                }
                tau = tau / 2;
            }
            double stepNorm = 0.0;
            for (int a = 0; a < d; a++) {
                double delta = -tau * step.get(a, 0);
                vbase[a] += delta;
                stepNorm += delta * delta;
            }
            // Newton converges to the root of the differenced gradient, whose own
            // O(hstep^2) bias puts any absolute gradient target out of reach
            if (FastMath.sqrt(stepNorm) < 1e-10) {
                converged = true;
                break;
            }
        }
        if (!converged) {
            InputOutput.line_warning(InputOutput.mfilename(new Object()),
                    "the saddle point search did not converge, the estimate may be inaccurate.\n");
        }

        Matrix vsad = new Matrix(1, d);
        for (int a = 0; a < d; a++) {
            vsad.set(0, a, vbase[a]);
        }

        // ---- cumulants of the tilted distribution at the saddle ----
        cgf.reset(vbase);
        double K0 = cgf.at(new int[d]);
        Matrix Sigma = new Matrix(d, d);
        for (int a = 0; a < d; a++) {
            for (int b = 0; b < d; b++) {
                Sigma.set(a, b, secondDiff(cgf, d, a, b));
            }
        }
        for (int a = 0; a < d; a++) {
            for (int b = a + 1; b < d; b++) {
                double sym = (Sigma.get(a, b) + Sigma.get(b, a)) / 2;
                Sigma.set(a, b, sym);
                Sigma.set(b, a, sym);
            }
        }
        double minEig = Double.POSITIVE_INFINITY;
        List<Complex> eigs = Sigma.eig();
        for (int i = 0; i < eigs.size(); i++) {
            double re = eigs.get(i).getReal();
            if (re < minEig) {
                minEig = re;
            }
        }
        if (minEig <= 0) {
            throw new RuntimeException("pfqn_nre: the tilted covariance is singular, a class has no demand at any station.");
        }

        double[][][] k3 = new double[d][d][d];
        for (int a = 0; a < d; a++) {
            for (int b = 0; b < d; b++) {
                for (int c = 0; c < d; c++) {
                    double acc = 0;
                    for (int s = 0; s < 8; s++) {
                        int s1 = 1 - 2 * ((s >> 0) & 1);
                        int s2 = 1 - 2 * ((s >> 1) & 1);
                        int s3 = 1 - 2 * ((s >> 2) & 1);
                        int[] off = new int[d];
                        off[a] += s1;
                        off[b] += s2;
                        off[c] += s3;
                        acc += s1 * s2 * s3 * cgf.at(off);
                    }
                    k3[a][b][c] = acc / (8 * HSTEP * HSTEP * HSTEP);
                }
            }
        }

        double[][][][] k4 = new double[d][d][d][d];
        for (int a = 0; a < d; a++) {
            for (int b = 0; b < d; b++) {
                for (int c = 0; c < d; c++) {
                    for (int e = 0; e < d; e++) {
                        double acc = 0;
                        for (int s = 0; s < 16; s++) {
                            int s1 = 1 - 2 * ((s >> 0) & 1);
                            int s2 = 1 - 2 * ((s >> 1) & 1);
                            int s3 = 1 - 2 * ((s >> 2) & 1);
                            int s4 = 1 - 2 * ((s >> 3) & 1);
                            int[] off = new int[d];
                            off[a] += s1;
                            off[b] += s2;
                            off[c] += s3;
                            off[e] += s4;
                            acc += s1 * s2 * s3 * s4 * cgf.at(off);
                        }
                        k4[a][b][c][e] = acc / (16 * HSTEP * HSTEP * HSTEP * HSTEP);
                    }
                }
            }
        }

        // ---- second-order Edgeworth factor, see the class header ----
        Matrix Sinv = Sigma.inv();
        double[][] S = new double[d][d];
        for (int a = 0; a < d; a++) {
            for (int b = 0; b < d; b++) {
                S[a][b] = Sinv.get(a, b);
            }
        }
        double rho4 = 0;
        for (int a = 0; a < d; a++) {
            for (int b = 0; b < d; b++) {
                for (int c = 0; c < d; c++) {
                    for (int e = 0; e < d; e++) {
                        rho4 += k4[a][b][c][e] * S[a][b] * S[c][e];
                    }
                }
            }
        }
        double[] u = new double[d];
        for (int c = 0; c < d; c++) {
            double acc = 0;
            for (int a = 0; a < d; a++) {
                for (int b = 0; b < d; b++) {
                    acc += S[a][b] * k3[a][b][c];
                }
            }
            u[c] = acc;
        }
        double rhoA = 0;
        for (int c = 0; c < d; c++) {
            for (int e = 0; e < d; e++) {
                rhoA += u[c] * S[c][e] * u[e];
            }
        }
        // staged contraction of k3 against three copies of Sigma^{-1}
        double[][][] T1 = new double[d][d][d];
        for (int i = 0; i < d; i++) {
            for (int b = 0; b < d; b++) {
                for (int c = 0; c < d; c++) {
                    double acc = 0;
                    for (int a = 0; a < d; a++) {
                        acc += S[i][a] * k3[a][b][c];
                    }
                    T1[i][b][c] = acc;
                }
            }
        }
        double[][][] T2 = new double[d][d][d];
        for (int i = 0; i < d; i++) {
            for (int j = 0; j < d; j++) {
                for (int c = 0; c < d; c++) {
                    double acc = 0;
                    for (int b = 0; b < d; b++) {
                        acc += S[j][b] * T1[i][b][c];
                    }
                    T2[i][j][c] = acc;
                }
            }
        }
        double rhoB = 0;
        for (int i = 0; i < d; i++) {
            for (int j = 0; j < d; j++) {
                for (int k = 0; k < d; k++) {
                    double acc = 0;
                    for (int c = 0; c < d; c++) {
                        acc += S[k][c] * T2[i][j][c];
                    }
                    rhoB += k3[i][j][k] * acc;
                }
            }
        }
        double corr = 1 + rho4 / 8 - (3 * rhoA + 2 * rhoB) / 24;
        if (corr <= 0) {
            InputOutput.line_warning(InputOutput.mfilename(new Object()),
                    "the Edgeworth correction is non-positive, falling back on the saddlepoint term.\n");
            corr = 1;
        }

        double lGs = K0 - dot(Nd, vbase) - (d / 2.0) * FastMath.log(2 * FastMath.PI)
                - 0.5 * FastMath.log(Sigma.det()) + lGscale;
        return new Result(lGs + FastMath.log(corr), lGs, vsad);
    }

    /** Mixed second difference of the cumulant generating function on the stencil. */
    private static double secondDiff(Cgf cgf, int d, int a, int b) {
        int[] pp = new int[d];
        pp[a] += 1;
        pp[b] += 1;
        int[] pm = new int[d];
        pm[a] += 1;
        pm[b] -= 1;
        int[] mp = new int[d];
        mp[a] -= 1;
        mp[b] += 1;
        int[] mm = new int[d];
        mm[a] -= 1;
        mm[b] -= 1;
        return (cgf.at(pp) - cgf.at(pm) - cgf.at(mp) + cgf.at(mm)) / (4 * HSTEP * HSTEP);
    }

    private static int[] unitoff(int d, int a, int sgn) {
        int[] off = new int[d];
        off[a] = sgn;
        return off;
    }

    private static double dot(double[] x, double[] y) {
        double acc = 0;
        for (int i = 0; i < x.length; i++) {
            acc += x[i] * y[i];
        }
        return acc;
    }

    /**
     * Cumulant generating function of the tilted single-class model, memoised on
     * the finite-difference stencil around the current expansion point.
     */
    private static final class Cgf {
        private final Matrix L;
        private final Matrix Nt;
        private final Matrix alpha;
        private final SolverOptions options;
        private final int d;
        private final Map<Long, Double> cache = new HashMap<Long, Double>();
        private double[] vbase;

        Cgf(Matrix L, Matrix Nt, Matrix alpha, SolverOptions options, int d) {
            this.L = L;
            this.Nt = Nt;
            this.alpha = alpha;
            this.options = options;
            this.d = d;
            this.vbase = new double[d];
        }

        /** Move the expansion point, which invalidates every memoised value. */
        void reset(double[] v) {
            this.vbase = v.clone();
            this.cache.clear();
        }

        /** Value at vbase+off*hstep. */
        double at(int[] off) {
            long key = 0;
            for (int r = d - 1; r >= 0; r--) {
                key = key * SPAN + (off[r] + 4);
            }
            Double cached = cache.get(key);
            if (cached != null) {
                return cached;
            }
            double[] v = new double[d];
            for (int r = 0; r < d; r++) {
                v[r] = vbase[r] + off[r] * HSTEP;
            }
            double y = atPoint(v);
            cache.put(key, y);
            return y;
        }

        /** log of the single-class LLD constant at the class tilt X=[exp(v),1]. */
        double atPoint(double[] v) {
            int M = L.getNumRows();
            int R = L.getNumCols();
            Matrix Lx = new Matrix(M, 1);
            for (int i = 0; i < M; i++) {
                double acc = 0;
                for (int r = 0; r < R; r++) {
                    double xr = r < d ? FastMath.exp(v[r]) : 1.0;
                    acc += L.get(i, r) * xr;
                }
                Lx.set(i, 0, acc);
            }
            return Pfqn_lldsingle.pfqn_lldsingle(Lx, Nt, alpha, options).lG;
        }
    }
}
