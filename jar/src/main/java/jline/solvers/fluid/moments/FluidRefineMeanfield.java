package jline.solvers.fluid.moments;

import jline.io.Ret;
import jline.util.matrix.Matrix;
import org.apache.commons.math3.util.FastMath;

import static jline.io.InputOutput.line_error;
import static jline.io.InputOutput.mfilename;

/**
 * Refined mean field correction of a fluid fixed point (Gast, POMACS 2017).
 *
 * <p>The mean-field fixed point x* is the leading term of an expansion of the true
 * stationary mean in powers of the system size. The next term is obtained by
 * carrying the second moment through the drift: writing A for the Jacobian at x*
 * and B for its Hessian tensor, the correction V solves the linear system
 *
 * <pre>
 *   A*V + (1/2) * sum_{j,k} Sigma_{jk} * d2F/dx_j dx_k = 0
 * </pre>
 *
 * with Sigma the stationary covariance from {@link FluidLyapunov}. Because Sigma
 * scales with the population, V is the O(1/N) term of the expansion written
 * directly in job counts, so no explicit density rescaling is needed. The Hessian
 * contraction is evaluated WITHOUT EVER FORMING THE TENSOR: writing
 * Sigma = sum_m lam_m*v_m*v_m' by eigendecomposition,
 *
 * <pre>
 *   sum_{jk} Sigma_{jk} d2F/dx_j dx_k = sum_m lam_m * d2F/dv_m^2
 * </pre>
 *
 * and each directional second derivative is one central second difference, so the
 * cost is O(rank(Sigma)) drift evaluations rather than O(n^2).
 *
 * <p>The drift must be TWICE DIFFERENTIABLE for this to mean anything. The
 * first-order closure is only piecewise linear -- its second derivative is zero
 * away from the kink and a delta at it -- so this must be called on the
 * Gaussian-closed drift, i.e. with the sigma2 that the moment closure converged to
 * under {@code options.method='refined'}. Passing sigma2 = 0 is rejected rather
 * than silently returning zero.
 *
 * <p>Port of {@code matlab/src/solvers/FLD/fluid_refine_meanfield.m}.
 *
 * @see FluidLyapunov
 * @see FluidMomentTerms
 */
public class FluidRefineMeanfield {

    /** Relative step of the second difference. */
    private static final double EPS_REL = 1e-4;

    /** The correction, plus what it took to compute it. */
    public static final class Result {
        /** The (n x 1) correction to be added to x. */
        public final double[] V;
        /** Number of covariance directions kept. */
        public final int rank;
        /** Step used by the second difference. */
        public final double stepsize;
        /** Residual norm of A*V + b. */
        public final double residual;
        /** Condition number of the Jacobian on the reachable subspace. */
        public final double condition;

        Result(double[] V, int rank, double stepsize, double residual, double condition) {
            this.V = V;
            this.rank = rank;
            this.stepsize = stepsize;
            this.residual = residual;
            this.condition = condition;
        }
    }

    private FluidRefineMeanfield() {
    }

    /**
     * The O(1/N) correction of a fluid fixed point.
     *
     * @param x      fluid fixed point
     * @param sigma2 per-station population variances defining the smooth drift
     * @param sigma  (n x n) stationary covariance
     * @param terms  representation from {@link FluidMomentTerms}
     * @param covblk per-station covariance blocks closing the DPS share ratio
     * @return the correction and its diagnostics
     */
    public static Result refine(double[] x, double[] sigma2, Matrix sigma,
                                FluidMomentTerms terms, Matrix[] covblk) {
        boolean anyPositive = false;
        for (double s : sigma2) {
            if (s > 0) {
                anyPositive = true;
                break;
            }
        }
        if (!anyPositive) {
            line_error(mfilename(new Object() {
            }), "The refined mean field expansion needs a twice-differentiable drift, but the "
                    + "first-order closure is only piecewise linear. Reach this function through "
                    + "options.method=\"refined\", which converges the Gaussian closure first.");
        }

        int n = x.length;

        // eigendecomposition of the covariance, dropping numerically null directions
        Matrix sym = new Matrix(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                sym.set(i, j, 0.5 * (sigma.get(i, j) + sigma.get(j, i)));
            }
        }
        Ret.Eigs eigs = sym.eigvec();
        double lamMax = 0;
        for (int m = 0; m < n; m++) {
            lamMax = FastMath.max(lamMax, eigs.values.get(m));
        }
        double eps = 2.220446049250313e-16;
        double keepAbove = lamMax * FastMath.sqrt(eps);

        double scale = FastMath.max(1, norm2(x));
        double step = EPS_REL * scale;

        Matrix F0 = terms.drift(x, sigma2, covblk);
        double[] b = new double[n];
        int rank = 0;
        for (int m = 0; m < n; m++) {
            double lam = eigs.values.get(m);
            if (!(lam > keepAbove) || !(lam > 0)) {
                continue;
            }
            rank++;
            double[] xp = new double[n];
            double[] xm = new double[n];
            for (int j = 0; j < n; j++) {
                double d = eigs.vectors.get(j, m);
                xp[j] = x[j] + step * d;
                xm[j] = x[j] - step * d;
            }
            Matrix Fp = terms.drift(xp, sigma2, covblk);
            Matrix Fm = terms.drift(xm, sigma2, covblk);
            for (int j = 0; j < n; j++) {
                b[j] += lam * (Fp.get(j) - 2 * F0.get(j) + Fm.get(j)) / (step * step);
            }
        }
        for (int j = 0; j < n; j++) {
            b[j] *= 0.5;
        }

        // solve A*V = -b on the REACHABLE subspace, where A is invertible
        Matrix A = terms.jacobian(x, sigma2, covblk);
        Matrix basis = FluidLyapunov.orth(terms.D);
        int k = basis.getNumCols();
        Matrix Ar = new Matrix(k, k);
        for (int a = 0; a < k; a++) {
            for (int c = 0; c < k; c++) {
                double acc = 0;
                for (int i = 0; i < n; i++) {
                    double vi = basis.get(i, a);
                    if (vi == 0) {
                        continue;
                    }
                    for (int j = 0; j < n; j++) {
                        acc += vi * A.get(i, j) * basis.get(j, c);
                    }
                }
                Ar.set(a, c, acc);
            }
        }
        double condAr = cond(Ar);
        if (!Double.isFinite(condAr) || condAr > 1.0 / FastMath.sqrt(eps)) {
            line_error(mfilename(new Object() {
            }), String.format("The fluid Jacobian is numerically singular on the reachable subspace "
                    + "(condition number %.3g), so the refinement equation A*V = -b has no meaningful "
                    + "solution. The fixed point sits at a drift kink or the model is marginally stable; "
                    + "use options.method=\"minnormal\", which resums the same correction without "
                    + "inverting A.", condAr));
        }
        Matrix br = new Matrix(k, 1);
        for (int a = 0; a < k; a++) {
            double acc = 0;
            for (int i = 0; i < n; i++) {
                acc += basis.get(i, a) * b[i];
            }
            br.set(a, 0, -acc);
        }
        Matrix Vr = new Matrix(k, 1);
        if (!Matrix.solveSafe(Ar, br, Vr)) {
            line_error(mfilename(new Object() {
            }), "The refinement equation A*V = -b could not be solved on the reachable subspace; "
                    + "use options.method=\"minnormal\", which resums the same correction without "
                    + "inverting A.");
        }
        double[] V = new double[n];
        for (int i = 0; i < n; i++) {
            double acc = 0;
            for (int a = 0; a < k; a++) {
                acc += basis.get(i, a) * Vr.get(a, 0);
            }
            V[i] = acc;
        }

        // The refinement is the NEXT TERM of an asymptotic expansion, so it is only
        // meaningful while it stays small against the leading term; a correction of
        // the same size as the fixed point means the expansion has not kicked in at
        // this population, and returning it would be worse than refusing.
        double normV = norm2(V);
        double normX = norm2(x);
        if (normV > 0.5 * FastMath.max(normX, FastMath.sqrt(eps))) {
            line_error(mfilename(new Object() {
            }), String.format("The 1/N refinement (norm %.3g) is not small against the mean-field fixed "
                    + "point (norm %.3g), so the asymptotic expansion is outside its range of validity at "
                    + "this population. Use options.method=\"minnormal\".", normV, normX));
        }

        double residual = 0;
        for (int i = 0; i < n; i++) {
            double acc = b[i];
            for (int j = 0; j < n; j++) {
                acc += A.get(i, j) * V[j];
            }
            residual += acc * acc;
        }
        return new Result(V, rank, step, FastMath.sqrt(residual), condAr);
    }

    private static double norm2(double[] v) {
        double acc = 0;
        for (double e : v) {
            acc += e * e;
        }
        return FastMath.sqrt(acc);
    }

    /** 2-norm condition number, from the singular values of a square matrix. */
    private static double cond(Matrix A) {
        int k = A.getNumRows();
        if (k == 0) {
            return 1.0;
        }
        org.apache.commons.math3.linear.RealMatrix rm =
                org.apache.commons.math3.linear.MatrixUtils.createRealMatrix(A.toArray2D());
        org.apache.commons.math3.linear.SingularValueDecomposition svd =
                new org.apache.commons.math3.linear.SingularValueDecomposition(rm);
        return svd.getConditionNumber();
    }
}
