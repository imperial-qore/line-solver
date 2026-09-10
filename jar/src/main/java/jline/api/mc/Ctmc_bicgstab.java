/**
 * Preconditioned stabilized biconjugate gradients for sparse linear systems
 *
 * Solves the sparse nonsymmetric system A*x = b by the BiCGSTAB method of van der Vorst
 * (1992), right-preconditioned by the same threshold incomplete LU factorization ILUT of
 * Saad (1994) that Ctmc_gmres uses, with the same Jacobi fallback on breakdown. It is the
 * short-recurrence counterpart of that kernel: work and storage per iteration are
 * constant rather than growing with the Krylov dimension, so it does not restart and does
 * not lose the optimality that restarting costs GMRES. Where GMRES(m) stagnates because
 * the useful subspace is wider than m, this converges; where it does not, GMRES(m) is the
 * more robust of the two, hence the order in which Ctmc_solve tries them.
 *
 * The equilibration, reverse Cuthill-McKee reordering and preconditioner are shared with
 * Ctmc_gmres rather than reimplemented, so both methods factorize the same matrix in the
 * same order and a switch between them cannot move a reported metric for a reason other
 * than the iteration itself.
 *
 * @since LINE 3.0
 */
package jline.api.mc;

import jline.util.matrix.Matrix;

public final class Ctmc_bicgstab {
    private Ctmc_bicgstab() {}

    /**
     * Default linear-solve residual, as in Ctmc_gmres. Much tighter than the fixed-point
     * tolerance options.iter_tol: switching solve method must not move a reported metric.
     */
    public static final double BICGSTAB_DEFAULT_TOL = 1e-12;

    /**
     * Default cap on complete iterations. BiCGSTAB storage is O(n) regardless of the
     * count, so the cap bounds time rather than memory.
     */
    public static final int BICGSTAB_DEFAULT_MAXIT = 200;

    /**
     * Threshold below which the rho and omega scalars are treated as a breakdown, taken
     * relative to the norms whose product formed them.
     */
    private static final double BREAKDOWN_TOL = 1e-14;

    /**
     * Outcome of a BiCGSTAB solve. The flag follows the MATLAB bicgstab convention: 0
     * converged, 1 iteration limit reached, 2 preconditioner ill-conditioned, 3
     * stagnation, 4 a scalar quantity became too small or too large to continue. Callers
     * must check it and fall back to another solve when it is nonzero.
     *
     * iter counts matrix-vector products with A: two per complete iteration, and one
     * when the iteration converges at its half step, so an ODD count is normal. Counting
     * products rather than iterations is what makes it comparable with the iter of
     * Ctmc_gmres and across the four codebases.
     */
    public static final class BicgstabResult {
        public final Matrix x;
        public final int flag;
        public final double relres;
        public final int iter;

        BicgstabResult(Matrix x, int flag, double relres, int iter) {
            this.x = x;
            this.flag = flag;
            this.relres = relres;
            this.iter = iter;
        }
    }

    /**
     * Solve A*x = b by preconditioned BiCGSTAB with the default tolerance and iteration
     * cap, starting from a uniform initial guess.
     *
     * @param A Coefficient matrix
     * @param b Right-hand side, as a column vector
     * @return Solution, convergence flag, relative residual and matrix-vector product count
     */
    public static BicgstabResult ctmc_bicgstab(Matrix A, Matrix b) {
        return ctmc_bicgstab(A, b, BICGSTAB_DEFAULT_TOL, 0, null);
    }

    /**
     * Solve A*x = b by preconditioned BiCGSTAB.
     *
     * @param A     Coefficient matrix
     * @param b     Right-hand side, as a column vector
     * @param tol   Relative residual tolerance, nonpositive for the default 1e-12
     * @param maxit Maximum number of complete iterations, nonpositive for min(n,200)
     * @param x0    Initial guess, null for the uniform vector 1/n
     * @return Solution, convergence flag, relative residual and matvec count
     */
    public static BicgstabResult ctmc_bicgstab(Matrix A, Matrix b, double tol, int maxit, Matrix x0) {
        int n = A.getNumRows();
        if (b.length() != n) {
            throw new IllegalArgumentException("Matrix dimensions incompatible: A is " + n + "x" + n
                    + " but b has " + b.length() + " entries");
        }
        double[] rhs = new double[n];
        for (int i = 0; i < n; i++) rhs[i] = b.get(i);
        double[] guess = null;
        if (x0 != null) {
            guess = new double[n];
            for (int i = 0; i < n; i++) guess[i] = x0.get(i);
        }
        Ctmc_gmres.Prepared prep = new Ctmc_gmres.Prepared(A);
        return solvePrepared(prep, rhs, guess, tol, maxit);
    }

    /**
     * Solve A*X = B for every column of B, reusing one incomplete factorization across
     * all of them and starting each column from the previous solution. This is the shape
     * of the stochastic complement, whose right-hand side is a whole block of the
     * generator: refactorizing per column would cost more than the direct solve it
     * replaces.
     *
     * @param A     Coefficient matrix
     * @param B     Right-hand sides, one per column
     * @param tol   Relative residual tolerance, nonpositive for the default 1e-12
     * @param maxit Maximum number of complete iterations, nonpositive for min(n,200)
     * @return The solution block, or null if any column failed to converge. Returning
     *         null rather than a partial block keeps the caller's fallback all-or-nothing.
     */
    public static Matrix ctmc_bicgstab(Matrix A, Matrix B, double tol, int maxit) {
        int n = A.getNumRows();
        if (B.getNumRows() != n) {
            throw new IllegalArgumentException("Matrix dimensions incompatible: A is " + n + "x" + n
                    + " but B has " + B.getNumRows() + " rows");
        }
        Ctmc_gmres.Prepared prep = new Ctmc_gmres.Prepared(A);
        int nrhs = B.getNumCols();
        Matrix X = new Matrix(n, nrhs);
        double[] guess = null;
        for (int c = 0; c < nrhs; c++) {
            double[] rhs = new double[n];
            for (int i = 0; i < n; i++) rhs[i] = B.get(i, c);
            BicgstabResult r = solvePrepared(prep, rhs, guess, tol, maxit);
            if (r.flag != 0) return null;
            guess = new double[n];
            for (int i = 0; i < n; i++) {
                guess[i] = r.x.get(i, 0);
                if (guess[i] != 0.0) X.set(i, c, guess[i]);
            }
        }
        return X;
    }

    /**
     * Run the iteration on an already equilibrated, reordered and preconditioned system.
     * The right-hand side and the initial guess are given in the ORIGINAL ordering and
     * unscaled; both are mapped here, as Ctmc_gmres does, so that no caller has to know
     * about the permutation.
     */
    private static BicgstabResult solvePrepared(Ctmc_gmres.Prepared prep, double[] rhsIn,
                                                double[] x0In, double tol, int maxit) {
        int n = prep.n;
        double[] scaled = new double[n];
        for (int i = 0; i < n; i++) scaled[i] = rhsIn[i] / prep.rowScale[i];
        double[] start = new double[n];
        if (x0In == null) {
            java.util.Arrays.fill(start, 1.0 / n);
        } else {
            System.arraycopy(x0In, 0, start, 0, n);
        }
        return bicgstab(prep, Ctmc_gmres.permute(scaled, prep.perm),
                Ctmc_gmres.permute(start, prep.perm), tol, maxit);
    }

    /**
     * Right-preconditioned BiCGSTAB. The recurrence is van der Vorst's, with the two
     * preconditioner applications placed on the search directions p and s, so the
     * residual r is the residual of the ORIGINAL system and the convergence test needs no
     * unpreconditioning.
     */
    private static BicgstabResult bicgstab(Ctmc_gmres.Prepared prep, double[] rhs, double[] x,
                                           double tol, int maxit) {
        int n = prep.n;
        Ctmc_gmres.Csr csr = prep.csr;
        Ctmc_gmres.Precond M = prep.M;
        int[] perm = prep.perm;

        if (tol <= 0) tol = BICGSTAB_DEFAULT_TOL;
        if (maxit <= 0) maxit = Math.min(n, BICGSTAB_DEFAULT_MAXIT);
        maxit = Math.max(1, Math.min(maxit, n));

        double bnorm = Ctmc_gmres.norm2(rhs);
        if (bnorm == 0.0) bnorm = 1.0;

        double[] r = new double[n];
        csr.mult(x, r);
        for (int i = 0; i < n; i++) r[i] = rhs[i] - r[i];

        double relres = Ctmc_gmres.norm2(r) / bnorm;
        if (relres <= tol) {
            return new BicgstabResult(Ctmc_gmres.asColumn(Ctmc_gmres.unpermute(x, perm)), 0, relres, 0);
        }

        // The shadow residual is fixed at the initial residual, the standard choice: any
        // vector not orthogonal to r would do, and this one cannot be orthogonal to it.
        double[] rhat = r.clone();
        double[] p = new double[n];
        double[] v = new double[n];
        double[] s = new double[n];
        double[] t = new double[n];
        double[] ph = new double[n];
        double[] sh = new double[n];

        double rho = 1.0;
        double alpha = 1.0;
        double omega = 1.0;
        int matvec = 0;
        int flag = 1;
        double prevrelres = relres;

        for (int it = 0; it < maxit; it++) {
            double rhoNew = Ctmc_gmres.dot(rhat, r);
            // rho vanishing is the biorthogonality breakdown of the underlying Lanczos
            // process, not slow convergence: restarting with a fresh shadow vector would
            // discard the iterate, so the caller is told to use another method instead.
            if (Math.abs(rhoNew) <= BREAKDOWN_TOL * Ctmc_gmres.norm2(rhat) * Ctmc_gmres.norm2(r)) {
                flag = 4;
                break;
            }
            if (it == 0) {
                System.arraycopy(r, 0, p, 0, n);
            } else {
                if (omega == 0.0) {
                    flag = 4;
                    break;
                }
                double beta = (rhoNew / rho) * (alpha / omega);
                for (int i = 0; i < n; i++) p[i] = r[i] + beta * (p[i] - omega * v[i]);
            }
            rho = rhoNew;

            M.apply(p, ph);
            csr.mult(ph, v);
            matvec++;

            double rhatv = Ctmc_gmres.dot(rhat, v);
            if (rhatv == 0.0 || !isFinite(rhatv)) {
                flag = 4;
                break;
            }
            alpha = rho / rhatv;

            for (int i = 0; i < n; i++) s[i] = r[i] - alpha * v[i];

            // Half-step convergence: s is the residual of x + alpha*ph, so a converged s
            // means the answer is reached without the second matvec of this iteration.
            double snorm = Ctmc_gmres.norm2(s);
            if (snorm / bnorm <= tol) {
                for (int i = 0; i < n; i++) x[i] += alpha * ph[i];
                relres = snorm / bnorm;
                flag = 0;
                break;
            }

            M.apply(s, sh);
            csr.mult(sh, t);
            matvec++;

            double tt = Ctmc_gmres.dot(t, t);
            if (tt == 0.0 || !isFinite(tt)) {
                flag = 4;
                break;
            }
            omega = Ctmc_gmres.dot(t, s) / tt;

            for (int i = 0; i < n; i++) x[i] += alpha * ph[i] + omega * sh[i];
            for (int i = 0; i < n; i++) r[i] = s[i] - omega * t[i];

            relres = Ctmc_gmres.norm2(r) / bnorm;
            if (relres <= tol) {
                flag = 0;
                break;
            }
            // omega vanishing stalls the update of x while leaving r finite, so the
            // iteration would spin without progress.
            if (Math.abs(omega) <= BREAKDOWN_TOL) {
                flag = 4;
                break;
            }
            // BiCGSTAB residuals are non-monotone by construction, so an increase is not
            // by itself stagnation and the test is against the BEST residual seen rather
            // than the previous one. Growing two orders of magnitude past that best is
            // divergence, and continuing from such an iterate is not worth the matvecs.
            if (relres > 1e2 * prevrelres) {
                flag = 3;
                break;
            }
            prevrelres = Math.min(prevrelres, relres);
        }

        for (int i = 0; i < n; i++) {
            if (!isFinite(x[i])) {
                return new BicgstabResult(Ctmc_gmres.asColumn(Ctmc_gmres.unpermute(x, perm)), 4,
                        Double.POSITIVE_INFINITY, matvec);
            }
        }
        if (relres <= tol) flag = 0;
        return new BicgstabResult(Ctmc_gmres.asColumn(Ctmc_gmres.unpermute(x, perm)), flag, relres, matvec);
    }

    private static boolean isFinite(double v) {
        return !Double.isNaN(v) && !Double.isInfinite(v);
    }
}
