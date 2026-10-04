/**
 * NLP Solver Utility for Linearly-Constrained Nonlinear Optimization.
 *
 * @since LINE 3.0
 */
package jline.api.mapqn;

import com.quantego.josqp.CSCMatrix;
import com.quantego.josqp.OSQP;

import java.util.ArrayList;
import java.util.List;

/**
 * Smooth minimisation over a polytope, by proximal projected gradient with an
 * OSQP projection oracle.
 *
 * <p>The QRF no-blocking family (MMI, MEM and their load-dependent variants)
 * minimises a smooth entropy-like objective over
 * {@code Aeq x = beq, Aub x <= bub, lb <= x <= ub}, with several hundred to
 * several thousand variables. Two earlier designs failed on it and are worth
 * recording:
 *
 * <ul>
 *   <li>A derivative-free inner solver (BOBYQA) needs O(n^2) interpolation
 *       points; on the M=3, N=4 instance (228 variables) it did not return
 *       within 21 CPU-minutes. An analytic gradient is therefore REQUIRED
 *       here, not optional.</li>
 *   <li>An augmented Lagrangian around a projected L-BFGS inner solver drove
 *       the penalty to its 1e8 ceiling with the constraint violation still
 *       oscillating around 1e-3: at that penalty the subproblem is too
 *       ill-conditioned for a limited-memory method to solve accurately, and
 *       an inaccurate subproblem solution makes the multiplier update
 *       meaningless.</li>
 * </ul>
 *
 * <p>The scheme below keeps every iterate ON the polytope instead of
 * penalising departures from it. Each step solves
 * {@code argmin_y  g'(y - x) + ||y - x||^2 / (2 t)} over the SAME polytope --
 * a strictly convex QP, which is where OSQP is at its most reliable -- and
 * then line-searches on the segment {@code [x, y]}, which is feasible by
 * convexity. The QP matrix {@code P = I / t} is constant, so only the linear
 * cost changes between iterations and no refactorisation is needed.
 *
 * <p>The objective is smooth but not convex (the MMI form carries
 * {@code -p log p'} cross terms), so this converges to a stationary point,
 * exactly as the MATLAB twin's {@code fmincon} and the native-Python twin's
 * {@code scipy.optimize.minimize(method='SLSQP')} do.
 */
public final class Mapqn_nlp_solver {
    private Mapqn_nlp_solver() {}

    /** Objective value at x. */
    public interface ObjectiveFn {
        double apply(double[] x);
    }

    /** Objective gradient at x, accumulated into gradOut (pre-zeroed by the caller). */
    public interface GradientFn {
        void apply(double[] x, double[] gradOut);
    }

    /** Armijo sufficient-decrease coefficient. */
    private static final double ARMIJO_C1 = 1.0e-4;

    /** Maximum backtracking steps per line search. */
    private static final int MAX_BACKTRACK = 30;

    /** Stationarity tolerance on the directional derivative of the prox step. */
    private static final double STATIONARY_TOL = 1.0e-10;

    /** Initial proximal step length. */
    private static final double STEP_INIT = 1.0;

    /** Proximal step below which no further progress is possible. */
    private static final double STEP_MIN = 1.0e-12;

    /** ADMM tolerance of the projection subproblem. */
    private static final double OSQP_EPS = 1.0e-9;

    public static double[] solve(ObjectiveFn objective, GradientFn gradient, int nVars,
                                 double[][] Aeq, double[] beq,
                                 double[][] Aub, double[] bub,
                                 double[] lb, double[] ub, double[] x0) {
        return solve(objective, gradient, nVars, Aeq, beq, Aub, bub, lb, ub, x0, 300);
    }

    /**
     * Minimise {@code objective} subject to {@code Aeq x = beq},
     * {@code Aub x <= bub} and {@code lb <= x <= ub}, by sequential quadratic
     * programming in the null space of the equalities.
     *
     * @param maxIter maximum SQP steps
     */
    public static double[] solve(ObjectiveFn objective, GradientFn gradient, int nVars,
                                 double[][] Aeq, double[] beq,
                                 double[][] Aub, double[] bub,
                                 double[] lb, double[] ub, double[] x0, int maxIter) {
        double[] x = x0.clone();
        for (int i = 0; i < nVars; i++) x[i] = Math.min(ub[i], Math.max(lb[i], x[i]));

        double[][] Z = nullSpace(Aeq, nVars);
        int d = (Z.length == 0) ? 0 : Z[0].length;
        if (d == 0) {
            // The equalities pin a single point; x IS the answer, not a failure
            // to optimise.
            return x;
        }

        // Reduced inequality system A t <= b, holding x0 + Z t inside the box
        // and inside the original inequality block. A box row whose Z row
        // vanishes is a coordinate the equalities PIN: in t it reads 0 <= 0,
        // which constrains nothing and only degrades the active-set solve.
        List<double[]> aRows = new ArrayList<double[]>();
        List<Double> bVals = new ArrayList<Double>();
        double zMax = 0.0;
        for (int i = 0; i < nVars; i++) {
            double nrm = 0.0;
            for (int k = 0; k < d; k++) nrm += Z[i][k] * Z[i][k];
            zMax = Math.max(zMax, Math.sqrt(nrm));
        }
        double zTol = 1.0e-12 * Math.max(1.0, zMax);
        for (int i = 0; i < nVars; i++) {
            double nrm = 0.0;
            for (int k = 0; k < d; k++) nrm += Z[i][k] * Z[i][k];
            if (Math.sqrt(nrm) <= zTol) continue;
            double[] up = new double[d];
            double[] lo = new double[d];
            for (int k = 0; k < d; k++) {
                up[k] = Z[i][k];
                lo[k] = -Z[i][k];
            }
            aRows.add(up);
            bVals.add(Double.valueOf(ub[i] - x[i]));
            aRows.add(lo);
            bVals.add(Double.valueOf(x[i] - lb[i]));
        }
        if (Aub != null && bub != null) {
            for (int r = 0; r < Aub.length; r++) {
                double[] row = new double[d];
                double nrm = 0.0;
                for (int k = 0; k < d; k++) {
                    double v = 0.0;
                    for (int i = 0; i < nVars; i++) {
                        if (Aub[r][i] != 0.0) v += Aub[r][i] * Z[i][k];
                    }
                    row[k] = v;
                    nrm += v * v;
                }
                double slack = bub[r];
                for (int i = 0; i < nVars; i++) {
                    if (Aub[r][i] != 0.0) slack -= Aub[r][i] * x[i];
                }
                if (Math.sqrt(nrm) <= 1.0e-12) continue;
                aRows.add(row);
                bVals.add(Double.valueOf(slack));
            }
        }
        int mRed = aRows.size();
        double[][] A = aRows.toArray(new double[0][]);
        double[] b = new double[mRed];
        for (int r = 0; r < mRed; r++) b[r] = bVals.get(r).doubleValue();

        double[] t = new double[d];
        double[] B = identity(d);
        double[] xt = new double[nVars];
        double[] gFull = new double[nVars];
        double[] gr = reducedGradient(gradient, objective, Z, x, t, nVars, d, gFull, xt);
        double f = objective.apply(expand(x, Z, t, nVars, d, xt));

        for (int iter = 0; iter < maxIter; iter++) {
            // SQP step: min 0.5 s'Bs + gr's subject to A (t + s) <= b. Every
            // iterate therefore stays feasible, and so does the whole segment
            // [t, t + s] by convexity, which is what lets the line search below
            // move on the objective alone.
            double[] rhs = new double[mRed];
            for (int r = 0; r < mRed; r++) {
                double at = 0.0;
                for (int k = 0; k < d; k++) at += A[r][k] * t[k];
                rhs[r] = b[r] - at;
            }
            double[] step = solveQp(B, gr, A, rhs, d, mRed);
            if (step == null) break;

            double slope = 0.0;
            for (int k = 0; k < d; k++) slope += gr[k] * step[k];
            if (slope > -STATIONARY_TOL) break;

            double alpha = 1.0;
            boolean accepted = false;
            double fTrial = f;
            double[] tTrial = new double[d];
            for (int bt = 0; bt < MAX_BACKTRACK; bt++) {
                for (int k = 0; k < d; k++) tTrial[k] = t[k] + alpha * step[k];
                fTrial = objective.apply(expand(x, Z, tTrial, nVars, d, xt));
                if (fTrial <= f + ARMIJO_C1 * alpha * slope) {
                    accepted = true;
                    break;
                }
                alpha *= 0.5;
            }
            if (!accepted) break;

            double[] grNew = reducedGradient(gradient, objective, Z, x, tTrial, nVars, d,
                    gFull, xt);
            double[] sVec = new double[d];
            double[] yVec = new double[d];
            for (int k = 0; k < d; k++) {
                sVec[k] = tTrial[k] - t[k];
                yVec[k] = grNew[k] - gr[k];
            }
            bfgsUpdate(B, sVec, yVec, d);

            System.arraycopy(tTrial, 0, t, 0, d);
            gr = grNew;
            f = fTrial;
            if (System.getProperty("jline.mapqn.debug") != null) {
                System.err.printf("[mapqn_nlp] iter=%d f=%.9f slope=%.3e alpha=%.3g%n",
                        iter, f, slope, alpha);
            }
        }
        return expand(x, Z, t, nVars, d, new double[nVars]).clone();
    }

    /** x0 + Z t, clipped to [0,1], written into out. */
    private static double[] expand(double[] x0, double[][] Z, double[] t, int nVars, int d,
                                   double[] out) {
        for (int i = 0; i < nVars; i++) {
            double v = x0[i];
            for (int k = 0; k < d; k++) v += Z[i][k] * t[k];
            out[i] = Math.min(1.0, Math.max(0.0, v));
        }
        return out;
    }

    /**
     * Z' grad f at x0 + Z t. The objective is p log p, so it is evaluated on
     * the CLIPPED point and a clipped coordinate carries no sensitivity.
     */
    private static double[] reducedGradient(GradientFn gradient, ObjectiveFn objective,
                                            double[][] Z, double[] x0, double[] t,
                                            int nVars, int d, double[] gFull, double[] xt) {
        double[] raw = new double[nVars];
        for (int i = 0; i < nVars; i++) {
            double v = x0[i];
            for (int k = 0; k < d; k++) v += Z[i][k] * t[k];
            raw[i] = v;
            xt[i] = Math.min(1.0, Math.max(0.0, v));
        }
        java.util.Arrays.fill(gFull, 0.0);
        gradient.apply(xt, gFull);
        for (int i = 0; i < nVars; i++) {
            if (raw[i] < 0.0 || raw[i] > 1.0) gFull[i] = 0.0;
        }
        double[] gr = new double[d];
        for (int k = 0; k < d; k++) {
            double v = 0.0;
            for (int i = 0; i < nVars; i++) {
                if (gFull[i] != 0.0) v += Z[i][k] * gFull[i];
            }
            gr[k] = v;
        }
        return gr;
    }

    /**
     * Orthonormal basis of the null space of {@code Aeq}, as an
     * {@code nVars x d} array, from the trailing right singular vectors.
     */
    private static double[][] nullSpace(double[][] Aeq, int nVars) {
        if (Aeq == null || Aeq.length == 0) {
            double[][] eye = new double[nVars][nVars];
            for (int i = 0; i < nVars; i++) eye[i][i] = 1.0;
            return eye;
        }
        org.apache.commons.math3.linear.RealMatrix M =
                new org.apache.commons.math3.linear.Array2DRowRealMatrix(Aeq, false);
        org.apache.commons.math3.linear.SingularValueDecomposition svd =
                new org.apache.commons.math3.linear.SingularValueDecomposition(M);
        double[] sv = svd.getSingularValues();
        double smax = (sv.length > 0) ? sv[0] : 0.0;
        double tol = Math.max(Aeq.length, nVars) * 2.220446049250313e-16 * smax;
        int rank = 0;
        for (int i = 0; i < sv.length; i++) if (sv[i] > tol) rank++;
        org.apache.commons.math3.linear.RealMatrix V = svd.getV();
        int d = nVars - rank;
        if (d <= 0) return new double[nVars][0];
        double[][] Z = new double[nVars][d];
        for (int i = 0; i < nVars; i++) {
            for (int k = 0; k < d; k++) {
                Z[i][k] = V.getEntry(i, rank + k);
            }
        }
        return Z;
    }

    /** min 0.5 s'Bs + g's subject to A s <= rhs, by OSQP. */
    private static double[] solveQp(double[] B, double[] g, double[][] A, double[] rhs,
                                    int d, int mRed) {
        int pnnz = d * (d + 1) / 2;
        int[] Pp = new int[d + 1];
        int[] Pi = new int[pnnz];
        double[] Px = new double[pnnz];
        int w = 0;
        for (int c = 0; c < d; c++) {
            Pp[c] = w;
            for (int r = 0; r <= c; r++) {
                Pi[w] = r;
                Px[w] = B[r * d + c];
                w++;
            }
        }
        Pp[d] = w;

        int[] colCount = new int[d];
        for (int r = 0; r < mRed; r++) {
            for (int c = 0; c < d; c++) if (A[r][c] != 0.0) colCount[c]++;
        }
        int nnz = 0;
        for (int c = 0; c < d; c++) nnz += colCount[c];
        int[] Ap = new int[d + 1];
        int pos = 0;
        for (int c = 0; c < d; c++) {
            Ap[c] = pos;
            pos += colCount[c];
        }
        Ap[d] = pos;
        int[] Ai = new int[nnz];
        double[] Ax = new double[nnz];
        int[] fill = new int[d];
        double[] lo = new double[mRed];
        double[] hi = new double[mRed];
        for (int r = 0; r < mRed; r++) {
            lo[r] = -OSQP.OSQP_INFTY;
            hi[r] = rhs[r];
            for (int c = 0; c < d; c++) {
                if (A[r][c] == 0.0) continue;
                int q = Ap[c] + fill[c]++;
                Ai[q] = r;
                Ax[q] = A[r][c];
            }
        }
        try {
            OSQP.Settings settings = new OSQP.Settings();
            settings.max_iter = 200000;
            settings.eps_abs = OSQP_EPS;
            settings.eps_rel = OSQP_EPS;
            // see Mapqn_qr_bounds_bas: josqp's adaptive rho stalls on these models
            settings.adaptive_rho = false;
            settings.polish = true;
            settings.verbose = false;
            OSQP solver = new OSQP(new OSQP.Data(d, mRed,
                    new CSCMatrix(d, d, pnnz, Pp, Pi, Px),
                    new CSCMatrix(mRed, d, nnz, Ap, Ai, Ax),
                    g.clone(), lo, hi, 0.0), settings);
            OSQP.Status status = solver.solve();
            if (status != OSQP.Status.SOLVED && status != OSQP.Status.SOLVED_INACCURATE) {
                return null;
            }
            return extractPrimal(solver, d);
        } catch (Exception e) {
            if (System.getProperty("jline.mapqn.debug") != null) {
                System.err.println("[mapqn_nlp] QP subproblem failed: " + e);
            }
            return null;
        }
    }

    private static double[] identity(int d) {
        double[] B = new double[d * d];
        for (int i = 0; i < d; i++) B[i * d + i] = 1.0;
        return B;
    }

    /** Powell-damped BFGS update, so B stays positive definite. */
    private static void bfgsUpdate(double[] B, double[] s, double[] y, int d) {
        double[] Bs = new double[d];
        for (int i = 0; i < d; i++) {
            double v = 0.0;
            for (int j = 0; j < d; j++) v += B[i * d + j] * s[j];
            Bs[i] = v;
        }
        double sBs = 0.0;
        double sy = 0.0;
        for (int i = 0; i < d; i++) {
            sBs += s[i] * Bs[i];
            sy += s[i] * y[i];
        }
        if (sBs <= 0.0) return;
        if (sy < 0.2 * sBs) {
            double theta = 0.8 * sBs / (sBs - sy);
            for (int i = 0; i < d; i++) y[i] = theta * y[i] + (1.0 - theta) * Bs[i];
            sy = 0.0;
            for (int i = 0; i < d; i++) sy += s[i] * y[i];
        }
        if (sy <= 1.0e-14) return;
        for (int i = 0; i < d; i++) {
            for (int j = 0; j < d; j++) {
                B[i * d + j] += y[i] * y[j] / sy - Bs[i] * Bs[j] / sBs;
            }
        }
    }

    /**
     * A point of the polytope {@code Aeq x = beq}, {@code Aub x <= bub},
     * {@code lb <= x <= ub} nearest to {@code x0}. Returns null when the
     * polytope solve fails.
     */
    public static double[] feasibleStart(double[][] Aeq, double[] beq,
                                         double[][] Aub, double[] bub, int nVars) {
        double[] lb = new double[nVars];
        double[] ub = new double[nVars];
        for (int i = 0; i < nVars; i++) ub[i] = 1.0;
        Polytope poly = new Polytope(Aeq, beq, Aub, bub, lb, ub, nVars);
        return new Projector(poly, 1.0).project(new double[nVars]);
    }

    /** Sparse CSC form of the polytope's constraint matrix, built once. */
    private static final class Polytope {
        final int nVars;
        final int mRows;
        final int nnz;
        final int[] Ap;
        final int[] Ai;
        final double[] Ax;
        final double[] rowLower;
        final double[] rowUpper;

        Polytope(double[][] Aeq, double[] beq, double[][] Aub, double[] bub,
                 double[] lb, double[] ub, int nVars) {
            int numEq = (Aeq != null && beq != null) ? Aeq.length : 0;
            int numIneq = (Aub != null && bub != null) ? Aub.length : 0;
            this.nVars = nVars;
            this.mRows = numEq + numIneq + nVars;

            int[] colCount = new int[nVars];
            for (int i = 0; i < numEq; i++) {
                for (int j = 0; j < nVars; j++) if (Aeq[i][j] != 0.0) colCount[j]++;
            }
            for (int i = 0; i < numIneq; i++) {
                for (int j = 0; j < nVars; j++) if (Aub[i][j] != 0.0) colCount[j]++;
            }
            int count = nVars;
            for (int j = 0; j < nVars; j++) count += colCount[j];
            this.nnz = count;

            this.Ap = new int[nVars + 1];
            int pos = 0;
            for (int j = 0; j < nVars; j++) {
                Ap[j] = pos;
                pos += colCount[j] + 1;
            }
            Ap[nVars] = pos;

            this.Ai = new int[nnz];
            this.Ax = new double[nnz];
            this.rowLower = new double[mRows];
            this.rowUpper = new double[mRows];
            int[] fill = new int[nVars];

            int r = 0;
            for (int i = 0; i < numEq; i++, r++) {
                rowLower[r] = beq[i];
                rowUpper[r] = beq[i];
                for (int j = 0; j < nVars; j++) {
                    if (Aeq[i][j] == 0.0) continue;
                    int w = Ap[j] + fill[j]++;
                    Ai[w] = r;
                    Ax[w] = Aeq[i][j];
                }
            }
            for (int i = 0; i < numIneq; i++, r++) {
                rowLower[r] = -OSQP.OSQP_INFTY;
                rowUpper[r] = bub[i];
                for (int j = 0; j < nVars; j++) {
                    if (Aub[i][j] == 0.0) continue;
                    int w = Ap[j] + fill[j]++;
                    Ai[w] = r;
                    Ax[w] = Aub[i][j];
                }
            }
            for (int j = 0; j < nVars; j++, r++) {
                int w = Ap[j] + fill[j]++;
                Ai[w] = r;
                Ax[w] = 1.0;
                rowLower[r] = lb[j];
                rowUpper[r] = ub[j];
            }
        }
    }

    /**
     * Proximal projection onto a {@link Polytope} at a fixed step length: the
     * QP {@code min 0.5 x'Px + q'x} with {@code P = I / t}. P is fixed, so the
     * OSQP instance is built once per step length and only its linear cost is
     * updated between iterations.
     */
    private static final class Projector {
        private final Polytope poly;
        private final double t;
        private OSQP solver;

        Projector(Polytope poly, double t) {
            this.poly = poly;
            this.t = t;
        }

        private void build() {
            int n = poly.nVars;
            int[] Pp = new int[n + 1];
            int[] Pi = new int[n];
            double[] Px = new double[n];
            for (int j = 0; j < n; j++) {
                Pp[j] = j;
                Pi[j] = j;
                Px[j] = 1.0 / t;
            }
            Pp[n] = n;
            CSCMatrix P = new CSCMatrix(n, n, n, Pp, Pi, Px);
            CSCMatrix A = new CSCMatrix(poly.mRows, n, poly.nnz, poly.Ap, poly.Ai, poly.Ax);
            OSQP.Settings settings = new OSQP.Settings();
            settings.max_iter = 200000;
            settings.eps_abs = OSQP_EPS;
            settings.eps_rel = OSQP_EPS;
            // see Mapqn_qr_bounds_bas: josqp's adaptive rho stalls on these models
            settings.adaptive_rho = false;
            settings.polish = true;
            settings.verbose = false;
            solver = new OSQP(new OSQP.Data(n, poly.mRows, P, A, new double[n],
                    poly.rowLower, poly.rowUpper, 0.0), settings);
        }

        /** Point of the polytope nearest to x0 (prox of the zero function). */
        double[] project(double[] x0) {
            double[] q = new double[poly.nVars];
            for (int i = 0; i < poly.nVars; i++) q[i] = -x0[i] / t;
            return solveWithCost(q, x0);
        }

        double[] solveWithCost(double[] q, double[] warmStart) {
            try {
                if (solver == null) build();
                solver.update_lin_cost(q);
                if (warmStart != null) solver.warm_start_x(warmStart);
                OSQP.Status status = solver.solve();
                if (status != OSQP.Status.SOLVED && status != OSQP.Status.SOLVED_INACCURATE) {
                    return null;
                }
                return extractPrimal(solver, poly.nVars);
            } catch (Exception e) {
                if (System.getProperty("jline.mapqn.debug") != null) {
                    System.err.println("[mapqn_nlp] projection failed: " + e);
                }
                return null;
            }
        }
    }

    /**
     * Read the primal solution out of a solved josqp workspace.
     *
     * <p>josqp 0.6.5 exposes {@code Workspace.solution} and {@code Solution.x}
     * as package-private, with no accessor, so reflection is the only route to
     * the primal point from outside {@code com.quantego.josqp}.
     */
    private static double[] extractPrimal(OSQP solver, int nVars) throws Exception {
        Object workspace = solver.getWorkspace();
        java.lang.reflect.Field solutionField = workspace.getClass().getDeclaredField("solution");
        solutionField.setAccessible(true);
        Object solution = solutionField.get(workspace);
        java.lang.reflect.Field xField = solution.getClass().getDeclaredField("x");
        xField.setAccessible(true);
        double[] x = (double[]) xField.get(solution);
        double[] copy = new double[nVars];
        System.arraycopy(x, 0, copy, 0, nVars);
        return copy;
    }
}
