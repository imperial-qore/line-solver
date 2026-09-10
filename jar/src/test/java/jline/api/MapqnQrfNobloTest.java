package jline.api;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.Timeout;
import static org.junit.jupiter.api.Assertions.*;

import java.util.concurrent.TimeUnit;

import jline.api.mapqn.Mapqn_nlp_solver;
import jline.api.mapqn.Mapqn_qrf_noblo_mmi;
import jline.api.mapqn.Mapqn_qrf_noblo_mmi.ConstraintSet;

/**
 * Constraint-set tests for the QRF no-blocking polytope.
 *
 * Ground truth is the AMPL model shipped with the QRF paper
 * (qrboundsbas_skel.mod, specialised to no blocking; the load-dependent form is
 * qrboundsrsrd_skel.mod), and an exact oracle that needs no solver at all.
 *
 * <p>The exact oracle: with M == 2 and N1 + N2 == N the pairwise joint p2 is
 * fully determined by the marginal P(N1 = n1), so the stationary distribution
 * of the underlying CTMC gives a point the polytope MUST contain. That checks
 * the sign, the q lookup and every index of THM30 and THM3, which is where the
 * marginal-balance families go wrong silently.
 *
 * <p>Regression for two defects in buildConstraints:
 * <ul>
 *   <li>THM30 and THM3 were not emitted at all. They and THM1 are the only
 *       blocks that mention q, and THM1 alone balances phases, so it is
 *       identically zero at K == 1: the polytope then carried no dependence on
 *       the service rates and the U bound was the vacuous [0, 1] on every
 *       instance (glpsol on the AMPL model gives 0.6666666667 for the first
 *       case below, not [0, 1]).
 *   <li>The load-dependent variants built a population-free 4D q and dropped
 *       alpha on the floor, so a load-dependent instance was indistinguishable
 *       from a constant-rate one.
 * </ul>
 */
public class MapqnQrfNobloTest {

    private static final int MR = 1;

    /** Two-station cyclic network with exponential service at both stations. */
    private static double[][][][] expQ(double mu1, double mu2) {
        int M = 2;
        double[][] rt = {{0.0, 1.0}, {1.0, 0.0}};
        double[] mu = {mu1, mu2};
        double[][][][] q = new double[M][M][1][1];
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < M; j++) {
                q[i][j][0][0] = (j != i) ? rt[i][j] * mu[i] : 0.0;
            }
        }
        return q;
    }

    private static double[][][][][] expQLd(double mu1, double mu2, int N, double[][] alpha) {
        int M = 2;
        double[][] rt = {{0.0, 1.0}, {1.0, 0.0}};
        double[] mu = {mu1, mu2};
        double[][][][][] q = new double[M][M][1][1][N + 1];
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < M; j++) {
                for (int n = 1; n <= N; n++) {
                    double a = alpha[i][n - 1];
                    q[i][j][0][0][n] = (j != i) ? rt[i][j] * mu[i] * a : 0.0;
                }
            }
        }
        return q;
    }

    /**
     * The exact pairwise joint of the 2-station cyclic M/M/1//N network, laid
     * out in the decision vector.
     */
    private static double[] exactPoint(double mu1, double mu2, int N) {
        int M = 2;
        int[] K = {1, 1};
        int Kmax = 1;
        double rho = mu2 / mu1;
        double[] p = new double[N + 1];
        double z = 0.0;
        for (int n1 = 0; n1 <= N; n1++) { p[n1] = Math.pow(rho, n1); z += p[n1]; }
        for (int n1 = 0; n1 <= N; n1++) p[n1] /= z;

        int numVars = Mapqn_qrf_noblo_mmi.computeNumVars(M, N, Kmax, MR);
        double[] x = new double[numVars];
        for (int n1 = 0; n1 <= N; n1++) {
            int n2 = N - n1;
            x[Mapqn_qrf_noblo_mmi.p2Index(0, n1, 0, 0, n1, 0, 0, M, N, Kmax, MR)] = p[n1];
            x[Mapqn_qrf_noblo_mmi.p2Index(1, n2, 0, 1, n2, 0, 0, M, N, Kmax, MR)] = p[n1];
            x[Mapqn_qrf_noblo_mmi.p2Index(0, n1, 0, 1, n2, 0, 0, M, N, Kmax, MR)] = p[n1];
            x[Mapqn_qrf_noblo_mmi.p2Index(1, n2, 0, 0, n1, 0, 0, M, N, Kmax, MR)] = p[n1];
        }
        x[Mapqn_qrf_noblo_mmi.eIndex(0, 0, M, N, Kmax, MR)] = 1.0 - p[0];
        x[Mapqn_qrf_noblo_mmi.eIndex(1, 0, M, N, Kmax, MR)] = 1.0 - p[N];
        return x;
    }

    private static ConstraintSet build(double[][][][] q, int N) {
        int M = 2;
        int[] K = {1, 1};
        int[] F = {N, N};
        double[][] BB = new double[1][M];
        int numVars = Mapqn_qrf_noblo_mmi.computeNumVars(M, N, 1, MR);
        return Mapqn_qrf_noblo_mmi.buildConstraints(q, M, MR, BB, F, N, K, 1, numVars);
    }

    private static double maxEqResidual(ConstraintSet cs, double[] x) {
        double worst = 0.0;
        for (int r = 0; r < cs.aeq.size(); r++) {
            double[] row = cs.aeq.get(r);
            double val = -cs.beq.get(r);
            for (int c = 0; c < row.length; c++) val += row[c] * x[c];
            worst = Math.max(worst, Math.abs(val));
        }
        return worst;
    }

    private static double maxIneqResidual(ConstraintSet cs, double[] x) {
        double worst = 0.0;
        for (int r = 0; r < cs.aub.size(); r++) {
            double[] row = cs.aub.get(r);
            double val = -cs.bub.get(r);
            for (int c = 0; c < row.length; c++) val += row[c] * x[c];
            worst = Math.max(worst, val);
        }
        return worst;
    }

    @Test
    public void testExactJointIsFeasible() {
        double[][] cases = {{1.0, 1.0, 2}, {4.0, 1.0, 3}, {0.5, 1.0, 3}, {1.0, 1.0, 4}};
        for (double[] c : cases) {
            double mu1 = c[0], mu2 = c[1];
            int N = (int) c[2];
            ConstraintSet cs = build(expQ(mu1, mu2), N);
            double[] x = exactPoint(mu1, mu2, N);
            assertEquals(0.0, maxEqResidual(cs, x), 1e-12,
                    "equalities violated at the exact joint, mu1=" + mu1 + " N=" + N);
            assertTrue(maxIneqResidual(cs, x) < 1e-12,
                    "inequalities violated at the exact joint, mu1=" + mu1 + " N=" + N);
        }
    }

    @Test
    public void testMarginalBalanceFamiliesAreEmitted() {
        int N = 3;
        ConstraintSet cs = build(expQ(1.0, 1.0), N);
        // THM30 contributes M*sum(K) = 2 rows and THM3 contributes sum(F) = 6.
        int thm30 = 2;
        int thm3 = 2 * N;
        assertEquals(133, cs.aeq.size(),
                "equality row count changed; THM30 (" + thm30 + " rows) and THM3 ("
                        + thm3 + " rows) are the families this pins");
    }

    @Test
    public void testPolytopeDependsOnServiceRates() {
        int N = 3;
        ConstraintSet a = build(expQ(1.0, 1.0), N);
        ConstraintSet b = build(expQ(4.0, 1.0), N);
        assertEquals(a.aeq.size(), b.aeq.size());
        boolean differs = false;
        for (int r = 0; r < a.aeq.size() && !differs; r++) {
            double[] ra = a.aeq.get(r);
            double[] rb = b.aeq.get(r);
            for (int c = 0; c < ra.length; c++) {
                if (Math.abs(ra[c] - rb[c]) > 1e-12) { differs = true; break; }
            }
        }
        assertTrue(differs, "the equality block does not depend on q at all: "
                + "THM30/THM3 are missing and the U bound degenerates to [0,1]");
    }

    @Test
    public void testLoadDependentPolytopeDependsOnAlpha() {
        int M = 2, N = 3;
        int[] K = {1, 1};
        int[] F = {N, N};
        double[][] BB = new double[1][M];
        int numVars = Mapqn_qrf_noblo_mmi.computeNumVars(M, N, 1, MR);
        double[][] one = {{1.0, 1.0, 1.0}, {1.0, 1.0, 1.0}};
        double[][] ramp = {{1.0, 2.0, 3.0}, {1.0, 1.0, 1.0}};

        ConstraintSet a = Mapqn_qrf_noblo_mmi.buildConstraintsLd(
                expQLd(1.0, 1.0, N, one), M, MR, BB, F, N, K, 1, numVars);
        ConstraintSet b = Mapqn_qrf_noblo_mmi.buildConstraintsLd(
                expQLd(1.0, 1.0, N, ramp), M, MR, BB, F, N, K, 1, numVars);

        assertEquals(a.aeq.size(), b.aeq.size());
        boolean differs = false;
        for (int r = 0; r < a.aeq.size() && !differs; r++) {
            double[] ra = a.aeq.get(r);
            double[] rb = b.aeq.get(r);
            for (int c = 0; c < ra.length; c++) {
                if (Math.abs(ra[c] - rb[c]) > 1e-12) { differs = true; break; }
            }
        }
        assertTrue(differs, "the load-dependent block ignores alpha, so a "
                + "load-dependent instance is indistinguishable from a constant one");
    }

    @Test
    public void testConstantAlphaReproducesThePopulationFreePolytope() {
        int M = 2, N = 3;
        int[] K = {1, 1};
        int[] F = {N, N};
        double[][] BB = new double[1][M];
        int numVars = Mapqn_qrf_noblo_mmi.computeNumVars(M, N, 1, MR);
        double[][] one = {{1.0, 1.0, 1.0}, {1.0, 1.0, 1.0}};

        ConstraintSet ld = Mapqn_qrf_noblo_mmi.buildConstraintsLd(
                expQLd(1.0, 1.0, N, one), M, MR, BB, F, N, K, 1, numVars);
        double[] x = exactPoint(1.0, 1.0, N);
        assertEquals(0.0, maxEqResidual(ld, x), 1e-12);
        assertTrue(maxIneqResidual(ld, x) < 1e-12);
    }

    @Test
    @Timeout(value = 5, unit = TimeUnit.MINUTES)
    public void testMmiEndToEndExponential() {
        int M = 2, N = 2;
        int[] K = {1, 1};
        double[][][] mu = {{{1.0}}, {{1.0}}};
        double[][][] v = {{{0.0}}, {{0.0}}};
        double[][] rt = {{0.0, 1.0}, {1.0, 0.0}};

        jline.api.mapqn.Mapqn_solution sol =
                Mapqn_qrf_noblo_mmi.solve(M, MR, K, N, mu, v, rt);
        double u1 = sol.getVariables().get("UN_1");
        double q1 = sol.getVariables().get("QN_1");
        double q2 = sol.getVariables().get("QN_2");
        // glpsol pins U1 to 2/3 exactly on this instance (min == max).
        assertEquals(2.0 / 3.0, u1, 1e-4);
        assertEquals(N, q1 + q2, 1e-3);
    }

    // ---- qrf.bethe: the tree-reweighted arm ------------------------------

    /** Three-station cyclic exponential chain, the smallest instance with slack. */
    private static double[][][][] expQ3(double[] mu) {
        int M = 3;
        double[][][][] q = new double[M][M][1][1];
        for (int i = 0; i < M; i++) {
            q[i][(i + 1) % M][0][0] = mu[i];
        }
        return q;
    }

    /**
     * betheGradient is the gradient of betheObjective.
     *
     * <p>Central-differenced coordinate by coordinate in the FULL variable
     * space, which is legitimate because both routines are defined there: the
     * polytope only decides which directions the solver may take, not what the
     * derivative is. The step is small against LOGTOL, so a coordinate sitting
     * at zero is differenced through log(LOGTOL) rather than through a
     * singularity -- which is exactly the sensitivity the restored n = 0 cells
     * introduce.
     */
    @Test
    public void testBetheGradientMatchesFiniteDifferences() {
        int M = 2, N = 3, Kmax = 1;
        int[] K = {1, 1};
        int[] F = {N, N};
        double[] x = exactPoint(1.5, 1.0, N);
        int numVars = Mapqn_qrf_noblo_mmi.computeNumVars(M, N, Kmax, MR);
        double[] g = new double[numVars];
        Mapqn_qrf_noblo_mmi.betheGradient(x, g, M, N, K, Kmax, F, MR);

        double h = 1e-9;
        for (int c = 0; c < numVars; c++) {
            double[] xp = x.clone();
            double[] xm = x.clone();
            xp[c] += h;
            xm[c] -= h;
            double num = (Mapqn_qrf_noblo_mmi.betheObjective(xp, M, N, K, Kmax, F, MR)
                    - Mapqn_qrf_noblo_mmi.betheObjective(xm, M, N, K, Kmax, F, MR)) / (2 * h);
            assertEquals(num, g[c], 1e-4 * Math.max(1.0, Math.abs(num)),
                    "coordinate " + c);
        }
    }

    /**
     * The population sums of the Bethe objective start at n = 0, which is what
     * the coded MMI and MEM bodies do not do.
     *
     * <p>Its difference against lambda*MMI + MEM must be exactly the n = 0
     * block, term for term. This is the assertion that would fail if the loop
     * bounds were copied from mmiObjective, and the block is a real number
     * rather than a structural zero: p_ii at n = 0 is P(n_i = 0) = 1 - U_i.
     */
    @Test
    public void testBetheSumsTheIdleCells() {
        // A SYNTHETIC point, not a feasible one: this is an identity between
        // three objective FUNCTIONS, and the point only has to make the n = 0
        // block visible.
        //
        // The three bodies do NOT share a lower bound, and each matches the
        // AMPL line it transcribes. The MI line reads `sum {ni, nj in 0..F}`,
        // so mmiObjective was corrected to n = 0 on 2026-09-02 (defect D1).
        // The MEM line reads `sum {ni in 1..F[i]}`, so memObjective keeps its
        // n = 1 bound -- it already matched its own spec. betheObjective sums
        // both terms from n = 0 by construction. The residual between
        // lambda*MI + MEM and bethe is therefore the MEM n = 0 block ALONE;
        // the MI n = 0 block is now inside mmiObjective.
        int M = 2, N = 2, Kmax = 2;
        int[] K = {2, 1};
        int[] F = {N, N};
        final double LOGTOL = 1e-6;
        final double lambda = 1.0 / M;
        int numVars = Mapqn_qrf_noblo_mmi.computeNumVars(M, N, Kmax, MR);
        double[] x = new double[numVars];
        for (int i = 0; i < numVars; i++) {
            x[i] = 0.05 + 0.4 * ((0.37 * i) % 1.0);
        }

        double idle = 0.0;
        for (int i = 0; i < M; i++) {
            for (int k = 0; k < K[i]; k++) {
                double pv = x[Mapqn_qrf_noblo_mmi.p2Index(i, 0, k, i, 0, k, 0, M, N, Kmax, MR)];
                idle += pv * Math.log(LOGTOL + pv);
            }
        }

        double coded = lambda * Mapqn_qrf_noblo_mmi.mmiObjective(x, M, N, K, Kmax, F, MR)
                + Mapqn_qrf_noblo_mmi.memObjective(x, M, N, K, Kmax, F, MR);
        assertTrue(Math.abs(idle) > 1e-3, "the n = 0 block is numerically empty");
        assertEquals(coded + idle,
                Mapqn_qrf_noblo_mmi.betheObjective(x, M, N, K, Kmax, F, MR), 1e-10);
    }

    /**
     * Where the polytope pins a single point every objective must report it,
     * and qrf.bethe is no exception: glpsol gives U1 = 2/3 with min == max on
     * this instance.
     */
    @Test
    @Timeout(value = 10, unit = TimeUnit.MINUTES)
    public void testBetheEndToEndExponential() {
        int M = 2, N = 2;
        int[] K = {1, 1};
        double[][][] mu = {{{1.0}}, {{1.0}}};
        double[][][] v = {{{0.0}}, {{0.0}}};
        double[][] rt = {{0.0, 1.0}, {1.0, 0.0}};

        jline.api.mapqn.Mapqn_solution sol =
                jline.api.mapqn.Mapqn_qrf_noblo_bethe.solve(M, MR, K, N, mu, v, rt);
        double u1 = sol.getVariables().get("UN_1");
        double q1 = sol.getVariables().get("QN_1");
        double q2 = sol.getVariables().get("QN_2");
        assertEquals(2.0 / 3.0, u1, 1e-4);
        assertEquals(N, q1 + q2, 1e-3);
    }

    /**
     * The property convexity buys: the answer is a property of the model, not
     * of where the descent started.
     *
     * <p>The Bethe objective at lambda = 1/M is the negative of a
     * tree-reweighted entropy whose edge weight rho_ij = 2/M is the uniform
     * point of the spanning tree polytope of K_M, so the program is convex and
     * every local optimum is global. The second start is the MEM optimum of
     * the SAME polytope -- a genuinely different feasible point (U1 differs in
     * the third digit here), not a perturbation of the first.
     *
     * <p>Three stations at N = 3 is the smallest instance where this can fail:
     * at M == 2, and at M == 3 with N == 2, the equality rows pin a single
     * point and no objective can move the answer.
     */
    @Test
    @Timeout(value = 20, unit = TimeUnit.MINUTES)
    public void testBetheOptimumDoesNotDependOnTheStartPoint() {
        final int M = 3, N = 3, Kmax = 1;
        final int[] K = {1, 1, 1};
        final int[] F = {N, N, N};
        final double[] rates = {1.0, 1.5, 2.0};
        double[][] BB = new double[1][M];
        final int numVars = Mapqn_qrf_noblo_mmi.computeNumVars(M, N, Kmax, MR);

        ConstraintSet cs = Mapqn_qrf_noblo_mmi.buildConstraints(
                expQ3(rates), M, MR, BB, F, N, K, Kmax, numVars);
        double[][] Aeq = cs.aeq.toArray(new double[0][]);
        double[] beq = new double[cs.beq.size()];
        for (int i = 0; i < beq.length; i++) beq[i] = cs.beq.get(i);
        double[][] Aub = cs.aub.isEmpty() ? null : cs.aub.toArray(new double[0][]);
        double[] bub = null;
        if (!cs.bub.isEmpty()) {
            bub = new double[cs.bub.size()];
            for (int i = 0; i < bub.length; i++) bub[i] = cs.bub.get(i);
        }
        double[] lb = new double[numVars];
        double[] ub = new double[numVars];
        for (int i = 0; i < numVars; i++) ub[i] = 1.0;

        double[] phase1 = Mapqn_nlp_solver.feasibleStart(Aeq, beq, Aub, bub, numVars);
        assertNotNull(phase1, "phase 1 found no feasible point");

        Mapqn_nlp_solver.ObjectiveFn memObj = new Mapqn_nlp_solver.ObjectiveFn() {
            public double apply(double[] x) {
                return Mapqn_qrf_noblo_mmi.memObjective(x, M, N, K, Kmax, F, MR);
            }
        };
        Mapqn_nlp_solver.GradientFn memGrad = new Mapqn_nlp_solver.GradientFn() {
            public void apply(double[] x, double[] gradOut) {
                Mapqn_qrf_noblo_mmi.memGradient(x, gradOut, M, N, K, Kmax, F, MR);
            }
        };
        double[] memOpt = Mapqn_nlp_solver.solve(memObj, memGrad, numVars,
                Aeq, beq, Aub, bub, lb, ub, phase1);

        Mapqn_nlp_solver.ObjectiveFn obj = new Mapqn_nlp_solver.ObjectiveFn() {
            public double apply(double[] x) {
                return Mapqn_qrf_noblo_mmi.betheObjective(x, M, N, K, Kmax, F, MR);
            }
        };
        Mapqn_nlp_solver.GradientFn grad = new Mapqn_nlp_solver.GradientFn() {
            public void apply(double[] x, double[] gradOut) {
                Mapqn_qrf_noblo_mmi.betheGradient(x, gradOut, M, N, K, Kmax, F, MR);
            }
        };
        double[] fromPhase1 = Mapqn_nlp_solver.solve(obj, grad, numVars,
                Aeq, beq, Aub, bub, lb, ub, phase1);
        double[] fromMem = Mapqn_nlp_solver.solve(obj, grad, numVars,
                Aeq, beq, Aub, bub, lb, ub, memOpt);

        // The two starts are genuinely different points of the same polytope.
        double startGap = 0.0;
        for (int i = 0; i < numVars; i++) {
            startGap = Math.max(startGap, Math.abs(phase1[i] - memOpt[i]));
        }
        assertTrue(startGap > 1e-4,
                "the MEM optimum coincides with the phase-1 point, so this asserts nothing");

        jline.api.mapqn.Mapqn_solution a =
                Mapqn_qrf_noblo_mmi.extractResults(fromPhase1, M, N, K, Kmax, F, MR);
        jline.api.mapqn.Mapqn_solution b =
                Mapqn_qrf_noblo_mmi.extractResults(fromMem, M, N, K, Kmax, F, MR);
        for (int i = 1; i <= M; i++) {
            // 1e-5 rather than solver precision: the restored n = 0 cells make
            // the DESCENT DIRECTION sensitive to LOGTOL even though the
            // objective value is not, so the two runs stop at slightly
            // different points of the same flat optimum.
            assertEquals(a.getVariables().get("UN_" + i),
                    b.getVariables().get("UN_" + i), 1e-5,
                    "station " + i + ": the optimum moved with the start point");
        }
        double f1 = Mapqn_qrf_noblo_mmi.betheObjective(fromPhase1, M, N, K, Kmax, F, MR);
        double f2 = Mapqn_qrf_noblo_mmi.betheObjective(fromMem, M, N, K, Kmax, F, MR);
        assertEquals(f1, f2, 1e-8, "the optimal value moved with the start point");

        // ... and the point it reports is feasible.
        assertTrue(maxEqResidual(cs, fromPhase1) < 1e-7);
        assertTrue(maxIneqResidual(cs, fromPhase1) < 1e-7);
    }
}
