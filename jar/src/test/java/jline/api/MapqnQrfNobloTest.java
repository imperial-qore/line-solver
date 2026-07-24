package jline.api;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.Timeout;
import static org.junit.jupiter.api.Assertions.*;

import java.util.concurrent.TimeUnit;

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
}
