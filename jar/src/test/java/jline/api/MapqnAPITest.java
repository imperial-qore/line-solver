package jline.api;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Timeout;
import static org.junit.jupiter.api.Assertions.*;
import java.util.concurrent.TimeUnit;

import jline.api.mapqn.*;
import jline.util.matrix.Matrix;

import java.util.Arrays;

/**
 * Test class for MAPQN API functionality.
 * Based on examples from 08_mapqn_sigme08 repository.
 *
 * Reference values computed using SolverMVA from LINE MATLAB.
 * The LP bounds should bracket the exact MVA solution.
 *
 * Also includes comparison tests validating LINE's MAPQN bound implementations
 * against reference values from the original MATLAB implementations in:
 * - qrf-revised/qrf_rsrd.m (QR Bounds with RS-RD blocking)
 * - qrf-revised/qrf_bas.m (QR Bounds with BAS blocking)
 * - 08_mapqn_sigme08/bnd_linearreduction_new.mod (LR Bounds)
 */
public class MapqnAPITest {

    // Maximum allowed MAPE for bounds to be considered valid
    private static final double MAX_MAPE = 0.20;  // 20% max error for bounds

    // Maximum allowed gap between upper and lower bounds
    private static final double MAX_BOUND_GAP = 0.50;  // 50% max gap

    // Tolerance for comparing with MATLAB reference values
    private static final double MATLAB_TOLERANCE = 0.01;  // 1% tolerance

    // Stricter tolerance for exact results (symmetric networks)
    private static final double EXACT_TOLERANCE = 1e-6;

    /**
     * Test QR Bounds for BAS (Blocking After Service) network.
     * Configuration validated against MATLAB qrf_bas.m from qrf-revised repository.
     *
     * Configuration: M=2, N=2, K=[1,1], MR=1 (no blocking)
     * Reference from MATLAB: LB = UB = 0.666667 (exact for symmetric tandem)
     */
    @Test
    @Timeout(value = 2, unit = TimeUnit.MINUTES)
    public void testQR_BAS() {
        // Reference value from MATLAB qrf_bas.m
        final double EXPECTED_U = 2.0 / 3.0;  // 0.666667
        final double TOLERANCE = 0.01;  // 1% tolerance

        int M = 2;
        int N = 2;
        int[] K = {1, 1};  // 1 phase per queue
        int[] F = {2, 2};  // Both queues can hold N jobs

        // Service rates mu{i}(k,h) - exponential rate 1
        Matrix[] mu = new Matrix[M];
        mu[0] = new Matrix(new double[][]{{1.0}});
        mu[1] = new Matrix(new double[][]{{1.0}});

        // Background transition rates (zeros)
        Matrix[] v = new Matrix[M];
        v[0] = new Matrix(new double[][]{{0.0}});
        v[1] = new Matrix(new double[][]{{0.0}});

        // Routing probabilities - tandem
        Matrix r = new Matrix(new double[][]{
            {0.0, 1.0},
            {1.0, 0.0}
        });

        // Blocking configurations - minimal: just unblocked state
        int MR = 1;  // Only unblocked configuration
        int f = 1;   // Finite capacity queue (1-based)

        // BB(m, i) = 1 if queue i is blocked in configuration m
        Matrix BB = new Matrix(new double[][]{{0, 0}});  // m=0: no blocking

        // MM(m, j) = blocking order matrix (MR x M)
        Matrix MM = new Matrix(new double[][]{{0, 0}});  // No blocking order

        // MM1(m, j) = extended blocking order info
        Matrix MM1 = new Matrix(new double[][]{{0, 0}});

        // ZZ(m) = number of blocked queues in configuration m
        int[] ZZ = {0};

        Mapqn_qr_bounds_bas_parameters params = new Mapqn_qr_bounds_bas_parameters(
            M, N, MR, f, K, F, MM, MM1, ZZ, BB, mu, v, r);

        // Get lower and upper bounds for queue 1
        Mapqn_solution solutionMin = Mapqn_qr_bounds_bas.INSTANCE.solve(params, 1, "min");
        Mapqn_solution solutionMax = Mapqn_qr_bounds_bas.INSTANCE.solve(params, 1, "max");

        assertNotNull(solutionMin, "Min solution should not be null");
        assertNotNull(solutionMax, "Max solution should not be null");

        double ULB = solutionMin.getObjectiveValue();
        double UUB = solutionMax.getObjectiveValue();

        assertFalse(Double.isNaN(ULB), "Lower bound should not be NaN");
        assertFalse(Double.isNaN(UUB), "Upper bound should not be NaN");

        // Bounds should be in valid range [0, 1]
        assertTrue(ULB >= 0, String.format("Lower bound (%.4f) should be non-negative", ULB));
        assertTrue(UUB <= 1, String.format("Upper bound (%.4f) should not exceed 1", UUB));

        // Lower bound <= Upper bound
        assertTrue(ULB <= UUB + 1e-6,
                  String.format("Lower bound (%.4f) should not exceed upper bound (%.4f)", ULB, UUB));

        // Validate against MATLAB reference: LB and UB should both be close to 2/3
        assertEquals(EXPECTED_U, ULB, TOLERANCE,
                    String.format("Lower bound (%.6f) should match MATLAB reference (%.6f)", ULB, EXPECTED_U));
        assertEquals(EXPECTED_U, UUB, TOLERANCE,
                    String.format("Upper bound (%.6f) should match MATLAB reference (%.6f)", UUB, EXPECTED_U));
    }

    /**
     * Test QR Bounds for RSRD network.
     * Configuration validated against MATLAB qrf_rsrd.m from qrf-revised repository.
     *
     * Configuration: M=2, N=2, K=[1,1], F=[2,2] (symmetric tandem)
     * Reference from MATLAB: LB = UB = 0.666667 (exact for symmetric tandem)
     */
    @Test
    @Timeout(value = 2, unit = TimeUnit.MINUTES)
    public void testQR_RSRD() {
        // Reference value from MATLAB qrf_rsrd.m
        final double EXPECTED_U = 2.0 / 3.0;  // 0.666667
        final double TOLERANCE = 0.01;  // 1% tolerance

        int M = 2;
        int N = 2;
        int[] F = {2, 2};  // Both can hold all jobs
        int[] K = {1, 1};  // Single phase per queue

        // Service rates mu{i}(k,h) - exponential rate 1
        Matrix[] mu = new Matrix[M];
        mu[0] = new Matrix(new double[][]{{1.0}});
        mu[1] = new Matrix(new double[][]{{1.0}});

        // Background transition rates (zeros)
        Matrix[] v = new Matrix[M];
        v[0] = new Matrix(new double[][]{{0.0}});
        v[1] = new Matrix(new double[][]{{0.0}});

        // Routing probabilities - tandem network
        Matrix r = new Matrix(new double[][]{
            {0.0, 1.0},
            {1.0, 0.0}
        });

        // Load-dependent rates (default = 1 for all)
        double[][] alpha = new double[M][];
        for (int i = 0; i < M; i++) {
            alpha[i] = new double[N];
            Arrays.fill(alpha[i], 1.0);
        }

        Mapqn_qr_bounds_rsrd_parameters params = new Mapqn_qr_bounds_rsrd_parameters(M, N, F, K, mu, v, alpha, r);

        // Get lower and upper bounds for queue 1
        Mapqn_solution solutionMin = Mapqn_qr_bounds_rsrd.INSTANCE.solve(params, 1, "min");
        Mapqn_solution solutionMax = Mapqn_qr_bounds_rsrd.INSTANCE.solve(params, 1, "max");

        assertNotNull(solutionMin, "Min solution should not be null");
        assertNotNull(solutionMax, "Max solution should not be null");

        double ULB = solutionMin.getObjectiveValue();
        double UUB = solutionMax.getObjectiveValue();

        assertFalse(Double.isNaN(ULB), "Lower bound should not be NaN");
        assertFalse(Double.isNaN(UUB), "Upper bound should not be NaN");

        // Bounds should be in valid range [0, 1]
        assertTrue(ULB >= 0, String.format("Lower bound (%.4f) should be non-negative", ULB));
        assertTrue(UUB <= 1, String.format("Upper bound (%.4f) should not exceed 1", UUB));

        // Lower bound <= Upper bound
        assertTrue(ULB <= UUB + 1e-6,
                  String.format("Lower bound (%.4f) should not exceed upper bound (%.4f)", ULB, UUB));

        // Validate against MATLAB reference: LB and UB should both be close to 2/3
        assertEquals(EXPECTED_U, ULB, TOLERANCE,
                    String.format("Lower bound (%.6f) should match MATLAB reference (%.6f)", ULB, EXPECTED_U));
        assertEquals(EXPECTED_U, UUB, TOLERANCE,
                    String.format("Upper bound (%.6f) should match MATLAB reference (%.6f)", UUB, EXPECTED_U));
    }

    /**
     * Test BndMVAVersion model.
     */
    @Test
    @Timeout(value = 2, unit = TimeUnit.MINUTES)
    public void testBndMVAVersion() {
        int M = 3;
        int K = 2;
        int N = 10;

        double[] muM = {35.71, 25.0};

        Matrix muMAP = new Matrix(new double[][]{
            {0.2428699491756, 0.0248957749874},
            {0.3693899392817, 52.2199871959660}
        });

        Matrix r = new Matrix(new double[][]{
            {0.2, 0.7, 0.1},
            {1.0, 0.0, 0.0},
            {1.0, 0.0, 0.0}
        });

        Matrix v = new Matrix(new double[][]{
            {0.0, 0.0},
            {0.0, 0.0}
        });

        MVAVersionParameters params = new MVAVersionParameters(M, N, K, muM, muMAP, r, v);

        Mapqn_solution solution = Mapqn_bnd_lr_mva.INSTANCE.solve(params, 1, 1);
        assertNotNull(solution, "Solution should not be null");

        double util = solution.getObjectiveValue();
        assertTrue(util >= 0, "Total utilization should be non-negative");
        assertTrue(util <= 1, "Utilization should not exceed 1");

        // Test individual level utilizations
        for (int k = 1; k <= K; k++) {
            Mapqn_solution levelSolution = Mapqn_bnd_lr_mva.INSTANCE.solve(params, 1, k);
            assertNotNull(levelSolution, String.format("Solution for queue 1, level %d should exist", k));

            double levelUtil = levelSolution.getVariable("UN_1_" + k);
            assertTrue(levelUtil >= 0, String.format("Utilization for level %d should be non-negative", k));
            assertTrue(levelUtil <= 1.0, String.format("Utilization for level %d should not exceed 1", k));
        }

        // Verify routing conservation
        for (int i = 0; i < M; i++) {
            double rowSum = 0;
            for (int j = 0; j < M; j++) {
                rowSum += r.get(i, j);
            }
            assertEquals(1.0, rowSum, 1e-6,
                        String.format("Routing probabilities from queue %d should sum to 1", i + 1));
        }
    }

    /**
     * Test Mapqn_bnd_lr with symmetric 2-queue tandem network.
     *
     * Reference (MVA exact): U1 = U2 = 0.375 for N=3, mu=[1,1], tandem routing
     * LR bounds give upper bounds on utilization (when maximizing).
     */
    @Test
    @Timeout(value = 2, unit = TimeUnit.MINUTES)
    public void testMapqn_bnd_lr_symmetric_tandem() {
        // Reference values from MVA
        final double EXACT_U1 = 0.375;
        final double EXACT_U2 = 0.375;

        int M = 2;
        int N = 3;
        int[] K = {1, 1};

        Matrix[] mu = new Matrix[M];
        mu[0] = new Matrix(new double[][]{{1.0}});
        mu[1] = new Matrix(new double[][]{{1.0}});

        Matrix[] v = new Matrix[M];
        v[0] = new Matrix(new double[][]{{0.0}});
        v[1] = new Matrix(new double[][]{{0.0}});

        Matrix r = new Matrix(new double[][]{
            {0.0, 1.0},
            {1.0, 0.0}
        });

        LinearReductionParameters params = new LinearReductionParameters(M, N, K, mu, r, v);

        // Get upper bounds for queue 1
        Mapqn_solution solQ1 = Mapqn_bnd_lr.INSTANCE.solve(params, 1, 1);
        assertNotNull(solQ1, "Solution for queue 1 should not be null");

        double U1_UB = solQ1.getUtilization(1, 1);
        assertTrue(U1_UB >= 0, "Queue 1 utilization should be non-negative");
        assertTrue(U1_UB <= 1, "Queue 1 utilization should not exceed 1");

        // Get upper bounds for queue 2
        Mapqn_solution solQ2 = Mapqn_bnd_lr.INSTANCE.solve(params, 2, 1);
        assertNotNull(solQ2, "Solution for queue 2 should not be null");

        double U2_UB = solQ2.getUtilization(2, 1);
        assertTrue(U2_UB >= 0, "Queue 2 utilization should be non-negative");
        assertTrue(U2_UB <= 1, "Queue 2 utilization should not exceed 1");

        // Upper bound must be >= exact value
        assertTrue(U1_UB >= EXACT_U1 - 0.01,
                  String.format("Upper bound U1=%.4f should be >= exact=%.4f", U1_UB, EXACT_U1));
        assertTrue(U2_UB >= EXACT_U2 - 0.01,
                  String.format("Upper bound U2=%.4f should be >= exact=%.4f", U2_UB, EXACT_U2));

        // For symmetric network, bounds should be similar
        double diff = Math.abs(U1_UB - U2_UB);
        assertTrue(diff < 0.1,
                  String.format("Symmetric network should have similar bounds: U1=%.4f, U2=%.4f", U1_UB, U2_UB));

        // LR bounds quality check - bound gap should be reasonable (within 2x of exact for this simple case)
        assertTrue(U1_UB <= EXACT_U1 * 3,
                  String.format("Upper bound U1=%.4f should not exceed 3x exact=%.4f", U1_UB, EXACT_U1));
    }

    /**
     * Test Mapqn_bnd_lr with asymmetric 3-queue network.
     *
     * Reference (MVA exact): U1=0.8620, U2=0.2586, U3=0.5172
     * for N=3, mu=[1,2,0.5], asymmetric routing
     * LR bounds give upper bounds on utilization.
     */
    @Test
    @Timeout(value = 2, unit = TimeUnit.MINUTES)
    public void testMapqn_bnd_lr_asymmetric_3queue() {
        // Reference values from MVA
        final double EXACT_U1 = 0.8620102215;
        final double EXACT_U2 = 0.2586030664;
        final double EXACT_U3 = 0.5172061329;

        int M = 3;
        int N = 3;
        int[] K = {1, 1, 1};

        Matrix[] mu = new Matrix[M];
        mu[0] = new Matrix(new double[][]{{1.0}});   // rate 1
        mu[1] = new Matrix(new double[][]{{2.0}});   // rate 2 (faster)
        mu[2] = new Matrix(new double[][]{{0.5}});   // rate 0.5 (slower)

        Matrix[] v = new Matrix[M];
        v[0] = new Matrix(new double[][]{{0.0}});
        v[1] = new Matrix(new double[][]{{0.0}});
        v[2] = new Matrix(new double[][]{{0.0}});

        // Asymmetric routing from exactvalid.m
        Matrix r = new Matrix(new double[][]{
            {0.1, 0.6, 0.3},
            {1.0, 0.0, 0.0},
            {1.0, 0.0, 0.0}
        });

        LinearReductionParameters params = new LinearReductionParameters(M, N, K, mu, r, v);

        // Test all queues - get upper bounds
        double[] exactU = {EXACT_U1, EXACT_U2, EXACT_U3};
        double[] upperBounds = new double[M];

        for (int i = 1; i <= M; i++) {
            Mapqn_solution sol = Mapqn_bnd_lr.INSTANCE.solve(params, i, 1);
            assertNotNull(sol, String.format("Solution for queue %d should not be null", i));

            upperBounds[i-1] = sol.getUtilization(i, 1);
            assertTrue(upperBounds[i-1] >= 0, String.format("Queue %d utilization should be non-negative", i));
            assertTrue(upperBounds[i-1] <= 1, String.format("Queue %d utilization should not exceed 1", i));

            // Upper bound must be >= exact value
            assertTrue(upperBounds[i-1] >= exactU[i-1] - 0.01,
                      String.format("Upper bound for queue %d (%.4f) should be >= exact (%.4f)",
                                   i, upperBounds[i-1], exactU[i-1]));
        }

        // Verify routing conservation
        for (int i = 0; i < M; i++) {
            double rowSum = 0;
            for (int j = 0; j < M; j++) {
                rowSum += r.get(i, j);
            }
            assertEquals(1.0, rowSum, 1e-6,
                        String.format("Routing probabilities from queue %d should sum to 1", i + 1));
        }
    }

    /**
     * Test Mapqn_bnd_lr with different service rates (2-queue).
     *
     * Reference (MVA exact): U1=0.4667, U2=0.9333
     * for N=3, mu=[1,0.5], tandem routing
     * LR bounds give upper bounds on utilization.
     */
    @Test
    @Timeout(value = 2, unit = TimeUnit.MINUTES)
    public void testMapqn_bnd_lr_different_rates() {
        // Reference values from MVA
        final double EXACT_U1 = 0.4666666667;
        final double EXACT_U2 = 0.9333333333;

        int M = 2;
        int N = 3;
        int[] K = {1, 1};

        Matrix[] mu = new Matrix[M];
        mu[0] = new Matrix(new double[][]{{1.0}});   // rate 1
        mu[1] = new Matrix(new double[][]{{0.5}});   // rate 0.5 (slower, bottleneck)

        Matrix[] v = new Matrix[M];
        v[0] = new Matrix(new double[][]{{0.0}});
        v[1] = new Matrix(new double[][]{{0.0}});

        Matrix r = new Matrix(new double[][]{
            {0.0, 1.0},
            {1.0, 0.0}
        });

        LinearReductionParameters params = new LinearReductionParameters(M, N, K, mu, r, v);

        // Get upper bounds
        Mapqn_solution solQ1 = Mapqn_bnd_lr.INSTANCE.solve(params, 1, 1);
        Mapqn_solution solQ2 = Mapqn_bnd_lr.INSTANCE.solve(params, 2, 1);

        assertNotNull(solQ1, "Solution for queue 1 should not be null");
        assertNotNull(solQ2, "Solution for queue 2 should not be null");

        double U1_UB = solQ1.getUtilization(1, 1);
        double U2_UB = solQ2.getUtilization(2, 1);

        assertTrue(U1_UB >= 0 && U1_UB <= 1, String.format("Queue 1 utilization (%.4f) should be in [0,1]", U1_UB));
        assertTrue(U2_UB >= 0 && U2_UB <= 1, String.format("Queue 2 utilization (%.4f) should be in [0,1]", U2_UB));

        // Upper bounds must be >= exact values
        assertTrue(U1_UB >= EXACT_U1 - 0.01,
                  String.format("Upper bound U1=%.4f should be >= exact=%.4f", U1_UB, EXACT_U1));
        assertTrue(U2_UB >= EXACT_U2 - 0.01,
                  String.format("Upper bound U2=%.4f should be >= exact=%.4f", U2_UB, EXACT_U2));
    }

    /**
     * Test Mapqn_bnd_lr with 2-phase MAP service at one queue.
     * Based on testsigme.m from 08_mapqn_sigme08.
     * LR bounds give upper bounds on utilization.
     */
    @Test
    @Timeout(value = 2, unit = TimeUnit.MINUTES)
    public void testMapqn_bnd_lr_with_map() {
        int M = 2;
        int N = 3;
        int[] K = {1, 2};  // Queue 1: 1 phase, Queue 2: 2 phases (MAP)

        Matrix[] mu = new Matrix[M];
        mu[0] = new Matrix(new double[][]{{1.0}});  // Queue 1: exponential rate 1

        // Queue 2: 2-phase MAP (D1 matrix)
        mu[1] = new Matrix(new double[][]{
            {0.5, 0.3},
            {0.2, 0.8}
        });

        Matrix[] v = new Matrix[M];
        v[0] = new Matrix(new double[][]{{0.0}});

        // Queue 2: background transitions (D0 off-diagonal)
        v[1] = new Matrix(new double[][]{
            {0.0, 0.1},
            {0.05, 0.0}
        });

        Matrix r = new Matrix(new double[][]{
            {0.0, 1.0},
            {1.0, 0.0}
        });

        LinearReductionParameters params = new LinearReductionParameters(M, N, K, mu, r, v);

        // Test queue 1 (single phase)
        Mapqn_solution q1Sol = Mapqn_bnd_lr.INSTANCE.solve(params, 1, 1);
        assertNotNull(q1Sol, "Solution for queue 1 should not be null");
        double q1Util = q1Sol.getUtilization(1, 1);
        assertTrue(q1Util >= 0 && q1Util <= 1,
                  String.format("Queue 1 utilization (%.4f) should be in [0,1]", q1Util));

        // Test queue 2, phase 1
        Mapqn_solution q2Phase1Sol = Mapqn_bnd_lr.INSTANCE.solve(params, 2, 1);
        assertNotNull(q2Phase1Sol, "Solution for queue 2, phase 1 should not be null");
        double q2Phase1Util = q2Phase1Sol.getUtilization(2, 1);
        assertTrue(q2Phase1Util >= 0 && q2Phase1Util <= 1,
                  String.format("Queue 2 phase 1 utilization (%.4f) should be in [0,1]", q2Phase1Util));

        // Test queue 2, phase 2
        Mapqn_solution q2Phase2Sol = Mapqn_bnd_lr.INSTANCE.solve(params, 2, 2);
        assertNotNull(q2Phase2Sol, "Solution for queue 2, phase 2 should not be null");
        double q2Phase2Util = q2Phase2Sol.getUtilization(2, 2);
        assertTrue(q2Phase2Util >= 0 && q2Phase2Util <= 1,
                  String.format("Queue 2 phase 2 utilization (%.4f) should be in [0,1]", q2Phase2Util));

        // With N=3 jobs, all queues should have positive utilization (upper bounds)
        assertTrue(q1Util > 0, String.format("Queue 1 should have positive utilization, got %.4f", q1Util));
        assertTrue(q2Phase1Util > 0 || q2Phase2Util > 0,
                  String.format("Queue 2 should have positive utilization, got phase1=%.4f, phase2=%.4f",
                               q2Phase1Util, q2Phase2Util));
    }

    // ========== Cross-validation: RSRD vs BAS ==========

    /**
     * Test that RSRD and BAS give consistent results for networks without blocking.
     *
     * When F >= N for all queues (no blocking possible), both methods should
     * produce the same bounds.
     */
    @Test
    @Timeout(value = 2, unit = TimeUnit.MINUTES)
    public void testRSRD_vs_BAS_consistency() {
        final double EXPECTED_U = 2.0 / 3.0;  // N/(N+1) for N=2

        int M = 2;
        int N = 2;
        int[] K = {1, 1};
        int[] F = {2, 2};

        Matrix[] mu = new Matrix[M];
        mu[0] = new Matrix(new double[][]{{1.0}});
        mu[1] = new Matrix(new double[][]{{1.0}});

        Matrix[] v = new Matrix[M];
        v[0] = new Matrix(new double[][]{{0.0}});
        v[1] = new Matrix(new double[][]{{0.0}});

        Matrix r = new Matrix(new double[][]{
            {0.0, 1.0},
            {1.0, 0.0}
        });

        // RSRD
        double[][] alpha = new double[M][];
        for (int i = 0; i < M; i++) {
            alpha[i] = new double[N];
            Arrays.fill(alpha[i], 1.0);
        }

        Mapqn_qr_bounds_rsrd_parameters rsrdParams = new Mapqn_qr_bounds_rsrd_parameters(
            M, N, F, K, mu, v, alpha, r);

        Mapqn_solution rsrdMin = Mapqn_qr_bounds_rsrd.INSTANCE.solve(rsrdParams, 1, "min");
        Mapqn_solution rsrdMax = Mapqn_qr_bounds_rsrd.INSTANCE.solve(rsrdParams, 1, "max");

        // BAS
        int MR = 1;
        int f = 1;
        Matrix BB = new Matrix(new double[][]{{0, 0}});
        Matrix MM = new Matrix(new double[][]{{0, 0}});
        Matrix MM1 = new Matrix(new double[][]{{0, 0}});
        int[] ZZ = {0};

        Mapqn_qr_bounds_bas_parameters basParams = new Mapqn_qr_bounds_bas_parameters(
            M, N, MR, f, K, F, MM, MM1, ZZ, BB, mu, v, r);

        Mapqn_solution basMin = Mapqn_qr_bounds_bas.INSTANCE.solve(basParams, 1, "min");
        Mapqn_solution basMax = Mapqn_qr_bounds_bas.INSTANCE.solve(basParams, 1, "max");

        // Both should converge
        assertFalse(Double.isNaN(rsrdMin.getObjectiveValue()), "RSRD min should converge");
        assertFalse(Double.isNaN(basMin.getObjectiveValue()), "BAS min should converge");

        // Both should match expected value
        assertEquals(EXPECTED_U, rsrdMin.getObjectiveValue(), MATLAB_TOLERANCE,
            "RSRD should match expected");
        assertEquals(EXPECTED_U, basMin.getObjectiveValue(), MATLAB_TOLERANCE,
            "BAS should match expected");

        // Both methods should give consistent results
        assertEquals(rsrdMin.getObjectiveValue(), basMin.getObjectiveValue(), MATLAB_TOLERANCE,
            "RSRD and BAS should give consistent lower bounds");
        assertEquals(rsrdMax.getObjectiveValue(), basMax.getObjectiveValue(), MATLAB_TOLERANCE,
            "RSRD and BAS should give consistent upper bounds");
    }

    // ========== RSRD Edge Cases ==========

    /**
     * Test RSRD bounds for asymmetric 3-queue network M=3, N=3.
     *
     * MATLAB Reference (qrf_rsrd.m):
     * - LB = 0.8565395986
     * - UB = 0.8648455262
     * - U1 = 0.856540, U2 = 0.256962, U3 = 0.513924 (at minimum)
     *
     * This tests the non-trivial case with different service rates and asymmetric routing.
     *
     * NOTE: Currently disabled - Apache Commons Math SimplexSolver cannot converge
     * for M>2 networks due to simplex cycling issues. The constraint system is
     * correct (matches MATLAB's qrf_rsrd.m exactly), but the SimplexSolver's
     * simplex algorithm cycles indefinitely. MATLAB's linprog with interior-point
     * algorithm handles this correctly. Use LR bounds (testMapqn_bnd_lr_asymmetric_3queue)
     * as an alternative for M>2 cases.
     */
    @Test
    @Timeout(value = 2, unit = TimeUnit.MINUTES)
    public void testRSRD_asymmetric_3queue_M3_N3() {
        // MATLAB reference values
        final double MATLAB_LB = 0.8565395986;
        final double MATLAB_UB = 0.8648455262;

        int M = 3;
        int N = 3;
        int[] K = {1, 1, 1};
        int[] F = {3, 3, 3};

        Matrix[] mu = new Matrix[M];
        mu[0] = new Matrix(new double[][]{{1.0}});
        mu[1] = new Matrix(new double[][]{{2.0}});
        mu[2] = new Matrix(new double[][]{{0.5}});

        Matrix[] v = new Matrix[M];
        v[0] = new Matrix(new double[][]{{0.0}});
        v[1] = new Matrix(new double[][]{{0.0}});
        v[2] = new Matrix(new double[][]{{0.0}});

        // Asymmetric routing
        Matrix r = new Matrix(new double[][]{
            {0.1, 0.6, 0.3},
            {1.0, 0.0, 0.0},
            {1.0, 0.0, 0.0}
        });

        double[][] alpha = new double[M][];
        for (int i = 0; i < M; i++) {
            alpha[i] = new double[N];
            Arrays.fill(alpha[i], 1.0);
        }

        Mapqn_qr_bounds_rsrd_parameters params = new Mapqn_qr_bounds_rsrd_parameters(
            M, N, F, K, mu, v, alpha, r);

        Mapqn_solution solMin = Mapqn_qr_bounds_rsrd.INSTANCE.solve(params, 1, "min");
        Mapqn_solution solMax = Mapqn_qr_bounds_rsrd.INSTANCE.solve(params, 1, "max");

        assertNotNull(solMin, "Min solution should not be null");
        assertNotNull(solMax, "Max solution should not be null");

        double lineLB = solMin.getObjectiveValue();
        double lineUB = solMax.getObjectiveValue();

        assertFalse(Double.isNaN(lineLB), "Lower bound should not be NaN");
        assertFalse(Double.isNaN(lineUB), "Upper bound should not be NaN");

        // Basic consistency checks
        assertTrue(lineLB >= 0, "Lower bound should be non-negative");
        assertTrue(lineUB <= 1, "Upper bound should not exceed 1");
        assertTrue(lineLB <= lineUB + 1e-6, "Lower bound should not exceed upper bound");

        // Validate against MATLAB reference (with tolerance for numerical differences)
        assertEquals(MATLAB_LB, lineLB, MATLAB_TOLERANCE,
            String.format("RSRD LB (%.6f) should match MATLAB reference (%.6f)", lineLB, MATLAB_LB));
        assertEquals(MATLAB_UB, lineUB, MATLAB_TOLERANCE,
            String.format("RSRD UB (%.6f) should match MATLAB reference (%.6f)", lineUB, MATLAB_UB));
    }

    /**
     * Test RSRD with MAP (2-phase) service process.
     *
     * This validates that the implementations handle non-exponential
     * (Markovian Arrival Process) service correctly.
     */
    @Test
    @Timeout(value = 2, unit = TimeUnit.MINUTES)
    public void testRSRD_with_MAP_service() {
        int M = 2;
        int N = 3;
        int[] K = {2, 1};  // Queue 1 has 2 phases (MAP), Queue 2 has 1 phase (Exp)
        int[] F = {3, 3};

        // Queue 1: 2-phase MAP service
        Matrix[] mu = new Matrix[M];
        mu[0] = new Matrix(new double[][]{
            {0.5, 0.3},
            {0.2, 0.8}
        });
        mu[1] = new Matrix(new double[][]{{1.0}});

        // Background transitions for MAP
        Matrix[] v = new Matrix[M];
        v[0] = new Matrix(new double[][]{
            {0.0, 0.1},
            {0.05, 0.0}
        });
        v[1] = new Matrix(new double[][]{{0.0}});

        Matrix r = new Matrix(new double[][]{
            {0.0, 1.0},
            {1.0, 0.0}
        });

        double[][] alpha = new double[M][];
        for (int i = 0; i < M; i++) {
            alpha[i] = new double[N];
            Arrays.fill(alpha[i], 1.0);
        }

        Mapqn_qr_bounds_rsrd_parameters params = new Mapqn_qr_bounds_rsrd_parameters(
            M, N, F, K, mu, v, alpha, r);

        Mapqn_solution solMin = Mapqn_qr_bounds_rsrd.INSTANCE.solve(params, 1, "min");
        Mapqn_solution solMax = Mapqn_qr_bounds_rsrd.INSTANCE.solve(params, 1, "max");

        assertNotNull(solMin, "Min solution should not be null");
        assertNotNull(solMax, "Max solution should not be null");

        double LB = solMin.getObjectiveValue();
        double UB = solMax.getObjectiveValue();

        assertFalse(Double.isNaN(LB), "Lower bound should not be NaN");
        assertFalse(Double.isNaN(UB), "Upper bound should not be NaN");

        // Basic validity
        assertTrue(LB >= 0, "Lower bound should be non-negative");
        assertTrue(UB <= 1, "Upper bound should not exceed 1");
        assertTrue(LB <= UB + 1e-6, "Lower bound should not exceed upper bound");

        // With N=3 jobs, utilization should be positive
        assertTrue(LB > 0, "With 3 jobs, utilization lower bound should be positive");
    }

    /**
     * Test BAS with actual blocking configurations (MR=2).
     * 2-queue tandem with blocking, validates MM indexing fix (1-based to 0-based).
     *
     * Configuration: M=2, N=3, f=1, F=[2,3], K=[2,2], MR=2
     * MATLAB reference (qrf_bas.m):
     *   Lower bound: 0.718653
     *   Upper bound: 0.823197
     */
    @Test
    @Timeout(value = 2, unit = TimeUnit.MINUTES)
    public void testBAS_blocking_M2_MR2() {
        final double MATLAB_LB = 0.718653;
        final double MATLAB_UB = 0.823197;

        int M = 2;
        int N = 3;
        int f = 1;  // Finite capacity queue (1-based)
        int[] K = {2, 2};
        int[] F = {2, 3};

        // Service rates mu{i}(k,h)
        Matrix[] mu = new Matrix[M];
        mu[0] = new Matrix(new double[][]{
            {1.0, 0.1},
            {0.1, 0.5}
        });
        mu[1] = new Matrix(new double[][]{
            {0.8, 0.2},
            {0.2, 0.6}
        });

        // Background transition rates (zeros)
        Matrix[] v = new Matrix[M];
        for (int i = 0; i < M; i++) {
            v[i] = new Matrix(new double[][]{{0.0, 0.0}, {0.0, 0.0}});
        }

        // Routing probabilities (tandem)
        Matrix r = new Matrix(new double[][]{
            {0.0, 1.0},
            {1.0, 0.0}
        });

        // Blocking configurations (MR=2)
        int MR = 2;

        // BB(m, i) = 1 if queue i is blocked in configuration m
        Matrix BB = new Matrix(new double[][]{
            {0, 0},   // m=0: no blocking
            {0, 1}    // m=1: queue 2 blocked
        });

        // MM(m, order) = queue index (1-based) at position 'order' in blocking list
        Matrix MM = new Matrix(new double[][]{
            {0, 0},   // m=0: no blocking
            {2, 0}    // m=1: queue 2 first (1-based)
        });

        // ZZ(m) = number of blocked queues in configuration m
        int[] ZZ = {0, 1};

        // MM1(m, j) = extended blocking order info
        Matrix MM1 = new Matrix(new double[][]{
            {0, 0},
            {0, 0}
        });

        Mapqn_qr_bounds_bas_parameters params = new Mapqn_qr_bounds_bas_parameters(
            M, N, MR, f, K, F, MM, MM1, ZZ, BB, mu, v, r);

        Mapqn_solution solMin = Mapqn_qr_bounds_bas.INSTANCE.solve(params, 1, "min");
        Mapqn_solution solMax = Mapqn_qr_bounds_bas.INSTANCE.solve(params, 1, "max");

        assertNotNull(solMin, "Min solution should not be null");
        assertNotNull(solMax, "Max solution should not be null");

        double lineLB = solMin.getObjectiveValue();
        double lineUB = solMax.getObjectiveValue();

        assertFalse(Double.isNaN(lineLB), "Lower bound should not be NaN");
        assertFalse(Double.isNaN(lineUB), "Upper bound should not be NaN");

        // Basic consistency
        assertTrue(lineLB >= 0, "Lower bound should be non-negative");
        assertTrue(lineUB <= 1, "Upper bound should not exceed 1");
        assertTrue(lineLB <= lineUB + 1e-6, "Lower bound should not exceed upper bound");

        // Validate against MATLAB reference
        assertEquals(MATLAB_LB, lineLB, MATLAB_TOLERANCE,
            String.format("BAS LB (%.6f) should match MATLAB reference (%.6f)", lineLB, MATLAB_LB));
        assertEquals(MATLAB_UB, lineUB, MATLAB_TOLERANCE,
            String.format("BAS UB (%.6f) should match MATLAB reference (%.6f)", lineUB, MATLAB_UB));
    }
}
