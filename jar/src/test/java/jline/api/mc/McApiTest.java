/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.mc;

import jline.util.Pair;
import jline.util.Triple;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Analytic validation of the Markov-chain APIs (jline.api.mc) on chains with
 * closed-form stationary distributions:
 * - 3-state birth-death CTMC (lambda=1, mu=2): pi = [4/7, 2/7, 1/7], reversible;
 * - 2-state DTMC [[0.9,0.1],[0.2,0.8]]: pi = [2/3, 1/3], reversible.
 */
public class McApiTest {

    private static final double TOL = 1e-8;
    private static final double ITER_TOL = 1e-2;

    private static Matrix birthDeathQ() {
        Matrix q = new Matrix(3, 3);
        q.set(0, 0, -1.0); q.set(0, 1, 1.0);
        q.set(1, 0, 2.0);  q.set(1, 1, -3.0); q.set(1, 2, 1.0);
        q.set(2, 1, 2.0);  q.set(2, 2, -2.0);
        return q;
    }

    private static final double[] PI_BD = {4.0 / 7.0, 2.0 / 7.0, 1.0 / 7.0};

    /**
     * Asymmetric 4-state chain. Asymmetry matters: a symmetric chain has a uniform
     * stationary vector, which hides every index-space error.
     */
    private static Matrix fourStateQ() {
        Matrix q = new Matrix(4, 4);
        q.set(0, 1, 2.0);
        q.set(1, 0, 1.0);
        q.set(1, 2, 3.0);
        q.set(2, 1, 1.5);
        q.set(2, 3, 0.7);
        q.set(3, 2, 2.2);
        for (int i = 0; i < 4; i++) {
            double s = 0.0;
            for (int j = 0; j < 4; j++) {
                if (i != j) {
                    s += q.get(i, j);
                }
            }
            q.set(i, i, -s);
        }
        return q;
    }

    private static double maxAbsDiff(Matrix a, Matrix b) {
        double m = 0.0;
        for (int i = 0; i < a.getNumElements(); i++) {
            m = Math.max(m, Math.abs(a.get(i) - b.get(i)));
        }
        return m;
    }

    private static Matrix twoStateP() {
        Matrix p = new Matrix(2, 2);
        p.set(0, 0, 0.9); p.set(0, 1, 0.1);
        p.set(1, 0, 0.2); p.set(1, 1, 0.8);
        return p;
    }

    private static final double[] PI_2ST = {2.0 / 3.0, 1.0 / 3.0};

    private static void assertProbabilityVector(Matrix pi, String context) {
        double sum = 0;
        for (int i = 0; i < pi.getNumElements(); i++) {
            double v = pi.get(i);
            assertTrue(v >= -TOL, context + ": negative probability " + v);
            sum += v;
        }
        assertEquals(1.0, sum, 1e-6, context + ": probabilities must sum to 1");
    }

    @Test
    public void ctmcSolveMatchesClosedForm() {
        Matrix pi = Ctmc_solve.ctmc_solve(birthDeathQ());
        assertProbabilityVector(pi, "ctmc_solve");
        for (int i = 0; i < 3; i++) {
            assertEquals(PI_BD[i], pi.get(i), TOL, "ctmc_solve pi[" + i + "]");
        }
    }

    @Test
    public void ctmcMakeinfgenRepairsDiagonal() {
        Matrix broken = birthDeathQ();
        broken.set(1, 1, 0.0); // corrupt a diagonal entry
        Matrix fixed = Ctmc_makeinfgen.ctmc_makeinfgen(broken);
        for (int i = 0; i < 3; i++) {
            double rowSum = 0;
            for (int j = 0; j < 3; j++) {
                rowSum += fixed.get(i, j);
                if (i != j) {
                    assertTrue(fixed.get(i, j) >= 0, "off-diagonal must be nonnegative");
                }
            }
            assertEquals(0.0, rowSum, TOL, "infinitesimal generator row sum");
        }
    }

    @Test
    public void ctmcTimereverseOfReversibleChainIsItself() {
        Matrix q = birthDeathQ();
        Matrix qRev = Ctmc_timereverse.ctmc_timereverse(q);
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                assertEquals(q.get(i, j), qRev.get(i, j), 1e-8,
                        "reversible chain must equal its time reverse at (" + i + "," + j + ")");
            }
        }
    }

    @Test
    public void ctmcKolmogorovCriterionHoldsForBirthDeath() {
        assertTrue(Ctmc_testpf_kolmogorov.ctmc_testpf_kolmogorov(birthDeathQ()),
                "birth-death chains are reversible");
    }

    @Test
    public void ctmcUniformizationConvergesToStationary() {
        Matrix pi0 = new Matrix(1, 3);
        pi0.set(0, 0, 1.0);
        Matrix piT = Ctmc_uniformization.ctmc_uniformization(pi0, birthDeathQ(), 50.0);
        for (int i = 0; i < 3; i++) {
            assertEquals(PI_BD[i], piT.get(i), 1e-6, "uniformization pi[" + i + "] at t=50");
        }
    }

    @Test
    public void ctmcTransientConvergesToStationary() {
        Matrix pi0 = new Matrix(1, 3);
        pi0.set(0, 0, 1.0);
        Pair<double[], List<double[]>> sol =
                Ctmc_transient.ctmc_transient(birthDeathQ(), pi0, 50.0);
        assertNotNull(sol, "ctmc_transient returned null");
        List<double[]> states = sol.getRight();
        assertTrue(states.size() > 1, "transient trajectory must have multiple points");
        double[] last = states.get(states.size() - 1);
        for (int i = 0; i < 3; i++) {
            assertEquals(PI_BD[i], last[i], 1e-4, "transient endpoint pi[" + i + "]");
        }
    }

    @Test
    public void ctmcRandomizationPreservesStationaryDistribution() {
        Pair<Matrix, Double> randomized = Ctmc_randomization.ctmc_randomization(birthDeathQ());
        Matrix p = randomized.getLeft();
        // Uniformized DTMC is stochastic and has the same stationary vector
        for (int i = 0; i < 3; i++) {
            double rowSum = 0;
            for (int j = 0; j < 3; j++) {
                assertTrue(p.get(i, j) >= -TOL, "randomization: negative transition prob");
                rowSum += p.get(i, j);
            }
            assertEquals(1.0, rowSum, TOL, "randomization: row must be stochastic");
        }
        for (int j = 0; j < 3; j++) {
            double v = 0;
            for (int i = 0; i < 3; i++) {
                v += PI_BD[i] * p.get(i, j);
            }
            assertEquals(PI_BD[j], v, TOL, "randomization: pi P != pi at state " + j);
        }
    }

    @Test
    public void ctmcTimeAverageConvergesToStationary() {
        Matrix pi0 = new Matrix(1, 3);
        pi0.set(0, 0, 1.0);
        Pair<Matrix, Matrix> ta = Ctmc_timeaverage.ctmc_timeaverage(pi0, birthDeathQ(), 1000.0);
        Matrix avg = ta.getLeft();
        for (int i = 0; i < 3; i++) {
            assertEquals(PI_BD[i], avg.get(i), ITER_TOL, "time average pi[" + i + "]");
        }
    }

    @Test
    public void ctmcSolveReducibleAgreesOnIrreducibleChain() {
        Pair<Matrix, List<List<Integer>>> sol =
                Ctmc_solve_reducible.ctmc_solve_reducible(birthDeathQ());
        Matrix pi = sol.getLeft();
        for (int i = 0; i < 3; i++) {
            assertEquals(PI_BD[i], pi.get(i), 1e-6, "solve_reducible pi[" + i + "]");
        }
    }

    /**
     * The aggregation methods work internally in permuted (macrostate-major) ordering
     * and must map their result back to the original state ordering. A partition whose
     * concatenation is the identity permutation cannot detect a missing or wrong
     * mapping, so pin a partition that reorders the states: listing the same physical
     * blocks in the reverse order yields v = [2, 0, 1], and the answer must be
     * unchanged because the partition itself is unchanged.
     */
    @Test
    public void ctmcKmsIsInvariantToMacroBlockOrdering() {
        List<List<Integer>> forward = new ArrayList<List<Integer>>();
        forward.add(Arrays.asList(0, 1));
        forward.add(Arrays.asList(2));

        List<List<Integer>> reversed = new ArrayList<List<Integer>>();
        reversed.add(Arrays.asList(2));
        reversed.add(Arrays.asList(0, 1));

        // Check at a SMALL step count as well as at convergence. The iteration
        // contracts to its fixed point from any starting vector, so seeding it in the
        // wrong index space is invisible after many steps and only shows up early.
        int[] stepCounts = {1, 3, 200};
        for (int s = 0; s < stepCounts.length; s++) {
            int steps = stepCounts[s];
            Triple<Matrix, Double, Double> kmsFwd =
                    Ctmc_kms.ctmc_kms(birthDeathQ(), forward, steps);
            Triple<Matrix, Double, Double> kmsRev =
                    Ctmc_kms.ctmc_kms(birthDeathQ(), reversed, steps);

            assertProbabilityVector(kmsRev.getFirst(), "kms reversed-block-order");
            for (int i = 0; i < 3; i++) {
                assertEquals(kmsFwd.getFirst().get(i), kmsRev.getFirst().get(i), TOL,
                        "kms must not depend on the order the macro-blocks are listed in"
                                + " (steps=" + steps + ") pi[" + i + "]");
            }
        }

        // At convergence both orderings must also land on the exact vector, in
        // ORIGINAL state ordering rather than permuted ordering.
        Triple<Matrix, Double, Double> kmsRevConverged =
                Ctmc_kms.ctmc_kms(birthDeathQ(), reversed, 200);
        for (int i = 0; i < 3; i++) {
            assertEquals(PI_BD[i], kmsRevConverged.getFirst().get(i), ITER_TOL,
                    "kms reversed-block-order pi[" + i + "]");
        }
    }

    /**
     * The macro-states are the SETS of states listed in MS, not merely their sizes.
     * Two partitions with identical block sizes but different members are different
     * partitions and must not produce identical iterates. This is the companion to
     * {@link #ctmcKmsIsInvariantToMacroBlockOrdering()}: a routine that block-indexes
     * contiguous offsets rather than the MS members reads the members nowhere, and so
     * silently solves for a contiguous partition whatever the caller asked for.
     */
    @Test
    public void ctmcAggregationHonoursMacroStateMembership() {
        // Two partitions of {0..3} with the same block sizes (2, 2), different members.
        List<List<Integer>> byHalves = new ArrayList<List<Integer>>();
        byHalves.add(Arrays.asList(0, 1));
        byHalves.add(Arrays.asList(2, 3));

        List<List<Integer>> interleaved = new ArrayList<List<Integer>>();
        interleaved.add(Arrays.asList(0, 2));
        interleaved.add(Arrays.asList(1, 3));

        Matrix q = fourStateQ();
        int steps = 3;

        Matrix tkHalves = Ctmc_takahashi.ctmc_takahashi(q, byHalves, steps).getFirst();
        Matrix tkInter = Ctmc_takahashi.ctmc_takahashi(q, interleaved, steps).getFirst();
        assertTrue(maxAbsDiff(tkHalves, tkInter) > 1e-6,
                "takahashi must honour macro-state membership, not just block sizes");

        Matrix kmsHalves = Ctmc_kms.ctmc_kms(q, byHalves, steps).getFirst();
        Matrix kmsInter = Ctmc_kms.ctmc_kms(q, interleaved, steps).getFirst();
        assertTrue(maxAbsDiff(kmsHalves, kmsInter) > 1e-6,
                "kms must honour macro-state membership, not just block sizes");
    }

    @Test
    public void ctmcAggregationMethodsRecoverStationaryVector() {
        List<List<Integer>> partition = new ArrayList<List<Integer>>();
        partition.add(Arrays.asList(0, 1));
        partition.add(Arrays.asList(2));

        Triple<Matrix, Double, Double> tk =
                Ctmc_takahashi.ctmc_takahashi(birthDeathQ(), partition, 200);
        assertProbabilityVector(tk.getFirst(), "takahashi");
        for (int i = 0; i < 3; i++) {
            assertEquals(PI_BD[i], tk.getFirst().get(i), ITER_TOL, "takahashi pi[" + i + "]");
        }

        Triple<Matrix, Double, Double> kms =
                Ctmc_kms.ctmc_kms(birthDeathQ(), partition, 200);
        assertProbabilityVector(kms.getFirst(), "kms");
        for (int i = 0; i < 3; i++) {
            assertEquals(PI_BD[i], kms.getFirst().get(i), ITER_TOL, "kms pi[" + i + "]");
        }

        Triple<Matrix, Double, Double> courtois =
                Ctmc_courtois.ctmc_courtois(birthDeathQ(), partition);
        assertProbabilityVector(courtois.getFirst(), "courtois");
    }

    @Test
    public void ctmcRandProducesValidGenerator() {
        Matrix q = Ctmc_rand.ctmc_rand(4);
        for (int i = 0; i < 4; i++) {
            double rowSum = 0;
            for (int j = 0; j < 4; j++) {
                rowSum += q.get(i, j);
                if (i != j) {
                    assertTrue(q.get(i, j) >= 0, "ctmc_rand: negative off-diagonal");
                }
            }
            assertEquals(0.0, rowSum, 1e-8, "ctmc_rand: generator row sum");
        }
    }

    // ------------------------------------------------------------------
    // DTMC
    // ------------------------------------------------------------------

    @Test
    public void dtmcSolveMatchesClosedForm() {
        Matrix pi = Dtmc_solve.dtmc_solve(twoStateP());
        for (int i = 0; i < 2; i++) {
            assertEquals(PI_2ST[i], pi.get(i), TOL, "dtmc_solve pi[" + i + "]");
        }
    }

    @Test
    public void dtmcTimereverseOfReversibleChainIsItself() {
        Matrix p = twoStateP();
        Matrix pRev = Dtmc_timereverse.dtmc_timereverse(p);
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                assertEquals(p.get(i, j), pRev.get(i, j), 1e-8,
                        "2-state chains are reversible; time reverse must match");
            }
        }
    }

    @Test
    public void dtmcMakestochasticNormalizesRows() {
        Matrix broken = twoStateP();
        broken.set(0, 1, 0.3); // row 0 now sums to 1.2
        Matrix fixed = Dtmc_makestochastic.dtmc_makestochastic(broken);
        for (int i = 0; i < 2; i++) {
            double rowSum = 0;
            for (int j = 0; j < 2; j++) {
                rowSum += fixed.get(i, j);
                assertTrue(fixed.get(i, j) >= 0, "dtmc_makestochastic: negative entry");
            }
            assertEquals(1.0, rowSum, 1e-8, "dtmc_makestochastic: row sum");
        }
    }

    @Test
    public void dtmcIsFeasibleGradesStochasticMatrices() {
        // Returns the finest tolerance exponent (1..15) at which P is
        // stochastic; an exactly stochastic matrix grades at the maximum
        assertEquals(15, Dtmc_isfeasible.dtmc_isfeasible(twoStateP()),
                "an exactly stochastic matrix must grade at 1e-15");
        Matrix broken = twoStateP();
        broken.set(0, 1, 0.3); // row 0 sums to 1.2: infeasible at every tolerance
        assertEquals(0, Dtmc_isfeasible.dtmc_isfeasible(broken),
                "a non-stochastic matrix must be graded infeasible");
    }

    @Test
    public void dtmcSolveReducibleAgreesOnIrreducibleChain() {
        Pair<Matrix, List<List<Integer>>> sol =
                Dtmc_solve_reducible.dtmc_solve_reducible(twoStateP());
        for (int i = 0; i < 2; i++) {
            assertEquals(PI_2ST[i], sol.getLeft().get(i), 1e-6,
                    "dtmc_solve_reducible pi[" + i + "]");
        }
    }

    @Test
    public void dtmcStochcompCensorsToConditionalDistribution() {
        // Censor the 3-state uniformized birth-death chain on states {0,1}
        Pair<Matrix, Double> randomized = Ctmc_randomization.ctmc_randomization(birthDeathQ());
        Matrix p3 = randomized.getLeft();
        List<Integer> keep = Arrays.asList(0, 1);
        Matrix censored = Dtmc_stochcomp.dtmc_stochcomp(p3, keep);
        Matrix piC = Dtmc_solve.dtmc_solve(censored);
        double norm = PI_BD[0] + PI_BD[1];
        assertEquals(PI_BD[0] / norm, piC.get(0), 1e-8, "censored pi[0]");
        assertEquals(PI_BD[1] / norm, piC.get(1), 1e-8, "censored pi[1]");
    }

    @Test
    public void dtmcRandProducesStochasticMatrix() {
        Matrix p = Dtmc_rand.dtmc_rand(4);
        for (int i = 0; i < 4; i++) {
            double rowSum = 0;
            for (int j = 0; j < 4; j++) {
                assertTrue(p.get(i, j) >= 0, "dtmc_rand: negative entry");
                rowSum += p.get(i, j);
            }
            assertEquals(1.0, rowSum, 1e-8, "dtmc_rand: row sum");
        }
    }

    @Test
    public void dtmcSimulateOccupancyMatchesStationary() {
        Matrix pi0 = new Matrix(1, 2);
        pi0.set(0, 0, 1.0);
        int n = 50000;
        int[] states = Dtmc_simulate.dtmc_simulate(twoStateP(), pi0, n);
        assertEquals(n, states.length, "dtmc_simulate must return n samples");
        double occ0 = 0;
        for (int s : states) {
            assertTrue(s == 0 || s == 1, "dtmc_simulate: invalid state " + s);
            if (s == 0) {
                occ0 += 1.0;
            }
        }
        occ0 /= n;
        assertEquals(PI_2ST[0], occ0, 0.05, "dtmc_simulate occupancy of state 0");
    }

    @Test
    public void ctmcRelsolveNormalizesToStationary() {
        // ctmc_relsolve returns the UNNORMALIZED stationary vector relative to
        // the reference state (entry 0 == 1). Normalizing must recover the
        // exact stationary distribution.
        Object[] r = Ctmc_relsolve.ctmc_relsolve(birthDeathQ());
        Matrix p = (Matrix) r[0];
        double sum = 0.0;
        for (int i = 0; i < 3; i++) {
            assertTrue(p.get(i) >= -TOL, "relative solution must be nonnegative");
            sum += p.get(i);
        }
        for (int i = 0; i < 3; i++) {
            assertEquals(PI_BD[i], p.get(i) / sum, 1e-8,
                    "normalized relsolve pi[" + i + "]");
        }
    }

    @Test
    public void ctmcTimereverseSatisfiesReversalIdentity() {
        // The time-reversed generator qr must be a proper generator sharing the
        // stationary distribution, with pi_i * qr_ij = pi_j * q_ji exactly.
        Matrix q = new Matrix(3, 3);
        q.set(0, 0, -3); q.set(0, 1, 2); q.set(0, 2, 1);
        q.set(1, 0, 1); q.set(1, 1, -4); q.set(1, 2, 3);
        q.set(2, 0, 2); q.set(2, 1, 1); q.set(2, 2, -3);
        Matrix pi = Ctmc_solve.ctmc_solve(q);
        Matrix qr = Ctmc_timereverse.ctmc_timereverse(q);
        for (int i = 0; i < 3; i++) {
            double rowSum = 0.0;
            for (int j = 0; j < 3; j++) {
                rowSum += qr.get(i, j);
                assertEquals(pi.get(i) * qr.get(i, j), pi.get(j) * q.get(j, i), 1e-9,
                        "reversal identity at (" + i + "," + j + ")");
            }
            assertEquals(0.0, rowSum, 1e-9, "reversed generator row sum");
        }
    }

    @Test
    public void ctmcStochcompMatchesConditionalStationary() {
        // The stochastic complement onto a subset S has stationary vector equal
        // to the original stationary distribution conditioned on S.
        Matrix q = new Matrix(4, 4);
        double[][] v = {{-3, 1, 1, 1}, {2, -5, 2, 1}, {1, 1, -4, 2}, {1, 2, 1, -4}};
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                q.set(i, j, v[i][j]);
            }
        }
        Matrix pi = Ctmc_solve.ctmc_solve(q);
        jline.solvers.ctmc.SolverCTMC.StochCompResult r =
                Ctmc_stochcomp.ctmc_stochcomp(q, java.util.Arrays.asList(0.0, 1.0));
        Matrix piS = Ctmc_solve.ctmc_solve(r.S);
        double norm = pi.get(0) + pi.get(1);
        assertEquals(pi.get(0) / norm, piS.get(0), 1e-8, "stochcomp conditional pi[0]");
        assertEquals(pi.get(1) / norm, piS.get(1), 1e-8, "stochcomp conditional pi[1]");
    }

    /**
     * Assembles the linear system ctmc_solve poses: the last column of Q is replaced by
     * ones to carry the normalization, and the transposed system is solved against e_n.
     */
    private static Matrix[] normalizedSystem(Matrix q) {
        int n = q.getNumRows();
        Matrix a = q.copy();
        for (int i = 0; i < n; i++) {
            a.set(i, n - 1, 1.0);
        }
        Matrix b = new Matrix(n, 1);
        b.set(n - 1, 0, 1.0);
        return new Matrix[]{a.transpose(), b};
    }

    /** Generator of an M/M/1/K queue, whose stationary vector is a truncated geometric. */
    private static Matrix mm1kQ(double lambda, double mu, int K) {
        int n = K + 1;
        Matrix q = new Matrix(n, n);
        for (int i = 0; i < n; i++) {
            double s = 0.0;
            if (i < n - 1) {
                q.set(i, i + 1, lambda);
                s += lambda;
            }
            if (i > 0) {
                q.set(i, i - 1, mu);
                s += mu;
            }
            q.set(i, i, -s);
        }
        return q;
    }

    @Test
    public void ctmcGmresMatchesTruncatedGeometric() {
        // Analytic oracle rather than a recorded baseline: pi(i) proportional to rho^i.
        // K = 20000 is past the size at which a direct factorization is the intended
        // method, and rho < 1 is the case in which an unpreconditioned or naturally
        // ordered elimination overflows.
        for (double rho : new double[]{0.7, 1.0}) {
            for (int K : new int[]{10, 20000}) {
                Matrix[] ab = normalizedSystem(mm1kQ(rho, 1.0, K));
                Ctmc_gmres.GmresResult r = Ctmc_gmres.ctmc_gmres(ab[0], ab[1]);
                assertEquals(0, r.flag, "gmres flag for rho=" + rho + ", K=" + K);

                double s = 0.0;
                double[] pi = new double[K + 1];
                for (int i = 0; i <= K; i++) {
                    pi[i] = Math.pow(rho, i);
                    s += pi[i];
                }
                double err = 0.0;
                for (int i = 0; i <= K; i++) {
                    err = Math.max(err, Math.abs(r.x.get(i, 0) - pi[i] / s));
                }
                assertTrue(err < 1e-9, "gmres error " + err + " for rho=" + rho + ", K=" + K);
            }
        }
    }

    @Test
    public void ctmcGmresAgreesWithDirectSolve() {
        // The dispatch in ctmc_solve must be numerically invisible, so the two methods
        // have to agree far below any tolerance a caller would notice.
        Matrix[] chains = {birthDeathQ(), fourStateQ(), mm1kQ(0.85, 1.0, 200)};
        for (Matrix q : chains) {
            Matrix[] ab = normalizedSystem(q);
            Matrix xd = new Matrix(ab[1].getNumRows(), 1);
            assertTrue(Matrix.solveDirect(ab[0], ab[1], xd), "direct solve succeeded");
            Ctmc_gmres.GmresResult r = Ctmc_gmres.ctmc_gmres(ab[0], ab[1]);
            assertEquals(0, r.flag, "gmres flag");
            assertTrue(maxAbsDiff(xd, r.x) < 1e-9, "gmres vs direct: " + maxAbsDiff(xd, r.x));
        }
    }

    @Test
    public void ctmcGmresReportsNonConvergence() {
        // A caller may only trust the answer when the flag is zero. One restart cycle of
        // dimension one cannot converge on a 500-state chain, and the kernel has to say
        // so rather than return the iterate it happens to hold.
        Matrix[] ab = normalizedSystem(mm1kQ(0.9, 1.0, 500));
        Ctmc_gmres.GmresResult r = Ctmc_gmres.ctmc_gmres(ab[0], ab[1], 1e-12, 1, 1, null);
        assertTrue(r.flag != 0, "gmres must report non-convergence");
        assertTrue(r.relres > 1e-12, "non-convergent residual above tolerance");
    }

    @Test
    public void ctmcSolveDispatchIsNumericallyInvisible() {
        // The two methods must agree far below any tolerance a caller would notice, and
        // the answer must not jump as a model grows past the dispatch threshold.
        jline.solvers.SolverOptions direct = new jline.solvers.SolverOptions();
        direct.method = "direct";
        jline.solvers.SolverOptions iterative = new jline.solvers.SolverOptions();
        iterative.method = "gmres";

        for (int K : new int[]{50, 400}) {
            Matrix q = mm1kQ(0.8, 1.0, K);
            Matrix pd = Ctmc_solve.ctmc_solve(q, direct);
            Matrix pg = Ctmc_solve.ctmc_solve(q, iterative);
            assertTrue(maxAbsDiff(pd, pg) < 1e-9, "dispatch deviation " + maxAbsDiff(pd, pg));

            double s = 0.0;
            double[] pi = new double[K + 1];
            for (int i = 0; i <= K; i++) {
                pi[i] = Math.pow(0.8, i);
                s += pi[i];
            }
            double err = 0.0;
            for (int i = 0; i <= K; i++) {
                err = Math.max(err, Math.abs(pg.get(i) - pi[i] / s));
            }
            assertTrue(err < 1e-9, "gmres dispatch error " + err);
        }
    }

    @Test
    public void ctmcGmresSurvivesZeroDiagonal() {
        // A generator with no diagonal entries breaks the incomplete factorization. The
        // kernel must fall back rather than propagate the breakdown, and the deterministic
        // 3-cycle it is given here has a stationary vector that is exactly known.
        Matrix z = new Matrix(3, 3);
        z.set(0, 1, 1.0);
        z.set(1, 2, 1.0);
        z.set(2, 0, 1.0);
        Matrix b = new Matrix(3, 1);
        b.set(0, 0, 1.0);
        Ctmc_gmres.GmresResult r = Ctmc_gmres.ctmc_gmres(z, b);
        assertNotNull(r.x, "breakdown path returns a vector");
        assertTrue(Double.isFinite(r.x.get(0, 0)), "breakdown path returns finite entries");
    }
}
