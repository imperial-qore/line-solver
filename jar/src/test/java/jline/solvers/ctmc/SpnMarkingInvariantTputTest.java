package jline.solvers.ctmc;

import jline.lang.Network;
import jline.lang.state.TestSPNModels;
import jline.solvers.NetworkAvgTable;
import org.junit.jupiter.api.Test;

import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Throughput regression tests for marking-invariant SPNs, i.e. nets whose firing
 * outcome returns exactly what the enabling condition consumed (P1 -> T1 -> P1).
 *
 * Such nets are the worst case for the CTMC firing accounting because the place
 * marking never changes: QLen is trivially correct while every firing can be
 * missed. Three separate defects hid here, each reporting a plausible-looking
 * throughput on a net whose exact answer is known in closed form:
 *
 * <ul>
 *   <li>the firing-completion flag was inferred from "some place marking changed",
 *       which is never true for these nets (reported Tput = 0);</li>
 *   <li>the PHASE branch dropped the space_fired columns, so no phase-advance row
 *       could ever match the state space and multi-phase modes never progressed;</li>
 *   <li>the firing rates were bucketed by mode index instead of job class, so every
 *       mode after the first fell outside 0..nclasses-1 and was discarded.</li>
 * </ul>
 *
 * The fixtures documented the expected throughput but nothing asserted it, which is
 * how all three survived. These tests assert it.
 *
 * Exactness argument: with one token and every mode's enabling condition satisfied,
 * all modes are permanently enabled and race. Each mode is then a renewal process of
 * mean 1 (it restarts immediately on completion), so each contributes throughput 1
 * regardless of its distribution, and the transition throughput is the mode count.
 * This matches MATLAB spn_basic_closed and the native Python implementation.
 */
public class SpnMarkingInvariantTputTest {

    private static final double TOL = 1e-9;

    /**
     * Throughput here is read off a generator that carries immediate transitions at
     * the GlobalConstants.Immediate scale (1/Immediate = 1e-8), so the stationary
     * solve perturbs it by a few 1e-8 -- the exact value is approached, not hit.
     * QLen keeps TOL: it is a marking invariant, not a solve result, and must stay
     * exact. This stays far tighter than any of the defects the file guards, which
     * miss whole modes or report Tput = 0.
     */
    private static final double TPUT_TOL = 1e-6;

    private static double tputOf(Network model) {
        SolverCTMC solver = new SolverCTMC(model);
        NetworkAvgTable table = solver.getAvgTable();
        List<Double> tput = table.getTput();
        return tput.get(0);
    }

    /**
     * P1 -> T1 -> P1 with 3 modes of mean 1 but different distributions
     * (Exp, Erlang-2, HyperExp with SCV 4). Each mode renews at mean 1, so the
     * throughput is 3 exactly and is insensitive to the phase structure.
     */
    @Test
    public void testThreeMixedModesSelfLoopTputIsThree() {
        assertEquals(3.0, tputOf(TestSPNModels.createTestModelWithThreeModes()), TPUT_TOL,
                "marking-invariant SPN with 3 mean-1 modes must have Tput = 3");
    }

    /**
     * The same net with 3 identical Exp(1) modes: the race of 3 exponentials has
     * exit rate 3, so the throughput is 3 exactly.
     */
    @Test
    public void testThreeExpModesSelfLoopTputIsThree() {
        assertEquals(3.0, tputOf(TestSPNModels.createTestModelWithThreeExpModes()), TPUT_TOL,
                "marking-invariant SPN with 3 Exp(1) modes must have Tput = 3");
    }

    /**
     * Marking invariance itself: the single token never leaves P1, so QLen is 1 in
     * every state. This is the metric that stayed correct while Tput was wrong, and
     * it guards against "fixing" the throughput by perturbing the token accounting.
     */
    @Test
    public void testThreeMixedModesSelfLoopQLenIsOne() {
        Network model = TestSPNModels.createTestModelWithThreeModes();
        SolverCTMC solver = new SolverCTMC(model);
        NetworkAvgTable table = solver.getAvgTable();
        assertEquals(1.0, table.getQLen().get(0), TOL,
                "the token never leaves P1, so QLen must be 1");
    }
}
