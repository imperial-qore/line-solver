package jline.examples.basic;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.examples.java.basic.StochPetriNetModel;
import jline.lang.Network;
import jline.solvers.NetworkAvgTable;
import jline.solvers.wrappers.jmt.SolverJMT;
import jline.solvers.ldes.LDESOptions;
import jline.solvers.ldes.SolverLDES;
import jline.util.Maths;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Validates the LDES simulator against SolverJMT on the stochastic Petri net
 * examples. Both engines simulate the same generalized stochastic Petri net, so
 * with a common seed/sample budget their per-place mean metrics must agree within
 * simulation tolerance. These tests exercise the atomic-firing race semantics of
 * the LDES Place/Transition engine (competing transitions, redundant firing modes,
 * inhibiting arcs, immediate transitions and non-Markovian firing times), which the
 * ground-truth SolverJMT/SolverCTMC realise as the underlying GSPN Markov chain.
 *
 * <p>The heavy-tailed {@code spn_pareto_service} model and the small-magnitude
 * places of the immediate-transition {@code spn_open_sevenplaces} model converge
 * slowly, so those two examples are validated by the dedicated harness rather than
 * here to keep the suite deterministic; all remaining examples agree at the default
 * sample budget.
 */
public class StochPetriNetLDESTest {

    private static final int SEED = 23000;
    private static final int SAMPLES = 200000;
    /** Relative tolerance for two independent simulators on non-trivial metrics. */
    private static final double REL_TOL = 5e-2;
    /** Metrics below this magnitude carry too much relative noise to compare. */
    private static final double ABS_FLOOR = 5e-2;

    @BeforeAll
    public static void setUp() {
        Maths.setRandomNumbersMatlab(true);
        GlobalConstants.setVerbose(VerboseLevel.SILENT);
    }

    private void validateLDESAgainstJMT(String name, Network jmtModel, Network ldesModel) {
        final NetworkAvgTable[] jmtHolder = new NetworkAvgTable[1];
        jline.TestTools.withSuppressedOutput(() -> {
            SolverJMT solver = new SolverJMT(jmtModel, "seed", SEED, "samples", SAMPLES,
                    "keep", 2, "verbose", 0);
            jmtHolder[0] = solver.getAvgTable();
        });
        NetworkAvgTable jmt = jmtHolder[0];

        LDESOptions opt = new LDESOptions();
        opt.verbose = VerboseLevel.SILENT;
        opt.seed = SEED;
        opt.samples = SAMPLES;
        NetworkAvgTable ldes = new SolverLDES(ldesModel, opt).getAvgTable();

        compare(name, "QLen", jmt.getQLen(), ldes.getQLen());
        compare(name, "Tput", jmt.getTput(), ldes.getTput());
        compare(name, "RespT", jmt.getRespT(), ldes.getRespT());
    }

    private void compare(String name, String metric, List<Double> jmt, List<Double> ldes) {
        assertEquals(jmt.size(), ldes.size(),
                name + ": " + metric + " row count mismatch (JMT " + jmt.size()
                        + " vs LDES " + ldes.size() + ")");
        for (int i = 0; i < jmt.size(); i++) {
            double ref = jmt.get(i);
            double got = ldes.get(i);
            if (Math.abs(ref) <= ABS_FLOOR) {
                continue;
            }
            double relErr = Math.abs(got - ref) / Math.abs(ref);
            assertTrue(relErr <= REL_TOL,
                    String.format("%s: %s[%d] LDES=%.6f deviates from JMT=%.6f by %.1f%% (> %.0f%%)",
                            name, metric, i, got, ref, relErr * 100, REL_TOL * 100));
        }
    }

    @Test
    public void testSpnBasicClosedLDESvsJMT() {
        validateLDESAgainstJMT("spn_basic_closed",
                StochPetriNetModel.spn_basic_closed(), StochPetriNetModel.spn_basic_closed());
    }

    @Test
    public void testSpnBasicOpenLDESvsJMT() {
        validateLDESAgainstJMT("spn_basic_open",
                StochPetriNetModel.spn_basic_open(), StochPetriNetModel.spn_basic_open());
    }

    @Test
    public void testSpnTwoModesLDESvsJMT() {
        validateLDESAgainstJMT("spn_twomodes",
                StochPetriNetModel.spn_twomodes(), StochPetriNetModel.spn_twomodes());
    }

    @Test
    public void testSpnFourModesLDESvsJMT() {
        validateLDESAgainstJMT("spn_fourmodes",
                StochPetriNetModel.spn_fourmodes(), StochPetriNetModel.spn_fourmodes());
    }

    @Test
    public void testSpnInhibitingLDESvsJMT() {
        validateLDESAgainstJMT("spn_inhibiting",
                StochPetriNetModel.spn_inhibiting(), StochPetriNetModel.spn_inhibiting());
    }

    @Test
    public void testSpnClosedTwoPlacesLDESvsJMT() {
        validateLDESAgainstJMT("spn_closed_twoplaces",
                StochPetriNetModel.spn_closed_twoplaces(), StochPetriNetModel.spn_closed_twoplaces());
    }

    @Test
    public void testSpnClosedFourPlacesLDESvsJMT() {
        validateLDESAgainstJMT("spn_closed_fourplaces",
                StochPetriNetModel.spn_closed_fourplaces(), StochPetriNetModel.spn_closed_fourplaces());
    }
}
