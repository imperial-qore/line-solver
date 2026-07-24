package jline.solvers.ctmc;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.io.Ret.ProbabilityResult;
import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.nodes.StatefulNode;
import jline.lang.processes.Exp;
import jline.solvers.Solver;
import jline.solvers.SolverOptions;
import jline.util.Maths;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;

/**
 * Value assertions for the CTMC state-probability getters.
 *
 * <p>These exist because {@code SolverCoverageTest} only asserts that the
 * probability getters return non-null, and {@code getProb} was not called at
 * all. Under that coverage {@code Solver_ctmc_marg} carried three stacked bugs
 * -- a state query matched against the full state space instead of the
 * station's columns, a state index read as a matrix row index, and an
 * exclusive-bound sum that normalised by zero -- each masking the next, so the
 * method had never returned a correct number. Calling a getter proves nothing;
 * only its value does.
 *
 * <p>Reference values are MATLAB, the ground-truth implementation:
 * {@code SolverCTMC.getProb} on this model returns 0.5294117647 for both
 * stations, and python-native agrees to the last digit printed.
 */
public class SolverCTMCProbTest {

    /** Probability of the initial state (both jobs at the Delay), from MATLAB. */
    private static final double P_INITIAL_STATE = 0.5294117647058824;
    private static final double TOL = 1e-10;

    @BeforeAll
    public static void setUp() {
        Maths.setRandomNumbersMatlab(true);
        GlobalConstants.setVerbose(VerboseLevel.SILENT);
    }

    /**
     * Two-station closed model, one class, two jobs: Delay Exp(1) -> FCFS
     * Queue Exp(3) -> Delay. Its chain has three states and the stationary
     * distribution is (1, 3, 9)/13 over (0,2), (1,1), (2,0) jobs at the Queue,
     * so the state both jobs sit at the Delay carries 9/17 once the FCFS phase
     * columns are accounted for. The model is small enough that the expected
     * value is checkable by hand and large enough that a wrong column slice
     * matches nothing.
     */
    private static Network closedDelayQueue() {
        Network model = new Network("cqn_delay_queue");
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue1", SchedStrategy.FCFS);
        ClosedClass cc = new ClosedClass(model, "Class1", 2, delay);
        delay.setService(cc, new Exp(1.0));
        queue.setService(cc, new Exp(3.0));
        model.link(model.serialRouting(delay, queue));
        return model;
    }

    private static SolverOptions ctmcOptions() {
        SolverOptions options = Solver.defaultOptions();
        options.verbose = VerboseLevel.SILENT;
        options.cutoff = Matrix.singleton(3);
        options.seed = 1;
        return options;
    }

    @Test
    public void getProbReturnsTheStationaryProbabilityOfTheSetState() {
        assertEquals(2, closedDelayQueue().getStatefulNodes().size(),
                "expected the Delay and the Queue to be stateful");

        // One solver per query on purpose: getProb rewrites the states of every
        // OTHER stateful node to -1 on the struct it shares with the model, so
        // querying two nodes from one solver would measure that side effect as
        // much as the getter.
        for (int i = 0; i < 2; i++) {
            Network model = closedDelayQueue();
            SolverCTMC solver = new SolverCTMC(model, ctmcOptions());
            solver.getAvg();
            StatefulNode node = model.getStatefulNodes().get(i);
            ProbabilityResult prob = solver.getProb(node);
            assertNotNull(prob, "getProb returned null for stateful node " + i);
            assertEquals(P_INITIAL_STATE, prob.getScalarProbability(), TOL,
                    "getProb(" + node.getName() + ") disagrees with the MATLAB reference");
        }
    }

    @Test
    public void getProbAggrAgreesWithGetProbOnAnExponentialModel() {
        // With exponential service there is one service phase, so the detailed
        // state and the per-class job counts describe the same event and the
        // two getters must return the same number. They are separate code paths
        // (Solver_ctmc_marg and Solver_ctmc_margaggr), which is what makes the
        // agreement worth asserting.
        //
        // Each getter gets its OWN solver, for the reason given above.
        Network m1 = closedDelayQueue();
        SolverCTMC s1 = new SolverCTMC(m1, ctmcOptions());
        s1.getAvg();
        double marg = s1.getProb(m1.getStatefulNodes().get(1)).getScalarProbability();

        Network m2 = closedDelayQueue();
        SolverCTMC s2 = new SolverCTMC(m2, ctmcOptions());
        s2.getAvg();
        double margAggr = s2.getProbAggr(m2.getStatefulNodes().get(1)).getScalarProbability();

        assertEquals(P_INITIAL_STATE, margAggr, TOL,
                "getProbAggr disagrees with the MATLAB reference");
        assertEquals(marg, margAggr, TOL,
                "getProb and getProbAggr disagree on an exponential model");
    }

    @Test
    public void getProbSysReturnsTheSameStateProbability() {
        // The system state is the join over stations, and every station is in
        // its initial state, so this is the same probability again.
        Network model = closedDelayQueue();
        SolverCTMC solver = new SolverCTMC(model, ctmcOptions());
        solver.getAvg();
        assertEquals(P_INITIAL_STATE, solver.getProbSys().getScalarProbability(), TOL,
                "getProbSys disagrees with the MATLAB reference");

        Network model2 = closedDelayQueue();
        SolverCTMC solver2 = new SolverCTMC(model2, ctmcOptions());
        solver2.getAvg();
        assertEquals(P_INITIAL_STATE, solver2.getProbSysAggr().getScalarProbability(), TOL,
                "getProbSysAggr disagrees with the MATLAB reference");
    }
}
