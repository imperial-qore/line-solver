package jline.solvers.ssa;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.constant.DropStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.solvers.ctmc.SolverCTMC;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * A closed job that finds no room must BLOCK, and the NRM engine must say so.
 * <p>
 * Until 2026-08-19 it did not (BUG-81). {@code capacityLoss} correctly declined
 * to DROP a closed job -- a closed network's population is an invariant -- but
 * nothing then stopped the reaction, so the firing went into the full station
 * anyway and the engine returned the UNCONSTRAINED product-form answer: on the
 * model below, QLen [1.96 2.07 1.97] against the exact [3.609 0.971 1.420], a
 * mean of 2.07 at a station that holds 2. The serial engine, whose producer
 * (AfterEventStation) already implements the open/closed contract, was right the
 * whole time, so the two SSA methods disagreed with each other.
 * </p><p>
 * Population conservation alone does not catch this: the broken engine conserved
 * it too. The oracle is the exact CTMC, plus the arithmetic fact that a station
 * capped at 2 cannot hold 2.07 on average.
 * </p>
 */
public class SolverSSANrmClosedCapacityTest {

    private static final int SAMPLES = 200000;
    private static final int[] SEEDS = {23000, 24000, 25000, 26000};

    /** Closed 3-queue tandem, Exp(1) FCFS everywhere, Q2 capped at cap. */
    private static Network cappedTandem(int cap, int njobs) {
        Network model = new Network("tandem");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.FCFS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.FCFS);
        Queue q3 = new Queue(model, "Q3", SchedStrategy.FCFS);
        ClosedClass c1 = new ClosedClass(model, "C1", njobs, q1, 0);
        q1.setService(c1, new Exp(1.0));
        q2.setService(c1, new Exp(1.0));
        q3.setService(c1, new Exp(1.0));
        q2.setClassCapacity(c1, cap);
        model.link(model.serialRouting(q1, q2, q3));
        return model;
    }

    /** M/M/1/K, lambda 0.8, mu 1, explicit DROP -- the open half of the contract. */
    private static Network mm1k(int K) {
        Network model = new Network("mm1k");
        Source source = new Source(model, "Source");
        Queue q = new Queue(model, "Q", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass c1 = new OpenClass(model, "C1");
        source.setArrival(c1, new Exp(0.8));
        q.setService(c1, new Exp(1.0));
        q.setClassCapacity(c1, K);
        q.setDropRule(c1, DropStrategy.Drop);
        model.link(model.serialRouting(source, q, sink));
        return model;
    }

    /** Mean over SEEDS of the per-station queue lengths under the NRM engine. */
    private static double[] nrmQLen(Network model, int nstations) {
        double[] acc = new double[nstations];
        for (int seed : SEEDS) {
            SolverSSA solver = new SolverSSA(model);
            solver.options.method = "nrm";
            solver.options.samples = SAMPLES;
            solver.options.seed = seed;
            Matrix q = solver.getAvgQLen();
            // The NRM must actually have run: Solver_ssa_analyzer downgrades an
            // ineligible model to the serial engine, which is CORRECT here, so a
            // silent fallback would pass this test without exercising the fix.
            assertTrue(solver.result != null && solver.result.method != null
                            && solver.result.method.contains("nrm"),
                    "NRM did not run: method was "
                            + (solver.result == null ? "<no result>" : solver.result.method));
            for (int i = 0; i < nstations; i++) {
                acc[i] += q.get(i, 0);
            }
        }
        for (int i = 0; i < nstations; i++) {
            acc[i] /= SEEDS.length;
        }
        return acc;
    }

    /** The capped station's mean cannot exceed its cap, and did (2.07 > 2). */
    @Test
    public void testNrmRespectsBindingClosedClassCapacity() {
        double[] q = nrmQLen(cappedTandem(2, 6), 3);
        assertTrue(q[1] <= 2.0,
                "Q2 holds at most 2 jobs, so its mean cannot be " + q[1]);
    }

    /** Against the exact chain: the constrained answer, not the free one. */
    @Test
    public void testNrmMatchesExactChainOnCappedClosedModel() {
        Network model = cappedTandem(2, 6);
        Matrix exact = new SolverCTMC(model).getAvgQLen();
        double[] q = nrmQLen(model, 3);
        double pop = 0.0;
        for (int i = 0; i < 3; i++) {
            double err = Math.abs(q[i] - exact.get(i, 0));
            assertTrue(err < 0.15, "station " + i + ": NRM " + q[i]
                    + ", exact " + exact.get(i, 0));
            pop += q[i];
        }
        // necessary, not sufficient -- the broken engine conserved it too
        assertTrue(Math.abs(pop - 6.0) < 1e-9, "population drifted to " + pop);
    }

    /**
     * The regression guard: an OPEN arrival that finds no room is still LOST.
     * The blocking gate is keyed on the class type, so M/M/1/K must keep the
     * exact loss behaviour it had; if it starts blocking instead, the source
     * stops offering and these means move. Station row 1 is the Queue.
     */
    @Test
    public void testOpenClassDropIsUnchanged() {
        double[] q1 = nrmQLen(mm1k(1), 2);
        assertTrue(Math.abs(q1[1] - 4.0 / 9.0) < 0.02,
                "M/M/1/1 QLen " + q1[1] + ", exact " + (4.0 / 9.0));
        double[] q2 = nrmQLen(mm1k(2), 2);
        assertTrue(Math.abs(q2[1] - 0.852459016393443) < 0.02,
                "M/M/1/2 QLen " + q2[1] + ", exact 0.852459016393443");
    }
}
