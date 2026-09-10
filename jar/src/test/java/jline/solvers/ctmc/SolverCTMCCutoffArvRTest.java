package jline.solvers.ctmc;

import jline.VerboseLevel;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Offered-vs-carried arrival rate (ArvR) at a bound in SolverCTMC.
 *
 * A finite bound on an open class is EITHER a physical capacity (setCapacity /
 * setClassCapacity) OR a state-space cutoff imposed only to keep the CTMC
 * enumeration finite. They must report ArvR differently, and the discriminator is
 * the station drop rule (State.isPhysicalCapacity / State.arrivalIsLost, applied
 * by AfterEventStation, which the CTMC state-space generator calls):
 *
 *   PHYSICAL cap  -> a job meeting a full buffer is really lost, so the loss is a
 *                    real event: ArvR reports the OFFERED rate (lambda), Tput < lambda.
 *   state-space CUTOFF -> the refused job never existed; a truncation artifact
 *                    that counts nowhere: ArvR reports the CARRIED rate (== Tput).
 *
 * This is the USER DECISION of 2026-07-17 (a cutoff is not a physical capacity),
 * documented in _kb/06-solver-catalog.md, and mirrors the SSA convention (BUG-85:
 * M/M/1/1 must report ArvR 0.800000, not the carried 0.444444).
 *
 * The BINDING control is that K=2 (physical) and cutoff=2 (truncation) yield an
 * IDENTICAL stationary distribution (M/M/1/2 == M/M/1 truncated at 2, both Tput
 * 0.630996) yet must report DIFFERENT ArvR: 0.9 (cap) vs 0.630996 (cutoff). A
 * solver that keyed ArvR off the state space alone -- or folded the cutoff into a
 * physical cap -- returns the same number for both and fails here.
 */
public class SolverCTMCCutoffArvRTest {

    private static final double LAM = 0.9;
    private static final double MU = 1.0;
    private static final double TOL = 1e-9;

    /** Source(Exp lambda) -> Queue FCFS(Exp mu) -> Sink; optional physical cap. */
    private static Network buildMM1(int capacity) {
        Network model = new Network("mm1");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass cls = new OpenClass(model, "C1");
        source.setArrival(cls, new Exp(LAM));
        queue.setService(cls, new Exp(MU));
        if (capacity > 0) {
            queue.setCapacity(capacity);
        }
        model.link(model.serialRouting(source, queue, sink));
        return model;
    }

    /** Returns {ArvR, Tput} at the Queue (station row 1, class 0). */
    private static double[] solve(int capacity, int cutoff) {
        SolverCTMC solver = new SolverCTMC(buildMM1(capacity), "cutoff", cutoff,
                "verbose", VerboseLevel.SILENT);
        Matrix AN = solver.getAvgArvR();
        Matrix TN = solver.getAvgTput();
        return new double[]{AN.get(1, 0), TN.get(1, 0)};
    }

    @Test
    public void cutoffReportsCarriedArvR() {
        // infinite capacity + tight cutoff: ArvR == carried Tput, both < lambda.
        int[] cutoffs = {2, 3, 5, 8};
        for (int cutoff : cutoffs) {
            double[] r = solve(-1, cutoff);
            assertEquals(r[1], r[0], TOL, "cutoff=" + cutoff + ": ArvR must equal carried Tput");
            assertTrue(r[0] < LAM - 1e-6, "cutoff=" + cutoff + ": carried ArvR must be below offered lambda");
        }
    }

    @Test
    public void physicalCapReportsOfferedArvR() {
        // finite physical capacity K + loose cutoff: ArvR == offered lambda, Tput < lambda.
        int[] Ks = {2, 3, 5};
        for (int K : Ks) {
            double[] r = solve(K, Math.max(K + 5, 20));
            assertEquals(LAM, r[0], TOL, "K=" + K + ": physical cap must report offered lambda");
            assertTrue(r[1] < LAM - 1e-6, "K=" + K + ": carried Tput must be below offered lambda");
        }
    }

    @Test
    public void cutoffAndCapShareDistributionButDifferInArvR() {
        // Decisive control: K=2 (physical) and cutoff=2 (truncation) share the same
        // stationary throughput but must report DIFFERENT ArvR.
        double[] cut = solve(-1, 2);
        double[] cap = solve(2, 20);
        assertEquals(cut[1], cap[1], TOL, "identical carried throughput");
        assertEquals(cut[1], cut[0], TOL, "cutoff -> carried ArvR");
        assertEquals(LAM, cap[0], TOL, "physical cap -> offered ArvR");
        assertTrue(cap[0] > cut[0] + 1e-3, "the two conventions must differ");
    }
}
