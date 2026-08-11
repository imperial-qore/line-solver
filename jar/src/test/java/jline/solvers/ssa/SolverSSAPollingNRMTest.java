package jline.solvers.ssa;

import jline.VerboseLevel;
import jline.lang.Network;
import jline.lang.ClosedClass;
import jline.lang.constant.PollingType;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.processes.Exp;
import jline.solvers.ctmc.SolverCTMC;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Validates the NRM's POLLING scheduling discipline against the exact CTMC on a
 * closed two-class cyclic-polling model.
 * <p>
 * A single-server polling station cycles through the class buffers it serves;
 * the controller decides who is in service, and the discipline (EXHAUSTIVE,
 * GATED, K-LIMITED, DECREMENTING) decides when a visit ends. Switchover between
 * buffers may be immediate (folded, never a state) or a timed exponential leg
 * (a genuine dwell). The NRM carries the controller [mode,pos,swk,ctr] in a
 * per-node store and gates the single-server service reaction on it; this test
 * checks that its steady-state queue lengths and throughputs match the CTMC,
 * which is exact for these finite closed models.
 * </p>
 * <p>
 * The NRM must actually have run: Solver_ssa_analyzer downgrades an ineligible
 * model to the serial engine, and the serial engine is CORRECT for polling too,
 * so a silent fallback would make this test pass while never exercising the NRM
 * polling code it exists to cover. Hence the method assertion.
 * </p>
 */
public class SolverSSAPollingNRMTest {

    private static final int N1 = 2;
    private static final int N2 = 2;
    private static final int SAMPLES = 300000;
    private static final int[] SEEDS = {12345, 23000, 31337, 42000};
    /** Relative tolerance on the seed-averaged simulated means. */
    private static final double RTOL = 0.06;

    private static Network model(PollingType ptype, boolean sw, int kpar) {
        Network model = new Network("poll");
        Delay delay = new Delay(model, "Think");
        Queue q = new Queue(model, "Poll", SchedStrategy.POLLING);
        q.setNumberOfServers(1);
        ClosedClass c1 = new ClosedClass(model, "C1", N1, delay);
        ClosedClass c2 = new ClosedClass(model, "C2", N2, delay);
        delay.setService(c1, new Exp(0.5));
        delay.setService(c2, new Exp(0.7));
        q.setService(c1, new Exp(1.0));
        q.setService(c2, new Exp(1.3));
        if (ptype == PollingType.KLIMITED) {
            q.setPollingType(ptype, kpar);
        } else {
            q.setPollingType(ptype);
        }
        if (sw) {
            q.setSwitchover(c1, new Exp(2.0));
            q.setSwitchover(c2, new Exp(3.0));
        }
        model.link(model.serialRouting(delay, q, delay));
        return model;
    }

    /** Exact CTMC queue length (station x class) and throughput. */
    private static Matrix[] ctmc(PollingType ptype, boolean sw, int kpar) {
        SolverCTMC solver = new SolverCTMC(model(ptype, sw, kpar), "verbose", VerboseLevel.SILENT);
        Matrix q = solver.getAvgQLen();
        Matrix t = solver.getAvgTput();
        return new Matrix[]{q, t};
    }

    /** Seed-averaged NRM queue length and throughput; asserts NRM actually ran. */
    private static Matrix[] nrm(PollingType ptype, boolean sw, int kpar, String label) {
        Matrix qAcc = null;
        Matrix tAcc = null;
        for (int seed : SEEDS) {
            SolverSSA solver = new SolverSSA(model(ptype, sw, kpar));
            solver.options.method = "nrm";
            solver.options.samples = SAMPLES;
            solver.options.seed = seed;
            solver.options.verbose = VerboseLevel.SILENT;
            Matrix q = solver.getAvgQLen();
            Matrix t = solver.getAvgTput();
            assertTrue(solver.result != null && solver.result.method != null
                            && solver.result.method.contains("nrm"),
                    label + ": NRM did not run, method was "
                            + (solver.result == null ? "<no result>" : solver.result.method));
            if (qAcc == null) {
                qAcc = new Matrix(q);
                tAcc = new Matrix(t);
            } else {
                qAcc = qAcc.add(1.0, q);
                tAcc = tAcc.add(1.0, t);
            }
        }
        qAcc = qAcc.scale(1.0 / SEEDS.length);
        tAcc = tAcc.scale(1.0 / SEEDS.length);
        return new Matrix[]{qAcc, tAcc};
    }

    private static void check(PollingType ptype, boolean sw, int kpar, String label) {
        Matrix[] exact = ctmc(ptype, sw, kpar);
        Matrix[] sim = nrm(ptype, sw, kpar, label);
        int pollRow = 1;   // station 0 = Think (Delay), station 1 = Poll
        for (int r = 0; r < 2; r++) {
            double qe = exact[0].get(pollRow, r);
            double qs = sim[0].get(pollRow, r);
            double te = exact[1].get(pollRow, r);
            double ts = sim[1].get(pollRow, r);
            assertClose(qs, qe, label + " QLen[Poll,c" + (r + 1) + "]");
            assertClose(ts, te, label + " Tput[Poll,c" + (r + 1) + "]");
        }
    }

    private static void assertClose(double got, double exact, String what) {
        double err = Math.abs(got - exact) / Math.max(1e-6, Math.abs(exact));
        assertTrue(err < RTOL, what + ": NRM=" + got + " CTMC=" + exact
                + " (" + (100.0 * err) + "% off)");
    }

    @Test
    public void exhaustiveNoSwitchover() {
        check(PollingType.EXHAUSTIVE, false, 0, "EXHAUSTIVE-nosw");
    }

    @Test
    public void exhaustiveSwitchover() {
        check(PollingType.EXHAUSTIVE, true, 0, "EXHAUSTIVE-sw");
    }

    /**
     * GATED with only immediate (folded) switchover. The JAR CTMC is a broken
     * oracle here: its state-space enumeration for a gated controller whose
     * budget column exists but whose switchover is immediate returns an EMPTY
     * reachable set, so getAvgQLen/getAvgTput are all zeros (verified directly:
     * exhaustive+immediate works, gated+immediate returns 0, gated+timed
     * switchover works). This is a pre-existing SolverCTMC defect, not an NRM
     * one. The NRM is therefore validated against the MATLAB CTMC's exact queue
     * lengths for the same model (Poll station, class 1 = 1.06106, class 2 =
     * 1.17936), hard-coded so the assertion cannot lean on the buggy oracle.
     */
    @Test
    public void gatedNoSwitchover() {
        Matrix[] sim = nrm(PollingType.GATED, false, 0, "GATED-nosw");
        assertClose(sim[0].get(1, 0), 1.06106, "GATED-nosw QLen[Poll,c1] (MATLAB-exact)");
        assertClose(sim[0].get(1, 1), 1.17936, "GATED-nosw QLen[Poll,c2] (MATLAB-exact)");
    }

    @Test
    public void gatedSwitchover() {
        check(PollingType.GATED, true, 0, "GATED-sw");
    }

    @Test
    public void kLimitedSwitchover() {
        check(PollingType.KLIMITED, true, 2, "KLIMITED-sw");
    }

    @Test
    public void decrementingSwitchover() {
        check(PollingType.DECREMENTING, true, 0, "DECREMENTING-sw");
    }
}
