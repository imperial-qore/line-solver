package jline.solvers.ssa;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.processes.Coxian;
import jline.lang.processes.Distribution;
import jline.lang.processes.Erlang;
import jline.lang.processes.Exp;
import jline.lang.processes.HyperExp;
import jline.solvers.ctmc.SolverCTMC;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Validates the NRM phase expansion at the non-preemptive BUFFERED family
 * (FCFS/LCFS/SIRO/HOL/SEPT/LEPT) against the exact CTMC.
 * <p>
 * Unlike the INF/PS family (see {@link SolverSSANrmPhaseTest}), a FCFS station is
 * NOT insensitive: its queue length depends on the service-time variability, so
 * the exact value is distribution-specific and is taken here from the CTMC. Only
 * the jobs actually in service carry a phase; the waiting jobs sit in the buffer
 * without one. The engine tracks the in-service phase multiset in the auxiliary
 * structure svcph, leaving nvec (the class total) unchanged, and this test is the
 * one that exercises that structure -- an all-exponential buffered model is blind
 * to it (nph == 1 makes svcph a single slot equal to the in-service count).
 * </p>
 * <p>
 * The sensitivity is itself the binding control: at a FCFS station a lower-SCV
 * service (Erlang) must give a SHORTER queue than exponential, and a higher-SCV
 * service (HyperExp/Coxian) a LONGER one. An implementation that silently
 * collapsed the distribution onto its mean, or fell back to a product-form
 * approximation, would return the exponential value for all of them and fail the
 * ordering assertion even if it passed a loose per-point band.
 * </p>
 */
public class SolverSSANrmBufferedPhaseTest {

    private static final int SAMPLES = 400000;
    private static final int SEED = 23000;
    private static final double RTOL = 0.02;   // several sd above the noise floor at SAMPLES
    private static final double MEAN = 0.5;
    private static final int NJOBS = 3;

    private static Distribution[] dists() {
        return new Distribution[]{
                new Exp(1.0 / MEAN),
                Erlang.fitMeanAndOrder(MEAN, 2),
                Erlang.fitMeanAndOrder(MEAN, 3),
                HyperExp.fitMeanAndSCV(MEAN, 4.0),
                Coxian.fitMeanAndSCV(MEAN, 3.0),
        };
    }

    private static String[] labels() {
        return new String[]{"Exp", "Erlang-2", "Erlang-3", "HyperExp-4", "Coxian-3"};
    }

    private static Network model(SchedStrategy sched, int nservers, Distribution svc) {
        Network model = new Network("bufph_" + sched);
        Queue queue = new Queue(model, "Q1", sched);
        queue.setNumberOfServers(nservers);
        Delay delay = new Delay(model, "Delay");
        ClosedClass c1 = new ClosedClass(model, "C1", NJOBS, delay);
        queue.setService(c1, svc);
        delay.setService(c1, new Exp(1.0));
        model.link(model.serialRouting(delay, queue));
        return model;
    }

    // The Queue node is constructed before the Delay, so it is station index 0 in
    // the getAvgQLen table (station order follows node construction, not routing).
    private static final int QUEUE_STATION = 0;

    private static double nrmQ(Network model) {
        SolverSSA solver = new SolverSSA(model);
        solver.options.method = "nrm";
        solver.options.samples = SAMPLES;
        solver.options.seed = SEED;
        return solver.getAvgQLen().get(QUEUE_STATION, 0);
    }

    private static double ctmcQ(Network model) {
        // Closed model with a finite population, so the CTMC state space is finite
        // and no cutoff is needed.
        SolverCTMC solver = new SolverCTMC(model);
        return solver.getAvgQLen().get(QUEUE_STATION, 0);
    }

    private void checkAgainstCtmc(SchedStrategy sched, int nservers) {
        Distribution[] ds = dists();
        String[] names = labels();
        double[] q = new double[ds.length];
        for (int i = 0; i < ds.length; i++) {
            double got = nrmQ(model(sched, nservers, ds[i]));
            double exact = ctmcQ(model(sched, nservers, ds[i]));
            q[i] = got;
            double err = Math.abs(got - exact) / exact;
            assertTrue(err < RTOL, sched + "(" + nservers + " srv)/" + names[i]
                    + ": NRM QLen " + got + " deviates from CTMC " + exact
                    + " by " + (100.0 * err) + "%");
        }
        if (nservers == 1 && sched == SchedStrategy.FCFS) {
            // Binding control: FCFS is variability-sensitive. Erlang-2/3 (SCV<1)
            // must be shorter than Exp, HyperExp/Coxian (SCV>1) longer. A
            // mean-collapsed or product-form impl would return the Exp value for
            // all five and fail here.
            assertTrue(q[1] < q[0] && q[2] < q[0],
                    "FCFS: Erlang queue should be shorter than Exp (got Erlang2=" + q[1]
                            + ", Erlang3=" + q[2] + ", Exp=" + q[0] + ")");
            assertTrue(q[3] > q[0] && q[4] > q[0],
                    "FCFS: HyperExp/Coxian queue should be longer than Exp (got HyperExp=" + q[3]
                            + ", Coxian=" + q[4] + ", Exp=" + q[0] + ")");
        }
    }

    @Test
    public void testFcfsPhaseSingleServer() {
        checkAgainstCtmc(SchedStrategy.FCFS, 1);
    }

    @Test
    public void testFcfsPhaseMultiServer() {
        checkAgainstCtmc(SchedStrategy.FCFS, 2);
    }

    @Test
    public void testHolPhaseSingleServer() {
        checkAgainstCtmc(SchedStrategy.HOL, 1);
    }

    @Test
    public void testSiroPhaseSingleServer() {
        checkAgainstCtmc(SchedStrategy.SIRO, 1);
    }
}
