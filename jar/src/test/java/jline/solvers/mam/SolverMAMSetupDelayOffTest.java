package jline.solvers.mam;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Distribution;
import jline.lang.processes.Erlang;
import jline.lang.processes.Exp;
import org.junit.jupiter.api.Test;

import static jline.api.mam.Qbd_setupdelayoff.qbd_setupdelayoff;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Setup/delay-off support in MAM.
 *
 * <p>An open station races the delay-off against the aggregate interarrival
 * Exp(lambda), the idle a Poisson stream sees, and is exact. It must be
 * dispatched ahead of the exact shortcuts, which model the station setup-free
 * and would return the setup-free rho/(1-rho).
 *
 * <p>Closed setup is deliberately not pinned here: MAM's closed analysis is an
 * approximation (it reads 0.7105 on this model against an exact 0.9 from
 * MVA/NC/CTMC), and the codebases disagree on whether the setup branch is
 * reached at all for a plain closed queue. See _kb/log.md.
 *
 * <p>Reference values are MATLAB's, which agree with a 4e6-sample JMT
 * simulation of the same M/M/1 with setup (1.2290 and 3.4162).
 */
public class SolverMAMSetupDelayOffTest {

    private static final double LAMBDA = 0.5;
    private static final double MU = 1.0;
    private static final double BETA = 4.0;

    private static Network openModel(Distribution setup, Distribution delayoff) {
        Network model = new Network("open_setup");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Q", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass oclass = new OpenClass(model, "C", 0);
        source.setArrival(oclass, new Exp(LAMBDA));
        queue.setService(oclass, new Exp(MU));
        if (setup != null) {
            queue.setDelayOff(oclass, setup, delayoff);
        }
        model.link(model.serialRouting(source, queue, sink));
        return model;
    }

    private static Network closedModel(Distribution setup, Distribution delayoff) {
        Network model = new Network("closed_setup");
        Delay think = new Delay(model, "Think");
        Queue queue = new Queue(model, "Q", SchedStrategy.FCFS);
        ClosedClass cclass = new ClosedClass(model, "C", 3, think);
        think.setService(cclass, new Exp(1.0));
        queue.setService(cclass, new Exp(2.0));
        if (setup != null) {
            queue.setDelayOff(cclass, setup, delayoff);
        }
        model.link(model.serialRouting(think, queue));
        return model;
    }

    @Test
    public void test_qbdSetupDelayOff_matchesMatlab() {
        // The level summation must take exactly n phases per level: summing n+1
        // while advancing by n overlaps the levels and reads high, by 1.8% for a
        // fast setup and 25% for a slow one.
        assertEquals(1.22727272, qbd_setupdelayoff(LAMBDA, MU, 2.0, 1.0, BETA, 1.0), 1e-6);
        assertEquals(3.41379310, qbd_setupdelayoff(LAMBDA, MU, 0.2, 1.0, BETA, 1.0), 1e-6);
        assertEquals(1.21022727, qbd_setupdelayoff(LAMBDA, MU, 2.0, 0.25, BETA, 1.0), 1e-6);
    }

    @Test
    public void test_qbdSetupDelayOff_monotoneInSetupTime() {
        double previous = 0.0;
        double[] means = {0.1, 0.5, 1.0, 2.0, 5.0};
        for (int i = 0; i < means.length; i++) {
            double qn = qbd_setupdelayoff(LAMBDA, MU, 1.0 / means[i], 1.0, BETA, 1.0);
            assertTrue(qn > previous, "queue length must grow with the setup time");
            previous = qn;
        }
    }

    @Test
    public void test_snIsFunction_marksSetupStations() {
        jline.lang.NetworkStruct sn = closedModel(new Exp(2.0), new Exp(BETA)).getStruct(false);
        // Station-indexed: the Delay is not a function station, the Queue is.
        assertEquals(0.0, sn.isfunction.get(0, 0), 1e-12);
        assertEquals(1.0, sn.isfunction.get(1, 0), 1e-12);
    }

    @Test
    public void test_openSetup_matchesMatlab() {
        assertEquals(1.227273,
                new SolverMAM(openModel(new Exp(2.0), new Exp(BETA))).getAvgQLen().get(1, 0), 1e-5);
        assertEquals(3.413793,
                new SolverMAM(openModel(new Exp(0.2), new Exp(BETA))).getAvgQLen().get(1, 0), 1e-5);
    }

    @Test
    public void test_openSetup_isNotIgnored() {
        // Without setup the station is a plain M/M/1 holding rho/(1-rho) = 1 job.
        double plain = new SolverMAM(openModel(null, null)).getAvgQLen().get(1, 0);
        assertEquals(1.0, plain, 1e-6);
        assertTrue(new SolverMAM(openModel(new Exp(2.0), new Exp(BETA))).getAvgQLen().get(1, 0) > plain);
    }

    @Test
    public void test_openSetup_erlangIsScvSensitive() {
        // An Erlang setup has SCV 0.25, so it must not be read as exponential.
        double erlang = new SolverMAM(
                openModel(Erlang.fitMeanAndOrder(0.5, 4), new Exp(BETA))).getAvgQLen().get(1, 0);
        double exponential = new SolverMAM(
                openModel(new Exp(2.0), new Exp(BETA))).getAvgQLen().get(1, 0);
        assertTrue(Math.abs(erlang - exponential) > 1e-4, "setup SCV must affect the queue length");
    }

    @Test
    public void test_setupDelayOff_featureGate() {
        // Setup/delay-off is honoured only by MAM, JMT and LDES. Without a
        // featset entry, CTMC/SSA/MVA/NC accepted the model and solved it
        // setup-free, silently returning the plain M/M/1 answer.
        jline.lang.Network withSetup = openModel(new Exp(2.0), new Exp(BETA));
        assertTrue(withSetup.getUsedLangFeatures().inspectFeature("SetupDelayOff"));

        assertTrue(new SolverMAM(withSetup).supports(withSetup));
        assertFalse(new jline.solvers.ctmc.SolverCTMC(withSetup).supports(withSetup));
        assertFalse(new jline.solvers.ssa.SolverSSA(withSetup).supports(withSetup));
        assertFalse(new jline.solvers.mva.SolverMVA(withSetup).supports(withSetup));
        assertFalse(new jline.solvers.nc.SolverNC(withSetup).supports(withSetup));

        // The gate must reject setup, not the whole model class.
        jline.lang.Network plain = openModel(null, null);
        assertTrue(new jline.solvers.ctmc.SolverCTMC(plain).supports(plain));
    }

    @Test
    public void test_closedWithoutSetup_isUnaffected() {
        // Single-class closed Delay+Queue now routes to the exact ldqbd method
        // (default MAM), which returns the exact 0.9 (matching MVA/NC/CTMC), not
        // the old dec.source approximation 0.7105.
        assertEquals(0.9,
                new SolverMAM(closedModel(null, null)).getAvgRespT().get(1, 0), 1e-5);
    }
}
