package jline.solvers.ag;

import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.OpenSignal;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SignalType;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.lang.processes.Immediate;
import jline.solvers.NetworkAvgTable;
import org.junit.jupiter.api.Test;

import jline.lang.FeatureSet;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * G-network (Gelenbe network) regression for the RCAT methods of SolverMAM.
 *
 * Tandem Source -> Queue1 -> Queue2 -> Sink with a signal class whose routing
 * chain is Source -> Queue1 -> Queue2 -> Sink. The signal is annihilated at
 * Queue1, the first station it reaches, so Queue2 sees no removals and its
 * throughput equals its arrival rate. Expected values are the exact
 * (untruncated) ones, which SolverCTMC reproduces to five digits.
 */
public class RcatGnetworkTest {

    private static final double TOL = 1e-4;

    private Network build(double lambdaPos, double lambdaNeg, double mu1, double mu2,
                          SignalType signalType) {
        Network model = new Network("GNet");
        Source source = new Source(model, "Source");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");

        OpenClass posClass = new OpenClass(model, "Positive", 0);
        source.setArrival(posClass, new Exp(lambdaPos));
        queue1.setService(posClass, new Exp(mu1));
        queue2.setService(posClass, new Exp(mu2));

        OpenSignal negClass = new OpenSignal(model, "Negative", signalType);
        source.setArrival(negClass, new Exp(lambdaNeg));
        queue1.setService(negClass, new Exp(mu1));
        queue2.setService(negClass, new Exp(mu2));

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(posClass, posClass, source, queue1, 1.0);
        P.set(posClass, posClass, queue1, queue2, 1.0);
        P.set(posClass, posClass, queue2, sink, 1.0);
        P.set(negClass, negClass, source, queue1, 1.0);
        P.set(negClass, negClass, queue1, queue2, 1.0);
        P.set(negClass, negClass, queue2, sink, 1.0);
        model.link(P);
        return model;
    }

    /** Queue length of the positive class at station index ist (0 = Source). */
    private double qlen(NetworkAvgTable t, int ist) {
        return t.getQLen().get(2 * ist);
    }

    /** Throughput of the positive class at station index ist. */
    private double tput(NetworkAvgTable t, int ist) {
        return t.getTput().get(2 * ist);
    }

    /**
     * Negative customers: a single job is removed per signal, so each isolated
     * process stays birth-death and the standard INAP reversed-rate estimator
     * is exact. q1 = lambda/(mu1 + lambdaNeg), q2 = T1/mu2.
     */
    @Test
    public void testNegativeSignalTandem() {
        for (String method : new String[]{"inap", "inapplus", "inapinf"}) {
            SolverAG solver = new SolverAG(build(1.0, 0.3, 2.0, 3.0, SignalType.NEGATIVE), method);
            NetworkAvgTable t = solver.getAvgTable();
            assertEquals(0.769231, qlen(t, 1), TOL, method + ": Queue1 queue length");
            assertEquals(0.408163, qlen(t, 2), TOL, method + ": Queue2 queue length");
            assertEquals(0.869565, tput(t, 1), TOL, method + ": Queue1 throughput");
            assertEquals(tput(t, 1), tput(t, 2), TOL,
                    method + ": Queue2 sees no signal, so it must not lose jobs");
        }
    }

    /**
     * Catastrophes empty the station, so the isolated process is no longer
     * birth-death and the mean-of-ratios INAP estimator overestimates the
     * departure rate (it returned Queue2 throughput 1.5995 against an arrival
     * rate of 0.8). The rate-conservation estimator is used instead.
     */
    @Test
    public void testCatastropheSignalTandem() {
        for (String method : new String[]{"inap", "inapplus", "inapinf"}) {
            SolverAG solver = new SolverAG(build(1.0, 0.3, 2.0, 3.0, SignalType.CATASTROPHE), method);
            NetworkAvgTable t = solver.getAvgTable();
            assertEquals(0.666667, qlen(t, 1), TOL, method + ": Queue1 queue length");
            assertEquals(0.363636, qlen(t, 2), TOL, method + ": Queue2 queue length");
            assertEquals(0.8, tput(t, 1), TOL, method + ": Queue1 throughput");
            assertEquals(tput(t, 1), tput(t, 2), TOL,
                    method + ": Queue2 throughput cannot exceed its arrival rate");
        }
    }

    /**
     * Solver_ag_build is the only analyzer in LINE that reads sn.issignal, and it
     * belongs to SolverAG. A G-network reaching any MAM method would be solved
     * with every signal turned into an ordinary customer -- a different model,
     * answered without a word -- so every MAM method refuses it by name, and the
     * MAM envelope no longer declares the G-network names at all.
     */
    @Test
    public void testMamMethodsRefuseSignalsAndAgAcceptsThem() {
        String[] mamMethods = new String[]{"default", "dec.source", "dec.mmap", "dec.poisson",
                "mna", "ldqbd", "bgchain"};
        for (String method : mamMethods) {
            jline.solvers.mam.SolverMAM solver = new jline.solvers.mam.SolverMAM(
                    build(1.0, 0.3, 2.0, 3.0, SignalType.NEGATIVE), method);
            String reason = solver.supportsModelMethod(method);
            assertTrue(reason != null && reason.contains("does not support G-network signals"),
                    method + " must refuse a G-network by name, got: " + reason);
            FeatureSet fs = jline.solvers.mam.SolverMAM.methodFeatureSet(method);
            assertFalse(fs.inspectFeature("OpenSignal"), method + " must not declare OpenSignal");
            assertFalse(fs.inspectFeature("ClosedSignal"), method + " must not declare ClosedSignal");
            assertFalse(fs.inspectFeature("SignalType_NEGATIVE"),
                    method + " must not declare SignalType_NEGATIVE");
        }
        for (String method : new String[]{"inap", "inapplus", "inapinf"}) {
            SolverAG solver = new SolverAG(build(1.0, 0.3, 2.0, 3.0, SignalType.NEGATIVE), method);
            assertEquals("", solver.supportsModelMethod(method),
                    method + " must still accept a G-network");
            assertTrue(SolverAG.methodFeatureSet(method).inspectFeature("OpenSignal"),
                    method + " must declare OpenSignal");
        }
    }

    /**
     * A model with no signal class is untouched by the gate: the same methods
     * that refuse the G-network above must still accept an ordinary tandem.
     */
    @Test
    public void testSignalGateLeavesOrdinaryModelsAlone() {
        Network model = new Network("Tandem");
        Source source = new Source(model, "Source");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass posClass = new OpenClass(model, "Positive", 0);
        source.setArrival(posClass, new Exp(1.0));
        queue1.setService(posClass, new Exp(2.0));
        queue2.setService(posClass, new Exp(3.0));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(posClass, posClass, source, queue1, 1.0);
        P.set(posClass, posClass, queue1, queue2, 1.0);
        P.set(posClass, posClass, queue2, sink, 1.0);
        model.link(P);

        // The MAM methods that refuse the G-network above must still accept an
        // ordinary tandem: the gate is on the signal marking, not on the shape.
        for (String method : new String[]{"default", "dec.source"}) {
            jline.solvers.mam.SolverMAM solver = new jline.solvers.mam.SolverMAM(model, method);
            assertEquals("", solver.supportsModelMethod(method),
                    method + " must accept a signal-free model");
        }
        // ... and so must SolverAG, on its own methods.
        for (String method : new String[]{"default", "inap"}) {
            SolverAG solver = new SolverAG(model, method);
            assertEquals("", solver.supportsModelMethod(method),
                    method + " must accept a signal-free model");
        }
    }

    /**
     * A signal is a trigger, not a job: the AG builder skips it in processMap
     * and reads only its Source arrival rate, so its service process at a queue
     * station is irrelevant. Models routinely declare it Immediate, which the
     * exponential-only RCAT gate used to reject.
     */
    @Test
    public void testImmediateSignalServiceIsAccepted() {
        Network model = new Network("GNetImmediateSignal");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");

        OpenClass posClass = new OpenClass(model, "Positive", 0);
        source.setArrival(posClass, new Exp(1.0));
        queue.setService(posClass, new Exp(2.0));

        OpenSignal negClass = new OpenSignal(model, "Negative", SignalType.NEGATIVE);
        source.setArrival(negClass, new Exp(0.2));
        queue.setService(negClass, new Immediate());

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(posClass, posClass, source, queue, 1.0);
        P.set(posClass, posClass, queue, sink, 1.0);
        P.set(negClass, negClass, source, queue, 1.0);
        P.set(negClass, negClass, queue, sink, 1.0);
        model.link(P);

        SolverAG solver = new SolverAG(model, "inap");
        assertEquals("", solver.supportsModelMethod("inap"),
                "an Immediate signal service must not be rejected");
        NetworkAvgTable t = solver.getAvgTable();
        // rho = lambda+/(mu + lambda-) = 1/2.2, QLen = rho/(1-rho).
        double rho = 1.0 / 2.2;
        assertEquals(rho / (1.0 - rho), t.getQLen().get(2), TOL, "Queue queue length");
    }

    /** Same, with the catastrophe rate raised to the positive arrival rate. */
    @Test
    public void testCatastropheSignalStrong() {
        for (String method : new String[]{"inap", "inapplus", "inapinf"}) {
            SolverAG solver = new SolverAG(build(1.0, 1.0, 2.0, 3.0, SignalType.CATASTROPHE), method);
            NetworkAvgTable t = solver.getAvgTable();
            assertEquals(0.414214, qlen(t, 1), TOL, method + ": Queue1 queue length");
            assertEquals(0.242641, qlen(t, 2), TOL, method + ": Queue2 queue length");
            assertEquals(0.585786, tput(t, 1), TOL, method + ": Queue1 throughput");
            assertEquals(tput(t, 1), tput(t, 2), TOL,
                    method + ": Queue2 throughput cannot exceed its arrival rate");
        }
    }
}
