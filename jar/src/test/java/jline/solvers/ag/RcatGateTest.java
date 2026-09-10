package jline.solvers.ag;

import jline.lang.*;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.*;
import jline.lang.processes.*;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

/**
 * RCAT (inap/inapplus/inapinf) gives every component a service-phase and an
 * arrival-phase dimension, so any law with a genuine (D0,D1) Markovian
 * representation is admissible. What is not must be rejected rather than
 * silently answered as if it were exponential.
 */
public class RcatGateTest {

    private static final String[] METHODS = {"inap", "inapplus", "inapinf"};

    private Network build(Distribution svc) {
        return build(new Exp(0.5), svc);
    }

    private Network build(Distribution arr, Distribution svc) {
        Network model = new Network("m");
        Source source = new Source(model, "S");
        Queue queue = new Queue(model, "Q", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "K");
        OpenClass oc = new OpenClass(model, "C", 0);
        source.setArrival(oc, arr);
        queue.setService(oc, svc);
        model.link(Network.serialRouting(source, queue, sink));
        return model;
    }

    @Test
    public void testRcatAcceptsExponential() {
        for (String m : METHODS) {
            SolverAG s = new SolverAG(build(new Exp(1.0)), m);
            assertEquals("", s.supportsModelMethod(m), m + " must accept an exponential model");
        }
    }

    /**
     * The service-phase dimension is what makes these admissible: a component is
     * a QBD over (queue length, phase), so an Erlang, a HyperExp or a Coxian is
     * represented exactly rather than collapsed to its mean rate.
     */
    @Test
    public void testRcatAcceptsPhaseTypeService() {
        for (String m : METHODS) {
            assertEquals("", new SolverAG(build(new Erlang(3.0, 3)), m).supportsModelMethod(m),
                "RCAT " + m + " must accept Erlang service");
            assertEquals("", new SolverAG(build(new HyperExp(0.5, 3.0, 10.0)), m).supportsModelMethod(m),
                "RCAT " + m + " must accept HyperExp service");
            assertEquals("", new SolverAG(build(new Cox2(1.0, 2.0, 0.5)), m).supportsModelMethod(m),
                "RCAT " + m + " must accept Coxian service");
            // Non-Markovian laws are admissible because sn_nonmarkov_toph fits
            // them to a phase-type before the analyzer sees them.
            assertEquals("", new SolverAG(build(new Det(1.0)), m).supportsModelMethod(m),
                "RCAT " + m + " must accept Det service");
        }
    }

    /** The arrival-phase dimension does the same for a non-Poisson Source. */
    @Test
    public void testRcatAcceptsPhaseTypeArrivals() {
        for (String m : METHODS) {
            SolverAG s = new SolverAG(build(new Erlang(1.0, 2), new Exp(1.0)), m);
            assertEquals("", s.supportsModelMethod(m), "RCAT " + m + " must accept Erlang arrivals");
        }
    }

    /**
     * A matrix-exponential is not a generator: its off-diagonal entries are not
     * rates, so a CTMC assembled from it is a rational generator whose
     * stationary solution is a signed vector. Refuse rather than mis-answer.
     */
    @Test
    public void testRcatRejectsMatrixExponential() {
        Matrix alpha = new Matrix(1, 2);
        alpha.set(0, 0, 1.0);
        alpha.set(0, 1, 0.0);
        Matrix A = new Matrix(2, 2);
        A.set(0, 0, -1.0);
        A.set(0, 1, 2.0);
        A.set(1, 0, 0.0);
        A.set(1, 1, -3.0);
        for (String m : METHODS) {
            String reason = new SolverAG(build(new ME(alpha, A)), m).supportsModelMethod(m);
            assertTrue(reason.contains("Markovian"),
                "RCAT " + m + " must reject an ME service law, got: " + reason);
        }
    }

    /**
     * build_rcat never reads sn.nservers, so a multiserver station is driven at
     * rho = lambda/mu instead of lambda/(c*mu). A stable M/M/2 (rho=0.6) then
     * looks unstable in isolation and pins at the maxStates truncation: 94.0
     * instead of the exact 1.875. Reject rather than mis-answer.
     */
    @Test
    public void testRcatRejectsMultiserver() {
        for (String m : METHODS) {
            Network model = new Network("mmc");
            Source source = new Source(model, "S");
            Queue queue = new Queue(model, "Q", SchedStrategy.FCFS);
            queue.setNumberOfServers(2);
            Sink sink = new Sink(model, "K");
            OpenClass oc = new OpenClass(model, "C", 0);
            source.setArrival(oc, new Exp(1.2));
            queue.setService(oc, new Exp(1.0));
            model.link(Network.serialRouting(source, queue, sink));

            SolverAG s = new SolverAG(model, m);
            String reason = s.supportsModelMethod(m);
            assertTrue(reason.contains("single-server stations only"),
                "RCAT " + m + " must reject a multiserver station, got: " + reason);
        }
    }
}
