package jline.solvers.mam;

import jline.lang.*;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.*;
import jline.lang.processes.*;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

/**
 * RCAT (inap/inapplus/inapinf) models each station-class by its mean rate only,
 * so a non-exponential process must be rejected rather than silently answered as
 * if it were exponential.
 */
public class RcatGateTest {

    private Network build(Distribution svc) {
        Network model = new Network("m");
        Source source = new Source(model, "S");
        Queue queue = new Queue(model, "Q", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "K");
        OpenClass oc = new OpenClass(model, "C", 0);
        source.setArrival(oc, new Exp(0.5));
        queue.setService(oc, svc);
        model.link(Network.serialRouting(source, queue, sink));
        return model;
    }

    @Test
    public void testRcatAcceptsExponential() {
        for (String m : new String[]{"inap", "inapplus", "inapinf"}) {
            SolverMAM s = new SolverMAM(build(new Exp(1.0)), m);
            assertEquals("", s.supportsModelMethod(m), m + " must accept an exponential model");
        }
    }

    @Test
    public void testRcatRejectsNonExponential() {
        for (String m : new String[]{"inap", "inapplus", "inapinf"}) {
            SolverMAM s = new SolverMAM(build(new Erlang(3.0, 3)), m);
            String reason = s.supportsModelMethod(m);
            assertTrue(reason.contains("exponential processes only"),
                "RCAT " + m + " must reject Erlang service, got: " + reason);
            assertTrue(reason.contains("Erlang"), "reason must name the offending process: " + reason);

            SolverMAM s2 = new SolverMAM(build(new HyperExp(0.5, 3.0, 10.0)), m);
            assertTrue(s2.supportsModelMethod(m).contains("exponential processes only"),
                "RCAT " + m + " must reject HyperExp service");
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
        for (String m : new String[]{"inap", "inapplus", "inapinf"}) {
            Network model = new Network("mmc");
            Source source = new Source(model, "S");
            Queue queue = new Queue(model, "Q", SchedStrategy.FCFS);
            queue.setNumberOfServers(2);
            Sink sink = new Sink(model, "K");
            OpenClass oc = new OpenClass(model, "C", 0);
            source.setArrival(oc, new Exp(1.2));
            queue.setService(oc, new Exp(1.0));
            model.link(Network.serialRouting(source, queue, sink));

            SolverMAM s = new SolverMAM(model, m);
            String reason = s.supportsModelMethod(m);
            assertTrue(reason.contains("single-server stations only"),
                "RCAT " + m + " must reject a multiserver station, got: " + reason);
        }
    }
}
