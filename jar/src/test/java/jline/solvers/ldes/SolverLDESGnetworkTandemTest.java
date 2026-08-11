package jline.solvers.ldes;

import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.Signal;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SignalType;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.solvers.NetworkAvgTable;
import org.junit.jupiter.api.Test;

import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Tandem G-network sample-path semantics for the LDES engine, checked against
 * the exact (SolverCTMC / SolverMAM) values.
 *
 * Two invariants are covered:
 *
 *  - A removal signal is annihilated at the station it reaches, so a signal
 *    routed Source -> Queue1 -> Queue2 -> Sink fires at Queue1 only. Routing it
 *    onward made a single signal remove one job per downstream station, which
 *    showed up as a Queue2 throughput deficit of lambdaNeg * P(Queue2 busy).
 *
 *  - A catastrophe empties the station, in-service job included. Tracking of
 *    in-service jobs used to be gated on the presence of a NEGATIVE signal, so
 *    a catastrophe-only model removed waiting jobs only and under-removed.
 */
public class SolverLDESGnetworkTandemTest extends SolverLDESTestFixtures {

    private Network build(double lambdaPos, double lambdaNeg, double mu1, double mu2,
                          SignalType signalType) {
        Network model = new Network("GNet");
        Source source = new Source(model, "Source");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");

        OpenClass posClass = new OpenClass(model, "Positive");
        source.setArrival(posClass, new Exp(lambdaPos));
        queue1.setService(posClass, new Exp(mu1));
        queue2.setService(posClass, new Exp(mu2));

        Signal negClass = new Signal(model, "Negative", signalType);
        source.setArrival(negClass, new Exp(lambdaNeg));
        queue1.setService(negClass, new Exp(mu1));
        queue2.setService(negClass, new Exp(mu2));

        RoutingMatrix routing = model.initRoutingMatrix();
        routing.set(posClass, posClass, source, queue1, 1.0);
        routing.set(posClass, posClass, queue1, queue2, 1.0);
        routing.set(posClass, posClass, queue2, sink, 1.0);
        routing.set(negClass, negClass, source, queue1, 1.0);
        routing.set(negClass, negClass, queue1, queue2, 1.0);
        routing.set(negClass, negClass, queue2, sink, 1.0);
        model.link(routing);
        return model;
    }

    private double metric(NetworkAvgTable t, List<Double> values, String station) {
        List<String> stations = t.getStationNames();
        List<String> classes = t.getClassNames();
        for (int i = 0; i < stations.size(); i++) {
            if (stations.get(i).equals(station) && classes.get(i).equals("Positive")) {
                return values.get(i);
            }
        }
        throw new IllegalStateException("no Positive row for station " + station);
    }

    private NetworkAvgTable solve(Network model) {
        LDESOptions options = createTestOptions(500000);
        return new SolverLDES(model, options).getAvgTable();
    }

    /**
     * Negative signal, lambdaPos = 1, lambdaNeg = 0.3, mu1 = 2, mu2 = 3.
     * Exact: Q1 = 0.769231, Q2 = 0.408163, and Queue2 loses no jobs.
     */
    @Test
    public void testNegativeSignalFiresOnceAtFirstStation() {
        NetworkAvgTable t = solve(build(1.0, 0.3, 2.0, 3.0, SignalType.NEGATIVE));
        double q1 = metric(t, t.getQLen(), "Queue1");
        double q2 = metric(t, t.getQLen(), "Queue2");
        double t1 = metric(t, t.getTput(), "Queue1");
        double t2 = metric(t, t.getTput(), "Queue2");
        assertEquals(0.769231, q1, 0.03, "Queue1 queue length");
        assertEquals(0.408163, q2, 0.03, "Queue2 queue length");
        assertEquals(t1, t2, 0.01,
                "the signal is annihilated at Queue1, so Queue2 must not lose jobs");
    }

    /**
     * Catastrophe, lambdaPos = 1, lambdaNeg = 0.3, mu1 = 2, mu2 = 3.
     * Exact: Q1 = 0.666667, Q2 = 0.363636, Queue1 throughput 0.8.
     */
    @Test
    public void testCatastropheEmptiesStationIncludingInServiceJob() {
        NetworkAvgTable t = solve(build(1.0, 0.3, 2.0, 3.0, SignalType.CATASTROPHE));
        double q1 = metric(t, t.getQLen(), "Queue1");
        double q2 = metric(t, t.getQLen(), "Queue2");
        double t1 = metric(t, t.getTput(), "Queue1");
        double t2 = metric(t, t.getTput(), "Queue2");
        assertEquals(0.666667, q1, 0.03, "Queue1 queue length");
        assertEquals(0.363636, q2, 0.03, "Queue2 queue length");
        assertEquals(0.8, t1, 0.02, "Queue1 throughput");
        assertEquals(t1, t2, 0.01,
                "the catastrophe is annihilated at Queue1, so Queue2 must not lose jobs");
    }

    /**
     * Catastrophe at lambdaNeg = lambdaPos = 1, mu1 = 2, mu2 = 3.
     * Exact: Q1 = 0.414214, Q2 = 0.242641, Queue1 throughput 0.585786. This is
     * the case where the AG builder used to drive the isolated Queue1 process
     * to its truncation bound (49.5) by treating the signal as extra load.
     */
    @Test
    public void testCatastropheStrong() {
        NetworkAvgTable t = solve(build(1.0, 1.0, 2.0, 3.0, SignalType.CATASTROPHE));
        double q1 = metric(t, t.getQLen(), "Queue1");
        double q2 = metric(t, t.getQLen(), "Queue2");
        double t1 = metric(t, t.getTput(), "Queue1");
        double t2 = metric(t, t.getTput(), "Queue2");
        assertEquals(0.414214, q1, 0.03, "Queue1 queue length");
        assertEquals(0.242641, q2, 0.03, "Queue2 queue length");
        assertEquals(0.585786, t1, 0.02, "Queue1 throughput");
        assertEquals(t1, t2, 0.01,
                "the catastrophe is annihilated at Queue1, so Queue2 must not lose jobs");
    }
}
