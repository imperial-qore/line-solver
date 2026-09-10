package jline.solvers.mva;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.RoutingStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Router;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.solvers.NetworkAvgTable;
import jline.solvers.mam.SolverMAM;
import org.junit.jupiter.api.Test;

/**
 * Round-robin dispatching in the two decomposition methods, MVA 'qna' and MAM
 * 'mna'. A round-robin node splits the departure stream one-in-k, so each
 * destination sees the k-fold convolution of the interarrival time: with a
 * Poisson source and k=2 the arrivals at each queue are Erlang-2 (SCV 1/2),
 * and the queue is shorter than under Bernoulli routing at the same rate.
 * Reference values from the MATLAB twin; RAND must return the exact M/M/1
 * value, which pins the k=1 branch to the original Markovian formula.
 */
public class RoundRobinSplitTest {

    private static final double TOL = 1e-3;

    private static Network rrModel(RoutingStrategy strategy) {
        Network model = new Network("rr_split");
        Source source = new Source(model, "Source");
        Router router = new Router(model, "Router");
        Queue q1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue q2 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass oclass = new OpenClass(model, "Class1");
        source.setArrival(oclass, new Exp(1.4));
        q1.setService(oclass, new Exp(1.0));
        q2.setService(oclass, new Exp(1.0));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(oclass, oclass, source, router, 1.0);
        P.set(oclass, oclass, router, q1, 1.0);
        P.set(oclass, oclass, router, q2, 1.0);
        P.set(oclass, oclass, q1, sink, 1.0);
        P.set(oclass, oclass, q2, sink, 1.0);
        model.link(P);
        router.setRouting(oclass, strategy);
        return model;
    }

    @Test
    public void qnaSeparatesRoundRobinFromBernoulli() {
        SolverMVA rr = new SolverMVA(rrModel(RoutingStrategy.RROBIN));
        rr.options.method = "qna";
        NetworkAvgTable trr = rr.getAvgTable();
        assertEquals(1.9250, trr.getQLen().get(1), TOL);
        assertEquals(1.9250, trr.getQLen().get(2), TOL);
        assertEquals(0.7, trr.getTput().get(1), TOL);
        assertEquals(0.7, trr.getTput().get(2), TOL);

        SolverMVA rand = new SolverMVA(rrModel(RoutingStrategy.RAND));
        rand.options.method = "qna";
        NetworkAvgTable trand = rand.getAvgTable();
        // Bernoulli split of a Poisson stream is Poisson: exact M/M/1
        assertEquals(2.3333, trand.getQLen().get(1), TOL);
        assertEquals(2.3333, trand.getQLen().get(2), TOL);
        assertTrue(trr.getQLen().get(1) < trand.getQLen().get(1));
    }

    @Test
    public void mnaSeparatesRoundRobinFromBernoulli() {
        SolverMAM rr = new SolverMAM(rrModel(RoutingStrategy.RROBIN));
        rr.options.method = "mna";
        NetworkAvgTable trr = rr.getAvgTable();
        // E_2/M/1 solved by MMAPPH1FCFS on the two-moment fit of the split flow
        assertEquals(1.8204, trr.getQLen().get(1), TOL);
        assertEquals(1.8204, trr.getQLen().get(2), TOL);
        assertEquals(0.7, trr.getTput().get(1), TOL);
        assertEquals(0.7, trr.getTput().get(2), TOL);

        SolverMAM rand = new SolverMAM(rrModel(RoutingStrategy.RAND));
        rand.options.method = "mna";
        NetworkAvgTable trand = rand.getAvgTable();
        assertEquals(2.3333, trand.getQLen().get(1), TOL);
        assertEquals(2.3333, trand.getQLen().get(2), TOL);
        assertTrue(trr.getQLen().get(1) < trand.getQLen().get(1));
    }
}
