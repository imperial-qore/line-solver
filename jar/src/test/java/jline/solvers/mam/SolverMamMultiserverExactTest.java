package jline.solvers.mam;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.Test;

import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.constant.SchedStrategy;
import jline.lang.processes.Erlang;
import jline.lang.processes.Exp;
import jline.lang.processes.HyperExp;
import jline.lang.processes.MMPP2;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.util.matrix.Matrix;

/**
 * SolverMAM must answer a multiserver FCFS station exactly when it can.
 *
 * Two fast paths are asserted here, both of which the generic path would
 * otherwise answer with the single-fast-server surrogate -- service divided by
 * the server count plus a surrogate delay, which discards arrival correlation
 * and the service shape alike:
 *
 *  - MAP/M/c, when the arrival stream is correlated and so refused by the
 *    renewal-gated PH/M/c path. On the MMPP2 model below the surrogate read
 *    1.757271 against the exact 1.521875, 15.5% high.
 *  - MAP/PH/c, whenever the service law is phase type rather than exponential.
 *
 * THE ORACLE IS THE MATLAB REFERENCE, which reaches these numbers through its
 * own Q-MAM and multiset-QBD routes.
 */
public class SolverMamMultiserverExactTest {

    private static Network mapMcModel(int c) {
        Network model = new Network("mapmc");
        Source src = new Source(model, "Src");
        Queue q = new Queue(model, "Q1", SchedStrategy.FCFS);
        q.setNumberOfServers(c);
        Sink snk = new Sink(model, "Snk");
        OpenClass cls = new OpenClass(model, "C1");
        src.setArrival(cls, new MMPP2(1.8, 0.4, 0.15, 0.25));
        q.setService(cls, new Exp(1.0));
        model.link(model.serialRouting(src, q, snk));
        return model;
    }

    @Test
    public void correlatedArrivalAtThreeServersIsExact() {
        SolverMAM solver = new SolverMAM(mapMcModel(3));
        Matrix Q = solver.getAvgQLen();
        assertEquals(1.521875, Q.get(1, 0), 1e-5);
    }

    @Test
    public void phaseTypeServiceAtTwoServersIsExact() {
        Network model = new Network("mphc");
        Source src = new Source(model, "S");
        Queue q = new Queue(model, "Q", SchedStrategy.FCFS);
        q.setNumberOfServers(2);
        Sink snk = new Sink(model, "K");
        OpenClass cls = new OpenClass(model, "C");
        src.setArrival(cls, new Exp(0.9));
        q.setService(cls, Erlang.fitMeanAndSCV(1.0, 0.5));
        model.link(model.serialRouting(src, q, snk));
        Matrix Q = new SolverMAM(model).getAvgQLen();
        assertEquals(1.075875, Q.get(1, 0), 1e-5);
    }

    @Test
    public void correlatedArrivalAndPhaseTypeServiceAtThreeServersIsExact() {
        Network model = new Network("mapphc");
        Source src = new Source(model, "S");
        Queue q = new Queue(model, "Q", SchedStrategy.FCFS);
        q.setNumberOfServers(3);
        Sink snk = new Sink(model, "K");
        OpenClass cls = new OpenClass(model, "C");
        src.setArrival(cls, new MMPP2(1.8, 0.4, 0.15, 0.25));
        q.setService(cls, HyperExp.fitMeanAndSCV(0.8, 4.0));
        model.link(model.serialRouting(src, q, snk));
        Matrix Q = new SolverMAM(model).getAvgQLen();
        assertEquals(1.145311, Q.get(1, 0), 1e-5);
    }
}
