package jline.solvers.mam;

import jline.lang.*;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.*;
import jline.lang.processes.Exp;
import jline.solvers.NetworkAvgTable;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

/**
 * dec.mmap indexes DEP/V by node while solver_mam_traffic reads them by
 * station. The two coincide only when every station precedes every non-station
 * node in node order. Declaring the Sink before the Queue breaks that, so
 * nodeToStation returns a non-station marker and the sweep blows up.
 * Both node orders must give the same (exact Jackson) answer.
 */
public class DecMmapNodeOrderTest {

    private Network build(boolean sinkFirst) {
        Network model = new Network("t");
        Source source = new Source(model, "S");
        Sink sink;
        Queue queue;
        if (sinkFirst) {
            sink = new Sink(model, "K");
            queue = new Queue(model, "Q1", SchedStrategy.FCFS);
        } else {
            queue = new Queue(model, "Q1", SchedStrategy.FCFS);
            sink = new Sink(model, "K");
        }
        OpenClass oc = new OpenClass(model, "C", 0);
        source.setArrival(oc, new Exp(0.5));
        queue.setService(oc, new Exp(1.0));
        model.link(Network.serialRouting(source, queue, sink));
        return model;
    }

    @Test
    public void testDecMmapLuckyOrder() {
        SolverMAM s = new SolverMAM(build(false), "dec.mmap");
        assertEquals(1.0, s.getAvgTable().getQLen().get(1), 1e-6, "Source,Queue,Sink");
    }

    @Test
    public void testDecMmapSinkFirst() {
        SolverMAM s = new SolverMAM(build(true), "dec.mmap");
        assertEquals(1.0, s.getAvgTable().getQLen().get(1), 1e-6, "Source,Sink,Queue");
    }

    @Test
    public void testDecMmapIndependentOfNodeOrder() {
        double[] q = new double[2];
        for (int i = 0; i < 2; i++) {
            SolverMAM s = new SolverMAM(build(i == 1), "dec.mmap");
            NetworkAvgTable t = s.getAvgTable();
            assertNotNull(t, "dec.mmap must solve regardless of node order");
            q[i] = t.getQLen().get(1);
        }
        // M/M/1 with rho=0.5 is exactly rho/(1-rho) = 1.0 (Burke)
        assertEquals(1.0, q[0], 1e-6, "dec.mmap Source,Queue,Sink");
        assertEquals(1.0, q[1], 1e-6, "dec.mmap Source,Sink,Queue (sink declared first)");
        assertEquals(q[0], q[1], 1e-12, "dec.mmap must not depend on node declaration order");
    }
}
