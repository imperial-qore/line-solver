package jline.solvers.mam;
import jline.lang.*;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.*;
import jline.lang.processes.Exp;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

public class MnaSlcGateTest {
    @Test
    public void mnaRejectsSelfLoopingClass() {
        Network m = new Network("slc");
        Delay d = new Delay(m, "Delay");
        Queue q1 = new Queue(m, "Q1", SchedStrategy.PS);
        Queue q2 = new Queue(m, "Q2", SchedStrategy.FCFS);
        ClosedClass c = new ClosedClass(m, "C", 1, d);
        SelfLoopingClass s1 = new SelfLoopingClass(m, "SLC1", 1, q1);
        SelfLoopingClass s2 = new SelfLoopingClass(m, "SLC2", 1, q2);
        d.setService(c, new Exp(1.0)); q1.setService(c, new Exp(1.5)); q2.setService(c, new Exp(1.5));
        q1.setService(s1, new Exp(1.5)); q2.setService(s2, new Exp(1.5));
        RoutingMatrix P = m.initRoutingMatrix();
        P.set(c, new Matrix(new double[][]{{0.7,0.15,0.15},{1,0,0},{1,0,0}}));
        P.set(s1, new Matrix(new double[][]{{0,0,0},{0,1,0},{0,0,0}}));
        P.set(s2, new Matrix(new double[][]{{0,0,0},{0,0,0},{0,0,1}}));
        m.link(P);
        String reason = new SolverMAM(m, "mna").supportsModelMethod("mna");
        assertTrue(reason.contains("self-looping"), "mna must reject SLC, got: [" + reason + "]");
    }

    @Test
    public void mnaAcceptsPlainMulticlass() {
        // Both classes reference the Delay and traverse Delay->Queue: there is
        // inter-station flow, so mna must NOT be rejected (guards against the
        // gate over-firing on ordinary closed classes).
        Network m = new Network("plain");
        Delay d = new Delay(m, "D");
        Queue q = new Queue(m, "Q", SchedStrategy.FCFS);
        ClosedClass a = new ClosedClass(m, "A", 1, d);
        ClosedClass b = new ClosedClass(m, "B", 1, d);
        d.setService(a, new Exp(1.0)); q.setService(a, new Exp(1.5));
        d.setService(b, new Exp(1.0)); q.setService(b, new Exp(1.5));
        RoutingMatrix P = m.initRoutingMatrix();
        P.set(a, m.serialRouting(d, q));
        P.set(b, m.serialRouting(d, q));
        m.link(P);
        assertEquals("", new SolverMAM(m, "mna").supportsModelMethod("mna"),
                "mna must accept a plain multiclass closed model");
    }
}
