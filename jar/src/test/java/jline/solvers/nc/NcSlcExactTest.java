package jline.solvers.nc;

import jline.lang.*;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.*;
import jline.lang.processes.Exp;
import jline.solvers.NetworkAvgTable;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

/**
 * NC exact on a single-server closed model with self-looping classes must use
 * the standard normalizing-constant path, not the load-dependent one (which
 * crashed in Pfqn_nc_sanitize). CTMC-exact Class1 QLen = [0.714286, 0.142857,
 * 0.142857]; the two SLC classes each hold 1 job at their queue.
 */
public class NcSlcExactTest {
    @Test
    public void ncExactSelfLoopingClosed() {
        Network m = new Network("slc");
        Delay d = new Delay(m, "D");
        Queue q1 = new Queue(m, "Q1", SchedStrategy.PS);
        Queue q2 = new Queue(m, "Q2", SchedStrategy.FCFS);
        ClosedClass c = new ClosedClass(m, "C", 1, d);
        SelfLoopingClass s1 = new SelfLoopingClass(m, "S1", 1, q1);
        SelfLoopingClass s2 = new SelfLoopingClass(m, "S2", 1, q2);
        d.setService(c, new Exp(1.0)); q1.setService(c, new Exp(1.5)); q2.setService(c, new Exp(1.5));
        q1.setService(s1, new Exp(1.5)); q2.setService(s2, new Exp(1.5));
        RoutingMatrix P = m.initRoutingMatrix();
        P.set(c, new Matrix(new double[][]{{0.7,0.15,0.15},{1,0,0},{1,0,0}}));
        P.set(s1, new Matrix(new double[][]{{0,0,0},{0,1,0},{0,0,0}}));
        P.set(s2, new Matrix(new double[][]{{0,0,0},{0,0,0},{0,0,1}}));
        m.link(P);
        NetworkAvgTable t = new SolverNC(m, new NCOptions().method("exact")).getAvgTable();
        // Rows are station-class; Class1 rows are Delay, Q1, Q2 (indices 0,1,3).
        assertEquals(0.714286, t.getQLen().get(0), 1e-4, "Class1 QLen at Delay");
        assertEquals(0.142857, t.getQLen().get(1), 1e-4, "Class1 QLen at Q1");
        assertEquals(0.142857, t.getQLen().get(3), 1e-4, "Class1 QLen at Q2");
    }
}
