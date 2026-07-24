package jline.solvers.mva;

import static org.junit.jupiter.api.Assertions.assertEquals;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Erlang;
import jline.lang.processes.Exp;
import jline.lang.processes.HyperExp;
import jline.solvers.NetworkAvgTable;
import org.junit.jupiter.api.Test;

/**
 * Tests for SolverMVA with method='sum' (summation method SUM/ESUM plus
 * closing method for open and mixed models). Reference values pinned to
 * the MATLAB solver_mva_sum implementation.
 */
public class SolverMVASumTest {

    private static final double TOL = 1e-3;

    @Test
    public void testClosedNonProductForm() {
        // Delay + FCFS(Erlang scv=0.5) + FCFS m=2 (HyperExp scv=3), K=10
        Network model = new Network("sum_closed");
        Delay think = new Delay(model, "Think");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.FCFS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.FCFS);
        q2.setNumberOfServers(2);
        ClosedClass c = new ClosedClass(model, "C", 10, think, 0);
        think.setService(c, new Exp(1));
        q1.setService(c, Erlang.fitMeanAndSCV(0.3, 0.5));
        q2.setService(c, HyperExp.fitMeanAndSCV(0.8, 3));
        model.link(Network.serialRouting(think, q1, q2));

        SolverMVA solver = new SolverMVA(model);
        solver.options.method = "sum";
        NetworkAvgTable t = solver.getAvgTable();
        // MATLAB: X=2.2234, Q=[2.2234 1.5257 6.2509], U(Q2)=0.88934
        assertEquals(2.2234, t.getTput().get(0), TOL);
        assertEquals(2.2234, t.getQLen().get(0), TOL);
        assertEquals(1.5257, t.getQLen().get(1), TOL);
        assertEquals(6.2509, t.getQLen().get(2), TOL);
        assertEquals(0.88934, t.getUtil().get(2), TOL);
    }

    @Test
    public void testOpenNonProductForm() {
        // Source(1) -> FCFS(Erlang scv=0.5) -> FCFS(HyperExp scv=4) -> Sink
        Network model = new Network("sum_open");
        Source src = new Source(model, "Src");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.FCFS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.FCFS);
        Sink snk = new Sink(model, "Snk");
        OpenClass o = new OpenClass(model, "O");
        src.setArrival(o, new Exp(1));
        q1.setService(o, Erlang.fitMeanAndSCV(0.5, 0.5));
        q2.setService(o, HyperExp.fitMeanAndSCV(0.6, 4));
        model.link(Network.serialRouting(src, q1, q2, snk));

        SolverMVA solver = new SolverMVA(model);
        solver.options.method = "sum";
        NetworkAvgTable t = solver.getAvgTable();
        // MATLAB: Q=[0.87494 2.8483], U=[0.5 0.6], X=1
        assertEquals(0.87494, t.getQLen().get(1), TOL);
        assertEquals(2.8483, t.getQLen().get(2), TOL);
        assertEquals(0.5, t.getUtil().get(1), TOL);
        assertEquals(0.6, t.getUtil().get(2), TOL);
    }

    @Test
    public void testMixedNonProductForm() {
        // closed class (Delay+Q1) + open class (Src->Q1->Sink)
        Network model = new Network("sum_mixed");
        Source src = new Source(model, "Src");
        Delay th = new Delay(model, "Th");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.FCFS);
        Sink snk = new Sink(model, "Snk");
        ClosedClass c = new ClosedClass(model, "C", 5, th, 0);
        OpenClass o = new OpenClass(model, "O");
        src.setArrival(o, new Exp(0.4));
        th.setService(c, new Exp(1));
        q1.setService(c, Erlang.fitMeanAndSCV(0.4, 0.5));
        q1.setService(o, Erlang.fitMeanAndSCV(0.5, 0.5));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(c, Network.serialRouting(th, q1));
        P.set(o, Network.serialRouting(src, q1, snk));
        model.link(P);

        SolverMVA solver = new SolverMVA(model);
        solver.options.method = "sum";
        NetworkAvgTable t = solver.getAvgTable();
        // MATLAB: C: Q(Th)=1.622, Q(Q1)=3.378, U(Q1)=0.64879, X=1.622
        //         O: Q(Q1)=1.0413, U(Q1)=0.2, X=0.4
        int rTh = -1, rQ1c = -1, rQ1o = -1;
        for (int i = 0; i < t.getStationNames().size(); i++) {
            String st = t.getStationNames().get(i);
            String cl = t.getClassNames().get(i);
            if ("Th".equals(st) && "C".equals(cl)) rTh = i;
            if ("Q1".equals(st) && "C".equals(cl)) rQ1c = i;
            if ("Q1".equals(st) && "O".equals(cl)) rQ1o = i;
        }
        assertEquals(1.622, t.getQLen().get(rTh), TOL);
        assertEquals(3.378, t.getQLen().get(rQ1c), TOL);
        assertEquals(0.64879, t.getUtil().get(rQ1c), TOL);
        assertEquals(1.0413, t.getQLen().get(rQ1o), TOL);
        assertEquals(0.2, t.getUtil().get(rQ1o), TOL);
    }
}
