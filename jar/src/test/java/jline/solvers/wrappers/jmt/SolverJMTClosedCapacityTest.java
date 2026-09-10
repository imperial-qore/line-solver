package jline.solvers.wrappers.jmt;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.DropStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertDoesNotThrow;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * BUG-81's JMT half: a binding station capacity a CLOSED class can reach must be
 * refused, not silently exported as something JMT does not enforce.
 * <p>
 * LINE blocks a closed job that finds no room -- the upstream departure is
 * disabled and the job stays where it is -- and no JMT drop strategy reproduces
 * that. JMT's queue section accepts exactly five ({@code drop}, {@code block},
 * {@code retrial}, {@code BAS blocking}, {@code waiting queue}), and neither of
 * the two a WAITQ station maps onto matches, measured on the fixture below:
 * </p><ul>
 * <li>{@code waiting queue}, which is what was written, does not enforce
 * {@code size} at all -- JSIM returned the UNCONSTRAINED [2.03 1.99 1.98],
 * X = 0.750, against the exact [3.6090 0.9711 1.4199], X = 0.6522.</li>
 * <li>{@code BAS blocking} enforces it but completes the service BEFORE
 * blocking, so the blocked job moves the instant room frees: [2.871 1.373 1.756],
 * X = 0.7126. A different queueing model, not a rounding.</li>
 * </ul>
 */
public class SolverJMTClosedCapacityTest {

    /** Closed 3-queue tandem, Exp(1) FCFS everywhere, Q2 capped at 2. */
    private static Network cappedTandem() {
        Network model = new Network("tandem");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.FCFS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.FCFS);
        Queue q3 = new Queue(model, "Q3", SchedStrategy.FCFS);
        ClosedClass c1 = new ClosedClass(model, "C1", 6, q1, 0);
        q1.setService(c1, new Exp(1.0));
        q2.setService(c1, new Exp(1.0));
        q3.setService(c1, new Exp(1.0));
        q2.setClassCapacity(c1, 2);
        model.link(model.serialRouting(q1, q2, q3));
        return model;
    }

    /** M/M/1/K: an OPEN refused arrival is LOST, which JMT does reproduce. */
    private static Network mm1k(int K) {
        Network model = new Network("mm1k");
        Source source = new Source(model, "Source");
        Queue q = new Queue(model, "Q", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass c1 = new OpenClass(model, "C1");
        source.setArrival(c1, new Exp(0.8));
        q.setService(c1, new Exp(1.0));
        q.setClassCapacity(c1, K);
        q.setDropRule(c1, DropStrategy.Drop);
        model.link(model.serialRouting(source, q, sink));
        return model;
    }

    @Test
    public void testBindingClosedCapacityIsRefused() {
        SolverJMT solver = new SolverJMT(cappedTandem());
        solver.options.samples = 10000;
        solver.options.seed = 23000;
        RuntimeException e = assertThrows(RuntimeException.class, solver::getAvgQLen);
        String msg = String.valueOf(e.getMessage());
        assertTrue(msg.contains("closed class") || msg.contains("Q2"),
                "refusal should name the station and the closed class, got: " + msg);
    }

    /**
     * The regression guard: the open half must still export and solve. If the
     * refusal widened to open classes, every M/M/1/K model would stop working.
     */
    @Test
    public void testOpenClassCapacityStillExports() {
        SolverJMT solver = new SolverJMT(mm1k(2));
        solver.options.samples = 20000;
        solver.options.seed = 23000;
        assertDoesNotThrow(solver::getAvgQLen);
    }

    /**
     * A Delay+Queue cycle whose two classes share ONE chain, so both are served
     * at the Queue and neither declares a capacity.
     */
    private static Network twoClassOneChain(int nEach) {
        Network model = new Network("twoclass");
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue1", SchedStrategy.PS);
        ClosedClass c1 = new ClosedClass(model, "Class1", nEach, delay, 0);
        ClosedClass c2 = new ClosedClass(model, "Class2", nEach, delay, 0);
        delay.setService(c1, new Exp(1.0));
        delay.setService(c2, new Exp(1.0));
        queue.setService(c1, new Exp(2.0));
        queue.setService(c2, new Exp(2.0));
        // Every (class, node) row sums to 1, and each class switches into the
        // other with probability 1/2 at both stations, which is what puts the
        // two of them in ONE chain.
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(c1, c1, delay, queue, 0.5);
        P.set(c1, c2, delay, queue, 0.5);
        P.set(c2, c1, delay, queue, 0.5);
        P.set(c2, c2, delay, queue, 0.5);
        P.set(c1, c1, queue, delay, 0.5);
        P.set(c1, c2, queue, delay, 0.5);
        P.set(c2, c1, queue, delay, 0.5);
        P.set(c2, c2, queue, delay, 0.5);
        model.link(P);
        return model;
    }

    /**
     * refreshCapacity DERIVES sn.cap for a station nobody capped, as
     * min(sum_c chaincap, sum_r classcap) -- and the classcap row is the CHAIN
     * population repeated per class, so a Queue serving two classes of one
     * 4-job chain carries 8, which those 4 jobs can never reach.
     * <p>
     * The reachability test in saveBufferCapacity was `cap == sum(njobs)`, an
     * EQUALITY only the single-class case satisfies, so every multi-class
     * station fell through to jmtStationCapAssert and was refused as a binding
     * buffer the model never declared. It is `>=` now.
     */
    @Test
    public void testDerivedMulticlassCapacityIsNotABuffer() {
        Network model = twoClassOneChain(2);
        NetworkStruct sn = model.getStruct(true);
        assertEquals(1, sn.nchains, "the two classes must share one chain");
        int ist = (int) sn.nodeToStation.get(model.getNodeIndex("Queue1"));
        assertEquals(8.0, sn.cap.get(ist), 0.0,
                "(2 classes served) x (chain population 4), all derived");

        SolverJMT solver = new SolverJMT(model);
        solver.options.samples = 10000;
        solver.options.seed = 23000;
        assertDoesNotThrow(solver::getAvgQLen);
    }

    /**
     * The same equality's other false positive: 4 jobs cannot fill 100 either,
     * so a DECLARED capacity above the population is not a buffer.
     */
    @Test
    public void testDeclaredCapacityAboveThePopulationIsNotABuffer() {
        Network model = twoClassOneChain(2);
        ((Queue) model.getNodeByName("Queue1")).setCapacity(100);

        SolverJMT solver = new SolverJMT(model);
        solver.options.samples = 10000;
        solver.options.seed = 23000;
        assertDoesNotThrow(solver::getAvgQLen);
    }
}
