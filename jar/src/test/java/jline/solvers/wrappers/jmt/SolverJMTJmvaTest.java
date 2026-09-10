package jline.solvers.wrappers.jmt;

import jline.VerboseLevel;
import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.processes.Exp;
import jline.solvers.NetworkAvgTable;
import jline.solvers.SolverOptions;
import jline.solvers.mva.SolverMVA;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import java.util.List;

import static jline.TestTools.LOOSE_FINE_TOL;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Coverage for the JMVA analytical engine behind SolverJMT.
 * <p>
 * The eleven {@code jmva*} method aliases had NO test in this suite, which is
 * how a fully broken JMVA path survived: JMVA reports per CHAIN
 * ({@code customerclass="Chain01"}), the reader stored that name verbatim, and
 * the consumer matched it against model class names, so every cell stayed at
 * its initialised zero under a "completed" status.
 * <p>
 * The second test is the load-bearing one. A model with one class per chain
 * makes the chain-to-class disaggregation an identity, so it would pass over a
 * broken implementation that merely renamed the chain; the class-switching
 * model puts two classes with DIFFERENT service rates in a single chain, where
 * the rescaling by ST/STchain, Vchain and alpha actually has to be right.
 * <p>
 * JMVA is exact on a product-form model, so exact MVA is the reference.
 */
public class SolverJMTJmvaTest {

    private static NetworkAvgTable solveJmva(Network model) {
        SolverOptions options = new SolverOptions(SolverType.JMT);
        options.method = "jmva";
        options.seed = 23000;
        options.verbose = VerboseLevel.SILENT;
        return new SolverJMT(model, options).getAvgTable();
    }

    private static NetworkAvgTable solveExactMva(Network model) {
        SolverOptions options = new SolverOptions(SolverType.MVA);
        options.verbose = VerboseLevel.SILENT;
        return new SolverMVA(model, options).getAvgTable();
    }

    /** Closed product-form Delay + FCFS Queue, single class. */
    private static Network singleClassModel() {
        Network model = new Network("jmva_single_class");
        Delay think = new Delay(model, "Think");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        ClosedClass jobClass = new ClosedClass(model, "Class1", 6, think, 0);
        think.setService(jobClass, Exp.fitRate(1.0));
        queue.setService(jobClass, Exp.fitRate(1.5));
        model.link(Network.serialRouting(think, queue));
        return model;
    }

    /**
     * Two classes in ONE chain via class switching, with different service
     * rates per class, so chains and classes genuinely differ.
     */
    private static Network switchingModel() {
        Network model = new Network("jmva_class_switch");
        Delay think = new Delay(model, "Think");
        Queue queue = new Queue(model, "Queue", SchedStrategy.PS);
        ClosedClass classA = new ClosedClass(model, "ClassA", 4, think, 0);
        ClosedClass classB = new ClosedClass(model, "ClassB", 0, think, 0);
        think.setService(classA, Exp.fitRate(1.0));
        think.setService(classB, Exp.fitRate(2.0));
        queue.setService(classA, Exp.fitRate(3.0));
        queue.setService(classB, Exp.fitRate(1.5));
        RoutingMatrix routing = model.initRoutingMatrix();
        routing.addConnection(classA, classB, think, queue, 1.0);
        routing.addConnection(classB, classA, queue, think, 1.0);
        model.link(routing);
        return model;
    }

    /**
     * Two classes in ONE chain visiting a MULTISERVER station, which is the only
     * configuration that reaches the finite-server branch of the Utilization
     * rescaling, {@code min(sum(finite NK), c_i) / c_i}. That factor is 1 at
     * c_i = 1, so every single-server model passes over a wrong multiserver
     * branch.
     * <p>
     * The switch is placed on the RETURN leg so that both classes visit the
     * multiserver queue; switching on the way in would leave only one class
     * there and the disaggregation would go untested at that station. Service
     * at the FCFS multiserver station is class-INDEPENDENT, which product form
     * requires; the class-dependent rates sit at the infinite-server Think,
     * where they are admissible and still exercise alpha.
     */
    private static Network multiserverSwitchingModel() {
        Network model = new Network("jmva_multiserver_chain");
        Delay think = new Delay(model, "Think");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        queue.setNumberOfServers(2);
        ClosedClass classA = new ClosedClass(model, "ClassA", 5, think, 0);
        ClosedClass classB = new ClosedClass(model, "ClassB", 0, think, 0);
        think.setService(classA, Exp.fitRate(1.0));
        think.setService(classB, Exp.fitRate(2.0));
        queue.setService(classA, Exp.fitRate(3.0));
        queue.setService(classB, Exp.fitRate(3.0));
        RoutingMatrix routing = model.initRoutingMatrix();
        routing.addConnection(classA, classA, think, queue, 1.0);
        routing.addConnection(classB, classB, think, queue, 1.0);
        routing.addConnection(classA, classB, queue, think, 1.0);
        routing.addConnection(classB, classA, queue, think, 1.0);
        model.link(routing);
        return model;
    }

    @Test
    @DisplayName("jmva matches exact MVA at a multiserver station with two classes in one chain")
    public void testJmvaMultiserverUtilizationBranch() {
        NetworkAvgTable jmva = solveJmva(multiserverSwitchingModel());
        NetworkAvgTable exact = solveExactMva(multiserverSwitchingModel());

        assertEquals(exact.getQLen().size(), jmva.getQLen().size(),
                "jmva and exact MVA must report the same station-class rows");

        List<Double> jQLen = jmva.getQLen(), eQLen = exact.getQLen();
        List<Double> jUtil = jmva.getUtil(), eUtil = exact.getUtil();
        List<Double> jRespT = jmva.getRespT(), eRespT = exact.getRespT();
        List<Double> jTput = jmva.getTput(), eTput = exact.getTput();

        for (int i = 0; i < eQLen.size(); i++) {
            assertEquals(eQLen.get(i), jQLen.get(i), LOOSE_FINE_TOL, "QLen row " + i);
            assertEquals(eUtil.get(i), jUtil.get(i), LOOSE_FINE_TOL, "Util row " + i);
            assertEquals(eRespT.get(i), jRespT.get(i), LOOSE_FINE_TOL, "RespT row " + i);
            assertEquals(eTput.get(i), jTput.get(i), LOOSE_FINE_TOL, "Tput row " + i);
        }
    }

    @Test
    @DisplayName("jmva returns a non-degenerate solution, not an all-zero table")
    public void testJmvaSingleClassIsNotAllZero() {
        NetworkAvgTable jmva = solveJmva(singleClassModel());
        List<Double> qlen = jmva.getQLen();

        double total = 0.0;
        for (int i = 0; i < qlen.size(); i++) {
            total += qlen.get(i);
        }
        assertTrue(total > 0.0,
                "jmva returned an all-zero queue length table; the chain-to-class"
                        + " mapping in getResultsJMVA has regressed");
        // Closed model: the queue lengths must account for the whole population.
        assertEquals(6.0, total, LOOSE_FINE_TOL,
                "jmva queue lengths must sum to the closed population");
    }

    @Test
    @DisplayName("jmva matches exact MVA per metric with two classes in one chain")
    public void testJmvaChainToClassDisaggregation() {
        NetworkAvgTable jmva = solveJmva(switchingModel());
        NetworkAvgTable exact = solveExactMva(switchingModel());

        assertEquals(exact.getQLen().size(), jmva.getQLen().size(),
                "jmva and exact MVA must report the same station-class rows");

        List<Double> jQLen = jmva.getQLen(), eQLen = exact.getQLen();
        List<Double> jUtil = jmva.getUtil(), eUtil = exact.getUtil();
        List<Double> jRespT = jmva.getRespT(), eRespT = exact.getRespT();
        List<Double> jTput = jmva.getTput(), eTput = exact.getTput();

        for (int i = 0; i < eQLen.size(); i++) {
            assertEquals(eQLen.get(i), jQLen.get(i), LOOSE_FINE_TOL, "QLen row " + i);
            assertEquals(eUtil.get(i), jUtil.get(i), LOOSE_FINE_TOL, "Util row " + i);
            assertEquals(eRespT.get(i), jRespT.get(i), LOOSE_FINE_TOL, "RespT row " + i);
            assertEquals(eTput.get(i), jTput.get(i), LOOSE_FINE_TOL, "Tput row " + i);
        }
    }
}
