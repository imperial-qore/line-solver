package jline.solvers.ln;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.constant.SchedStrategy;
import jline.lang.layered.Activity;
import jline.lang.layered.Entry;
import jline.lang.layered.LayeredNetwork;
import jline.lang.layered.Processor;
import jline.lang.layered.Task;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.processes.Exp;
import jline.solvers.LayeredNetworkAvgTable;
import jline.solvers.wrappers.jmt.SolverJMT;
import jline.solvers.ldes.SolverLDES;
import jline.solvers.mva.SolverMVA;
import jline.solvers.nc.SolverNC;
import jline.solvers.ssa.SolverSSA;

import org.junit.jupiter.api.Test;

import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Tests for the stochastic iteration support: Solver.isStochastic
 * classification and the SolverLN Robbins-Monro controller (convergedStoch).
 */
public class StochIterTest extends SolverLNTestBase {

    private Network buildQN() {
        Network model = new Network("qn");
        Delay node1 = new Delay(model, "D");
        Queue node2 = new Queue(model, "Q", SchedStrategy.PS);
        ClosedClass cclass = new ClosedClass(model, "C", 3, node1, 0);
        node1.setService(cclass, new Exp(1));
        node2.setService(cclass, new Exp(2));
        model.link(Network.serialRouting(node1, node2));
        return model;
    }

    private LayeredNetwork buildLQN() throws Exception {
        LayeredNetwork model = new LayeredNetwork("lqntest");

        Processor P1 = new Processor(model, "P1", 1, SchedStrategy.PS);
        Processor P2 = new Processor(model, "P2", 1, SchedStrategy.PS);

        Task T1 = new Task(model, "T1", 5, SchedStrategy.REF);
        T1.on(P1);
        T1.setThinkTime(Exp.fitMean(10));

        Task T2 = new Task(model, "T2", 1, SchedStrategy.FCFS);
        T2.on(P2);

        Entry E1 = new Entry(model, "E1");
        E1.on(T1);
        Entry E2 = new Entry(model, "E2");
        E2.on(T2);

        Activity A1 = new Activity(model, "A1", Exp.fitMean(1));
        A1.on(T1);
        A1.boundTo(E1);
        A1.synchCall(E2, 1);

        Activity A2 = new Activity(model, "A2", Exp.fitMean(0.5));
        A2.on(T2);
        A2.repliesTo(E2);
        A2.boundTo(E2);

        return model;
    }

    @Test
    public void testIsStochasticClassification() {
        Network model = buildQN();
        assertFalse(new SolverMVA(model).isStochastic(), "MVA should be deterministic");
        assertFalse(new SolverNC(model).isStochastic(), "NC default should be deterministic");
        assertTrue(new SolverNC(model, "imci").isStochastic(), "NC imci should be stochastic");
        assertTrue(new SolverNC(model, "ls").isStochastic(), "NC ls should be stochastic");
        assertTrue(new SolverSSA(model).isStochastic(), "SSA should be stochastic");
        assertTrue(new SolverLDES(model).isStochastic(), "LDES should be stochastic");
        SolverJMT jmt = new SolverJMT(model);
        assertTrue(jmt.isStochastic(), "JMT default should be stochastic");
        assertFalse(jmt.isStochasticMethod("jmva.mva"), "JMT jmva.mva should be deterministic");
        assertTrue(jmt.isStochasticMethod("jmva.ls"), "JMT jmva.ls should be stochastic");
        assertTrue(jmt.isStochasticMethod("default/jsim"), "runtime-resolved jsim should be stochastic");
        assertFalse(new SolverNC(model).isStochasticMethod("default/comom"),
                "runtime-resolved comom should be deterministic");
        assertTrue(new SolverNC(model).isStochasticMethod("default/imci"),
                "runtime-resolved imci should be stochastic");
    }

    @Test
    public void testDeterministicLayersResolveOff() throws Exception {
        SolverLN solver = new SolverLN(buildLQN());
        suppressOutput(() -> solver.getEnsembleAvg());
        assertEquals("off", solver.stochiterMode, "MVA layers must resolve stochiter mode off");
    }

    @Test
    public void testForcedRobbinsMonroMatchesDeterministicFixedPoint() throws Exception {
        // Reference: deterministic Picard iteration
        SolverLN reference = new SolverLN(buildLQN());
        LayeredNetworkAvgTable refTable =
                suppressOutput(() -> (LayeredNetworkAvgTable) reference.getEnsembleAvg());
        double tputRef = tputOf(refTable, "T1");

        // Forced Robbins-Monro mode on the same (deterministic) layers: the
        // controller must start Polyak averaging and converge to the same
        // fixed point, since the noise is zero.
        SolverLN solver = new SolverLN(buildLQN());
        solver.options.config.stochiter = "rm";
        LayeredNetworkAvgTable rmTable =
                suppressOutput(() -> (LayeredNetworkAvgTable) solver.getEnsembleAvg());
        assertEquals("rm", solver.stochiterMode);
        assertTrue(solver.stochAvgCount > 0, "Polyak averaging never started");
        assertTrue(solver.hasconverged, "RM iteration did not converge");
        double tputRM = tputOf(rmTable, "T1");
        assertEquals(tputRef, tputRM, 0.02 * tputRef,
                "RM solution deviates from the deterministic fixed point");
    }

    private static double tputOf(LayeredNetworkAvgTable table, String nodeName) {
        List<String> names = table.getNodeNames();
        List<Double> tputs = table.getTput();
        for (int i = 0; i < names.size(); i++) {
            if (names.get(i).equals(nodeName)) {
                return tputs.get(i);
            }
        }
        throw new IllegalArgumentException("Node not found: " + nodeName);
    }
}
