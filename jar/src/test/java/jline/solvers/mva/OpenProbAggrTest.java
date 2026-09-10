package jline.solvers.mva;

import static org.junit.jupiter.api.Assertions.assertEquals;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.solvers.fluid.SolverFluid;
import org.junit.jupiter.api.Test;

/**
 * The mixed branch of getProbAggr/getProbSysAggr, which used to throw
 * "not yet implemented for models with open classes" in the JAR while MATLAB,
 * native Python and the C++ all served it.
 * <p>
 * A Source carries no population of its own (EXT), so its aggregate probability
 * is 1. An M/M/1 queue in the empty state has P = 1 - rho exactly, which is what
 * the multinomial-geometric product form reduces to for one open class; the
 * system probability is the same number, the Source contributing nothing.
 * Reference values agree with MATLAB and with native Python.
 */
public class OpenProbAggrTest {

    private static final double TOL = 1e-9;

    private static Network openModel() {
        Network model = new Network("mm1");
        Source source = new Source(model, "S");
        Queue queue = new Queue(model, "Q1", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass oclass = new OpenClass(model, "C1");
        source.setArrival(oclass, new Exp(0.5));
        queue.setService(oclass, new Exp(1.0));
        jline.lang.RoutingMatrix P = model.initRoutingMatrix();
        P.set(oclass, oclass, Network.serialRouting(source, queue, sink));
        model.link(P);
        return model;
    }

    private static Network closedModel() {
        Network model = new Network("cqn");
        Delay delay = new Delay(model, "Think");
        Queue queue = new Queue(model, "Q1", SchedStrategy.PS);
        ClosedClass cclass = new ClosedClass(model, "C1", 2, delay);
        delay.setService(cclass, new Exp(1.0));
        queue.setService(cclass, new Exp(2.0));
        jline.lang.RoutingMatrix P = model.initRoutingMatrix();
        P.set(cclass, cclass, Network.serialRouting(delay, queue));
        model.link(P);
        return model;
    }

    @Test
    public void mvaServesTheMixedBranch() {
        assertEquals(1.0, new SolverMVA(openModel()).getProbAggr(0).probability.get(0, 0), TOL);
        assertEquals(0.5, new SolverMVA(openModel()).getProbAggr(1).probability.get(0, 0), TOL);
        assertEquals(0.5, new SolverMVA(openModel()).getProbSysAggr().probability.get(0, 0), TOL);
    }

    @Test
    public void fluidServesTheMixedBranch() {
        // the fluid reaches rho by integration, so its 1 - rho carries the ODE
        // solver's residual (measured 1.8e-9) rather than the exact 0.5
        assertEquals(1.0, new SolverFluid(openModel()).getProbAggr(0).probability.get(0, 0), 1e-6);
        assertEquals(0.5, new SolverFluid(openModel()).getProbAggr(1).probability.get(0, 0), 1e-6);
    }

    /**
     * Both stations on a FRESH solver. Station 1 used to answer 0.16, the value
     * for station 0's state: the state row was fetched by walking
     * {@code sn.state.entrySet()} and counting to the stateful index, which
     * assumes an iteration order a HashMap does not provide.
     */
    @Test
    public void theClosedBranchIsUnchanged() {
        assertEquals(0.36, new SolverMVA(closedModel()).getProbAggr(0).probability.get(0, 0), 1e-6);
        assertEquals(0.36, new SolverMVA(closedModel()).getProbAggr(1).probability.get(0, 0), 1e-6);
        assertEquals(0.36, new SolverMVA(closedModel()).getProbSysAggr().probability.get(0, 0), 1e-6);
    }
}
