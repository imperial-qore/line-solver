package jline.opt;

import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.opt.objectives.MinimizeCost;
import jline.opt.objectives.UtilizationConstraint;
import jline.opt.results.OptimizationResult;
import jline.opt.solver.BisectionSolver;
import jline.opt.variables.ServerAllocation;
import org.junit.jupiter.api.Test;

import java.util.LinkedHashMap;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Exact server sizing via bisection: matches the native Python
 * {@code opt_bisection_sizing.py} reference (6 servers, cost 60, 4 evaluations).
 */
public class BisectionE2ETest {

    @Test
    public void testBisectionSizing() {
        Network model = new Network("MMc");
        Source source = new Source(model, "Arrivals");
        Queue queue = new Queue(model, "Server", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Departures");
        OpenClass jobs = new OpenClass(model, "Jobs");
        source.setArrival(jobs, new Exp(3.0));
        queue.setService(jobs, new Exp(1.0));
        model.link(Network.serialRouting(source, queue, sink));

        OptimizationProblem problem = new OptimizationProblem(model);
        Map<jline.lang.nodes.Station, Double> serverCost =
                new LinkedHashMap<jline.lang.nodes.Station, Double>();
        serverCost.put(queue, 10.0);
        problem.addVariable(new ServerAllocation(queue, 1, 10));
        problem.setObjective(new MinimizeCost(serverCost, null, null, null));
        problem.addConstraint(new UtilizationConstraint(queue, 0.5));

        OptimizationResult result = new BisectionSolver(problem).solve();

        assertEquals(6, ((Number) result.getVariableValue("Server_servers")).intValue());
        assertEquals(60.0, result.objectiveValue, 1e-9);
        assertTrue(result.feasible);
        assertEquals(4, result.modelEvaluations, "O(log n) probe count");
    }
}
