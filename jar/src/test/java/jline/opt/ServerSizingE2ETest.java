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
import jline.opt.solver.LineOptSolverOptions;
import jline.opt.variables.ServerAllocation;
import org.junit.jupiter.api.Test;

import java.util.LinkedHashMap;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * End-to-end optimization: minimum servers of an M/M/c queue under a 50%
 * utilization SLA. lambda=3, mu=1 -> optimum c*=6 (cost 60). Matches the native
 * Python {@code opt_server_sizing.py} reference (6 servers, cost 60,
 * 10 evaluations at seed 42).
 */
public class ServerSizingE2ETest {

    @Test
    public void testServerSizing() {
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

        OptimizationResult result = problem.solve(new LineOptSolverOptions().setSeed(42));

        assertEquals(6, ((Number) result.getVariableValue("Server_servers")).intValue(),
                "optimal servers");
        assertEquals(60.0, result.objectiveValue, 1e-9, "cost");
        assertTrue(result.feasible, "feasible");
    }
}
