package jline.opt;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Node;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.opt.objectives.MinimizeCost;
import jline.opt.objectives.MinimizeSystemResponseTime;
import jline.opt.objectives.ResponseTimeConstraint;
import jline.opt.objectives.UtilizationConstraint;
import jline.opt.results.OptimizationResult;
import jline.opt.solver.BisectionSolver;
import jline.opt.solver.LineOptSolverOptions;
import jline.opt.variables.JobPopulation;
import jline.opt.variables.RoutingProbabilities;
import jline.opt.variables.ServerAllocation;
import jline.opt.variables.ServiceRate;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * End-to-end validation of the remaining native line-opt examples against their
 * Python reference outputs, exercising every solver and variable type.
 */
public class ExamplesE2ETest {

    /** opt_population_sizing: max users under an RT SLA. Native: 5 users. */
    @Test
    public void testPopulationSizing() {
        Network model = new Network("Interactive");
        Delay think = new Delay(model, "Think");
        Queue server = new Queue(model, "AppServer", SchedStrategy.PS);
        ClosedClass users = new ClosedClass(model, "Users", 1, think);
        think.setService(users, new Exp(1.0));
        server.setService(users, new Exp(2.0));
        model.link(Network.serialRouting(think, server));

        OptimizationProblem problem = new OptimizationProblem(model);
        problem.addVariable(new JobPopulation(users, 1, 50));
        problem.setObjective(new MinimizeCost(null, null, null, null));
        problem.addConstraint(new ResponseTimeConstraint(server, users, 2.0));

        OptimizationResult result = new BisectionSolver(problem, "max_feasible").solve();
        assertEquals(5, ((Number) result.getVariableValue("Users_population")).intValue());
        assertTrue(result.feasible);
    }

    /** opt_service_rate: cheapest rate under an RT SLA. Native: rate ~5.001. */
    @Test
    public void testServiceRate() {
        Network model = new Network("MM1");
        Source source = new Source(model, "Arrivals");
        Queue queue = new Queue(model, "Server", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Departures");
        OpenClass jobs = new OpenClass(model, "Jobs");
        source.setArrival(jobs, new Exp(3.0));
        queue.setService(jobs, new Exp(4.0));
        model.link(Network.serialRouting(source, queue, sink));

        OptimizationProblem problem = new OptimizationProblem(model);
        Map<jline.lang.nodes.Station, Double> rateCost =
                new LinkedHashMap<jline.lang.nodes.Station, Double>();
        rateCost.put(queue, 20.0);
        problem.addVariable(new ServiceRate(queue, jobs, 3.5, 8.0));
        problem.setObjective(new MinimizeCost(null, rateCost, null, null));
        problem.addConstraint(new ResponseTimeConstraint(queue, jobs, 0.5));

        OptimizationResult result = problem.solve(
                new LineOptSolverOptions().setSeed(42).setMaxIterations(60));
        double rate = ((Number) result.getVariableValue("Server_Jobs_rate")).doubleValue();
        assertEquals(5.0, rate, 0.05, "optimal rate (theory 5.0)");
        assertTrue(result.feasible);
    }

    /** opt_load_balancing: routing split minimizing system RT. Native fast ~0.736. */
    @Test
    public void testLoadBalancing() {
        Network model = new Network("LoadBalance");
        Source source = new Source(model, "S");
        Queue fast = new Queue(model, "Fast", SchedStrategy.PS);
        Queue slow = new Queue(model, "Slow", SchedStrategy.PS);
        Sink sink = new Sink(model, "K");
        OpenClass jobs = new OpenClass(model, "Jobs");
        source.setArrival(jobs, new Exp(2.0));
        fast.setService(jobs, new Exp(4.0));
        slow.setService(jobs, new Exp(2.5));
        model.addLink(source, fast);
        model.addLink(source, slow);
        model.addLink(fast, sink);
        model.addLink(slow, sink);

        OptimizationProblem problem = new OptimizationProblem(model);
        List<Node> targets = new ArrayList<Node>();
        targets.add(fast);
        targets.add(slow);
        problem.addVariable(new RoutingProbabilities(jobs, source, targets));
        problem.setObjective(new MinimizeSystemResponseTime(jobs));

        OptimizationResult result = problem.solve(
                new LineOptSolverOptions().setSeed(42).setMaxIterations(60));
        double[] probs = (double[]) result.getVariableValue("Jobs_routing_from_S");
        assertEquals(0.736, probs[0], 0.02, "fraction to fast (theory 0.743)");
        assertEquals(1.0, probs[0] + probs[1], 1e-9);
    }

    /** opt_robust_sizing: sizing under a peak scenario. Native: 9 servers. */
    @Test
    public void testRobustSizing() {
        Network base = buildMMc(3.0);
        Network peak = buildMMc(4.5);
        Queue queue = (Queue) base.getNodes().get(1);

        OptimizationProblem problem = new OptimizationProblem(base);
        Map<jline.lang.nodes.Station, Double> serverCost =
                new LinkedHashMap<jline.lang.nodes.Station, Double>();
        serverCost.put(queue, 10.0);
        problem.addVariable(new ServerAllocation(queue, 1, 12));
        problem.setObjective(new MinimizeCost(serverCost, null, null, null));
        problem.addConstraint(new UtilizationConstraint(queue, 0.5));
        problem.addScenario(peak, 1.0);

        OptimizationResult result = new BisectionSolver(problem).solve();
        assertEquals(9, ((Number) result.getVariableValue("Server_servers")).intValue());
        assertEquals(90.0, result.objectiveValue, 1e-9);
        assertTrue(result.feasible);
    }

    private Network buildMMc(double arrivalRate) {
        Network model = new Network("MMc");
        Source source = new Source(model, "Arrivals");
        Queue queue = new Queue(model, "Server", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Departures");
        OpenClass jobs = new OpenClass(model, "Jobs");
        source.setArrival(jobs, new Exp(arrivalRate));
        queue.setService(jobs, new Exp(1.0));
        model.link(Network.serialRouting(source, queue, sink));
        return model;
    }
}
