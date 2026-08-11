package jline.opt;

import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.opt.objectives.Constraint;
import jline.opt.objectives.MinimizeCost;
import jline.opt.objectives.UtilizationConstraint;
import jline.opt.pareto.ParetoPoint;
import jline.opt.pareto.ParetoSweep;
import jline.opt.variables.ServerAllocation;
import org.junit.jupiter.api.Test;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Epsilon-constraint cost frontier over a utilization SLA, matching the native
 * Python {@code opt_pareto_frontier.py} reference:
 * 0.3->100, 0.4->80, 0.5->60, 0.6->50, 0.75->40 (0.9 dominated).
 */
public class ParetoE2ETest {

    @Test
    public void testParetoFrontier() {
        Network model = new Network("MMc");
        Source source = new Source(model, "Arrivals");
        final Queue queue = new Queue(model, "Server", SchedStrategy.FCFS);
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

        double[] epsilons = {0.3, 0.4, 0.5, 0.6, 0.75, 0.9};
        ParetoSweep sweep = new ParetoSweep(problem, new ParetoSweep.ConstraintFactory() {
            public Constraint make(double epsilon) {
                return new UtilizationConstraint(queue, epsilon);
            }
        }, epsilons, "bisection");
        sweep.solve();

        List<ParetoPoint> frontier = sweep.getFrontier();
        double[][] expected = {{0.3, 100}, {0.4, 80}, {0.5, 60}, {0.6, 50}, {0.75, 40}};
        assertEquals(expected.length, frontier.size(), "frontier size");
        for (int i = 0; i < expected.length; i++) {
            assertEquals(expected[i][0], frontier.get(i).epsilon, 1e-9, "epsilon " + i);
            assertEquals(expected[i][1], frontier.get(i).objectiveValue, 1e-9, "cost " + i);
        }
    }
}
