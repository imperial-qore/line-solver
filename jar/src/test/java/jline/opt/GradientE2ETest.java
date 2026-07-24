package jline.opt;

import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.opt.objectives.MinimizeCost;
import jline.opt.objectives.ResponseTimeConstraint;
import jline.opt.results.EvaluationResult;
import jline.opt.results.OptimizationResult;
import jline.opt.results.SensitivityData;
import jline.opt.sensitivity.SensitivityTable;
import jline.opt.solver.LineOptSolverOptions;
import jline.opt.variables.DecisionVariable;
import jline.opt.variables.ServiceRate;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Validates the analytic sensitivity + gradient optimizer: the cheapest service
 * rate meeting an RT SLA on an M/M/1 (theory rate = 5.0), and the exact
 * open-network d(RespT)/d(rate) = -1 at mu=4, lambda=3.
 */
public class GradientE2ETest {

    private Network mm1(double mu) {
        Network model = new Network("MM1");
        Source source = new Source(model, "Arrivals");
        Queue queue = new Queue(model, "Server", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Departures");
        OpenClass jobs = new OpenClass(model, "Jobs");
        source.setArrival(jobs, new Exp(3.0));
        queue.setService(jobs, new Exp(mu));
        model.link(Network.serialRouting(source, queue, sink));
        return model;
    }

    @Test
    public void testOpenSensitivity() {
        Network model = mm1(4.0);
        SensitivityData sens = SensitivityTable.compute(model);
        assertNotNull(sens, "open product-form sensitivity should be available");
        Map<String, Map<String, Double>> rt = sens.forKind("RespT");
        String mkey = SensitivityData.metricKey("Server", "Jobs");
        String pkey = SensitivityData.paramKey("Server", "Jobs");
        // d[1/(mu-lambda)]/dmu = -1/(mu-lambda)^2 = -1 at mu=4, lambda=3
        assertEquals(-1.0, rt.get(mkey).get(pkey), 1e-6);
    }

    @Test
    public void testGradientServiceRate() {
        Network model = mm1(4.0);
        Queue queue = (Queue) model.getNodes().get(1);
        OpenClass jobs = (OpenClass) model.getClasses().get(0);

        OptimizationProblem problem = new OptimizationProblem(model);
        Map<jline.lang.nodes.Station, Double> rateCost =
                new LinkedHashMap<jline.lang.nodes.Station, Double>();
        rateCost.put(queue, 20.0);
        problem.addVariable(new ServiceRate(queue, jobs, 3.5, 8.0));
        problem.setObjective(new MinimizeCost(null, rateCost, null, null));
        problem.addConstraint(new ResponseTimeConstraint(queue, jobs, 0.5));

        OptimizationResult result = problem.solve(
                new LineOptSolverOptions().setSeed(42).setOptimizer("gradient").setMaxIterations(60));
        double rate = ((Number) result.getVariableValue("Server_Jobs_rate")).doubleValue();
        assertEquals(5.0, rate, 0.1, "gradient optimum (theory 5.0)");
        assertTrue(result.feasible);
    }
}
