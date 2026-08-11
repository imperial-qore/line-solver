package jline.opt;

import jline.opt.objectives.MaximizePerformance;
import jline.opt.objectives.MinimizeCost;
import jline.opt.objectives.ResponseTimeConstraint;
import jline.opt.objectives.ThroughputConstraint;
import jline.opt.objectives.UtilizationConstraint;
import jline.opt.results.EvaluationResult;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.assertEquals;

/** Verifies constraint violations and objective/penalty arithmetic. */
public class ObjectivesTest {

    private EvaluationResult sample() {
        EvaluationResult er = new EvaluationResult();
        er.setResponseTime("Server", "Jobs", 0.8);
        er.setThroughput("Server", "Jobs", 5.0);
        er.utilizations.put("Server", 0.75);
        return er;
    }

    @Test
    public void testConstraints() {
        EvaluationResult er = sample();
        Map<String, Object> vv = new LinkedHashMap<String, Object>();
        vv.put("Server_servers", 3);
        assertEquals(0.3, new ResponseTimeConstraint("Server", "Jobs", 0.5, null).evaluate(er, vv), 1e-12);
        assertEquals(0.0, new UtilizationConstraint("Server", 0.9, null).evaluate(er, vv), 1e-12);
        assertEquals(5.0, new ThroughputConstraint("Server", null, 10.0, null).evaluate(er, vv), 1e-12);
    }

    @Test
    public void testNonFiniteMetricIsInfeasible() {
        EvaluationResult er = new EvaluationResult();
        // no response time recorded -> getResponseTime returns +inf -> strongly infeasible
        Map<String, Object> vv = new LinkedHashMap<String, Object>();
        double viol = new ResponseTimeConstraint("Server", "Jobs", 0.5, null).evaluate(er, vv);
        assertEquals(Double.POSITIVE_INFINITY, viol, 0.0);
    }

    @Test
    public void testMinimizeCostWithPenalty() {
        EvaluationResult er = sample();
        Map<String, Object> vv = new LinkedHashMap<String, Object>();
        vv.put("Server_servers", 3);
        Map<jline.lang.nodes.Station, Double> serverCost =
                new LinkedHashMap<jline.lang.nodes.Station, Double>();
        // build a station named "Server"
        jline.lang.Network m = new jline.lang.Network("m");
        jline.lang.nodes.Queue q = new jline.lang.nodes.Queue(m, "Server",
                jline.lang.constant.SchedStrategy.FCFS);
        serverCost.put(q, 100.0);
        List<jline.opt.objectives.Constraint> subj = new ArrayList<jline.opt.objectives.Constraint>();
        subj.add(new ResponseTimeConstraint("Server", "Jobs", 0.5, null));
        MinimizeCost mc = new MinimizeCost(serverCost, null, null, subj);
        assertEquals(300.0, mc.evaluate(er, vv), 1e-9);
        assertEquals(300.0 + 0.3 * 1e6, mc.evaluateWithPenalty(er, vv, 1e6), 1e-3);
    }

    @Test
    public void testMaximizePerformance() {
        EvaluationResult er = sample();
        Map<String, Object> vv = new LinkedHashMap<String, Object>();
        MaximizePerformance mp = new MaximizePerformance(1.0, 2.0);
        // -(5 + 2*(1/0.8)) = -7.5
        assertEquals(-7.5, mp.evaluate(er, vv), 1e-12);
    }
}
