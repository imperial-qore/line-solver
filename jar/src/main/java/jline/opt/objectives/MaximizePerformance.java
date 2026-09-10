package jline.opt.objectives;

import jline.lang.nodes.Station;
import jline.opt.results.EvaluationResult;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Maximize a weighted combination of throughput, 1/response-time and
 * 1/queue-length, subject to an optional budget constraint. Because the DE
 * solver minimizes, {@link #evaluate} returns the negated performance. Mirrors
 * native-Python {@code MaximizePerformance}.
 */
public class MaximizePerformance extends Objective {

    private final double throughputWeight;
    private final double responseTimeWeight;
    private final double queueLengthWeight;
    private final List<String> stations;   // null = all stations in result

    public MaximizePerformance(double throughputWeight, double responseTimeWeight,
                               double queueLengthWeight, List<Station> stations,
                               Double budget, Map<String, Double> budgetTerms) {
        this.throughputWeight = throughputWeight;
        this.responseTimeWeight = responseTimeWeight;
        this.queueLengthWeight = queueLengthWeight;
        if (stations != null) {
            this.stations = new ArrayList<String>();
            for (Station s : stations) {
                this.stations.add(s.getName());
            }
        } else {
            this.stations = null;
        }
        if (budget != null) {
            this.constraints.add(new BudgetConstraint(budget,
                    budgetTerms != null ? budgetTerms : new LinkedHashMap<String, Double>()));
        }
    }

    public MaximizePerformance(double throughputWeight, double responseTimeWeight) {
        this(throughputWeight, responseTimeWeight, 0.0, null, null, null);
    }

    public boolean isMinimization() {
        return true;
    }

    public double evaluate(EvaluationResult result, Map<String, Object> variableValues) {
        double performance = 0.0;
        List<String> sts = stations;
        if (sts == null) {
            sts = new ArrayList<String>(result.throughputs.keySet());
        }
        for (String station : sts) {
            if (throughputWeight > 0) {
                performance += throughputWeight * result.getThroughput(station);
            }
            if (responseTimeWeight > 0) {
                double rt = result.getResponseTime(station);
                if (rt > 0 && !Double.isInfinite(rt)) {
                    performance += responseTimeWeight * (1.0 / rt);
                }
            }
            if (queueLengthWeight > 0) {
                double qlen = result.getQueueLength(station);
                if (qlen > 0) {
                    performance += queueLengthWeight * (1.0 / qlen);
                }
            }
        }
        return -performance;
    }
}
