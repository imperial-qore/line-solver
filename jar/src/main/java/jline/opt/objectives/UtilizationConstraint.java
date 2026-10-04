package jline.opt.objectives;

import jline.lang.nodes.Station;
import jline.opt.results.EvaluationResult;

import java.util.Map;

/** Utilization constraint: U &lt;= maxValue. */
public class UtilizationConstraint extends Constraint {

    private final String station;
    private final double maxValue;

    public UtilizationConstraint(String station, double maxValue, String name) {
        super(name);
        this.station = station;
        this.maxValue = maxValue;
    }

    public UtilizationConstraint(Station station, double maxValue) {
        this(station.getName(), maxValue, null);
    }

    protected String generateName() {
        return "Util_" + station + "_le_" + maxValue;
    }

    public String getStation() {
        return station;
    }

    public double getMaxValue() {
        return maxValue;
    }

    public double evaluate(EvaluationResult result, Map<String, Object> variableValues) {
        return upperBoundViolation(result.getUtilization(station), maxValue);
    }
}
