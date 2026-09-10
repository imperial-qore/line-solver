package jline.opt.objectives;

import jline.lang.JobClass;
import jline.lang.nodes.Station;
import jline.opt.results.EvaluationResult;

import java.util.Map;

/** Per-station response time constraint: RT &lt;= maxValue. */
public class ResponseTimeConstraint extends Constraint {

    private final String station;
    private final String jobclass;   // null for aggregate
    private final double maxValue;

    public ResponseTimeConstraint(String station, String jobclass, double maxValue, String name) {
        super(name);
        this.station = station;
        this.jobclass = jobclass;
        this.maxValue = maxValue;
    }

    public ResponseTimeConstraint(Station station, JobClass jobclass, double maxValue) {
        this(station.getName(), jobclass == null ? null : jobclass.getName(), maxValue, null);
    }

    public ResponseTimeConstraint(Station station, double maxValue) {
        this(station.getName(), null, maxValue, null);
    }

    protected String generateName() {
        if (jobclass != null) {
            return "RT_" + station + "_" + jobclass + "_le_" + maxValue;
        }
        return "RT_" + station + "_le_" + maxValue;
    }

    public String getStation() {
        return station;
    }

    public String getJobClass() {
        return jobclass;
    }

    public double getMaxValue() {
        return maxValue;
    }

    public double evaluate(EvaluationResult result, Map<String, Object> variableValues) {
        double actual = (jobclass != null)
                ? result.getResponseTime(station, jobclass)
                : result.getResponseTime(station);
        return upperBoundViolation(actual, maxValue);
    }
}
