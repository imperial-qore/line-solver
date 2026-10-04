package jline.opt.objectives;

import jline.lang.JobClass;
import jline.lang.nodes.Station;
import jline.opt.results.EvaluationResult;

import java.util.Map;

/** Throughput constraint: Tput &gt;= minValue. */
public class ThroughputConstraint extends Constraint {

    private final String station;
    private final String jobclass;   // null for aggregate
    private final double minValue;

    public ThroughputConstraint(String station, String jobclass, double minValue, String name) {
        super(name);
        this.station = station;
        this.jobclass = jobclass;
        this.minValue = minValue;
    }

    public ThroughputConstraint(Station station, JobClass jobclass, double minValue) {
        this(station.getName(), jobclass == null ? null : jobclass.getName(), minValue, null);
    }

    public ThroughputConstraint(Station station, double minValue) {
        this(station.getName(), null, minValue, null);
    }

    protected String generateName() {
        if (jobclass != null) {
            return "Tput_" + station + "_" + jobclass + "_ge_" + minValue;
        }
        return "Tput_" + station + "_ge_" + minValue;
    }

    public String getStation() {
        return station;
    }

    public double getMinValue() {
        return minValue;
    }

    public double evaluate(EvaluationResult result, Map<String, Object> variableValues) {
        double actual = (jobclass != null)
                ? result.getThroughput(station, jobclass)
                : result.getThroughput(station);
        return lowerBoundViolation(actual, minValue);
    }
}
