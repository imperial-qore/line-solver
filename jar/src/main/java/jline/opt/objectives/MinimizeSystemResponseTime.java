package jline.opt.objectives;

import jline.lang.JobClass;
import jline.opt.results.EvaluationResult;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Minimize the end-to-end (system) response time. Typical uses are load
 * balancing and routing optimization. Mirrors native-Python
 * {@code MinimizeSystemResponseTime}.
 */
public class MinimizeSystemResponseTime extends Objective {

    private final String jobclass;   // null = throughput-weighted aggregate

    public MinimizeSystemResponseTime(String jobclass, List<Constraint> subjectTo) {
        this.jobclass = jobclass;
        if (subjectTo != null) {
            this.constraints = new ArrayList<Constraint>(subjectTo);
        }
    }

    public MinimizeSystemResponseTime(JobClass jobclass) {
        this(jobclass == null ? null : jobclass.getName(), null);
    }

    public MinimizeSystemResponseTime() {
        this((String) null, null);
    }

    public String getJobClass() {
        return jobclass;
    }

    public boolean isMinimization() {
        return true;
    }

    public double evaluate(EvaluationResult result, Map<String, Object> variableValues) {
        return (jobclass != null)
                ? result.getSystemResponseTime(jobclass)
                : result.getSystemResponseTime();
    }
}
