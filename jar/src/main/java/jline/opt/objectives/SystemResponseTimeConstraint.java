package jline.opt.objectives;

import jline.lang.JobClass;
import jline.opt.results.EvaluationResult;

import java.util.Map;

/** End-to-end (chain-level) response time constraint: SysRespT &lt;= maxValue. */
public class SystemResponseTimeConstraint extends Constraint {

    private final String jobclass;   // null for throughput-weighted aggregate
    private final double maxValue;

    public SystemResponseTimeConstraint(String jobclass, double maxValue, String name) {
        super(name);
        this.jobclass = jobclass;
        this.maxValue = maxValue;
    }

    public SystemResponseTimeConstraint(JobClass jobclass, double maxValue) {
        this(jobclass == null ? null : jobclass.getName(), maxValue, null);
    }

    public SystemResponseTimeConstraint(double maxValue) {
        this((String) null, maxValue, null);
    }

    protected String generateName() {
        if (jobclass != null) {
            return "SysRT_" + jobclass + "_le_" + maxValue;
        }
        return "SysRT_le_" + maxValue;
    }

    public String getJobClass() {
        return jobclass;
    }

    public double getMaxValue() {
        return maxValue;
    }

    public double evaluate(EvaluationResult result, Map<String, Object> variableValues) {
        double actual = (jobclass != null)
                ? result.getSystemResponseTime(jobclass)
                : result.getSystemResponseTime();
        return upperBoundViolation(actual, maxValue);
    }
}
