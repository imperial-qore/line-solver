package jline.opt.objectives;

import jline.opt.results.EvaluationResult;

import java.util.Map;

/**
 * Abstract base for line-opt constraints. Each constraint computes a
 * non-negative violation amount (0 if satisfied) given an
 * {@link EvaluationResult} and the current variable values. Mirrors
 * native-Python {@code line_solver.opt.objectives.Constraint}.
 */
public abstract class Constraint {

    protected final String name;

    protected Constraint(String name) {
        this.name = (name != null) ? name : "";
    }

    public String getName() {
        return name.isEmpty() ? generateName() : name;
    }

    protected abstract String generateName();

    /** Violation amount: 0 if satisfied, positive if violated. */
    public abstract double evaluate(EvaluationResult result, Map<String, Object> variableValues);

    public boolean isSatisfied(EvaluationResult result, Map<String, Object> variableValues,
                               double tolerance) {
        return evaluate(result, variableValues) <= tolerance;
    }

    /**
     * Violation for an upper-bound constraint (actual &lt;= bound). A non-finite
     * metric signals divergence and is reported as strongly infeasible (+inf),
     * mirroring the native guard that avoids NaN silently satisfying bounds.
     */
    protected static double upperBoundViolation(double actual, double bound) {
        if (!isFiniteNumber(actual)) {
            return Double.POSITIVE_INFINITY;
        }
        return Math.max(0.0, actual - bound);
    }

    /** Violation for a lower-bound constraint (actual &gt;= bound). */
    protected static double lowerBoundViolation(double actual, double bound) {
        if (!isFiniteNumber(actual)) {
            return Double.POSITIVE_INFINITY;
        }
        return Math.max(0.0, bound - actual);
    }

    private static boolean isFiniteNumber(double v) {
        return !Double.isNaN(v) && !Double.isInfinite(v);
    }
}
