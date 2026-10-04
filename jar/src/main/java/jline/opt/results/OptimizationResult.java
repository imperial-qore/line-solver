package jline.opt.results;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Result from a single optimization run. Mirrors native-Python
 * {@code line_solver.opt.results.OptimizationResult}: the optimal variable
 * values (by name), objective value, constraint status, and solver statistics.
 * Variable values are stored as {@code Object} because a variable decodes to a
 * scalar (server count, rate, population, priority) or a vector (routing).
 */
public class OptimizationResult {

    public double objectiveValue = Double.POSITIVE_INFINITY;
    public final Map<String, Object> variableValues = new LinkedHashMap<String, Object>();
    public final Map<String, Double> constraintViolations = new LinkedHashMap<String, Double>();
    public boolean feasible = false;
    public int iterations = 0;
    public double solveTime = 0.0;
    public int modelEvaluations = 0;
    public final List<Double> convergenceHistory = new ArrayList<Double>();
    public String terminatedBy = "";

    public boolean isFeasible() {
        return feasible;
    }

    public double getObjectiveValue() {
        return objectiveValue;
    }

    public Object getVariableValue(String name) {
        return variableValues.get(name);
    }

    public double getConstraintViolation(String name) {
        Double v = constraintViolations.get(name);
        return v == null ? 0.0 : v;
    }

    public double getTotalViolation() {
        double sum = 0.0;
        for (double v : constraintViolations.values()) {
            sum += v;
        }
        return sum;
    }

    public String toString() {
        String status = feasible ? "feasible" : "infeasible";
        return String.format("OptimizationResult(obj=%.4f, %s, iters=%d, evals=%d, time=%.2fs)",
                objectiveValue, status, iterations, modelEvaluations, solveTime);
    }
}
