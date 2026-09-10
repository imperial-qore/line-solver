package jline.opt.objectives;

import jline.opt.results.EvaluationResult;

import java.util.LinkedHashMap;
import java.util.Map;

/** Budget constraint: total cost &lt;= budget, cost from variable values. */
public class BudgetConstraint extends Constraint {

    private final double budget;
    private final Map<String, Double> costCoefficients;

    public BudgetConstraint(double budget, Map<String, Double> costCoefficients, String name) {
        super(name);
        this.budget = budget;
        this.costCoefficients = (costCoefficients != null)
                ? costCoefficients : new LinkedHashMap<String, Double>();
    }

    public BudgetConstraint(double budget, Map<String, Double> costCoefficients) {
        this(budget, costCoefficients, null);
    }

    protected String generateName() {
        return "Budget_le_" + budget;
    }

    public double getBudget() {
        return budget;
    }

    public Map<String, Double> getCostCoefficients() {
        return costCoefficients;
    }

    public double computeCost(Map<String, Object> variableValues) {
        double total = 0.0;
        for (Map.Entry<String, Double> e : costCoefficients.entrySet()) {
            Object value = variableValues.get(e.getKey());
            total += e.getValue() * ObjectiveUtil.numericValue(value);
        }
        return total;
    }

    public double evaluate(EvaluationResult result, Map<String, Object> variableValues) {
        return Math.max(0.0, computeCost(variableValues) - budget);
    }
}
