package jline.opt.objectives;

import jline.opt.results.EvaluationResult;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Abstract base for line-opt objectives. An objective defines the scalar to
 * minimize, optionally with attached constraints folded in as penalties by
 * {@link #evaluateWithPenalty}. Mirrors native-Python
 * {@code line_solver.opt.objectives.Objective}.
 */
public abstract class Objective {

    protected List<Constraint> constraints = new ArrayList<Constraint>();

    public List<Constraint> getConstraints() {
        return constraints;
    }

    /** The scalar to minimize (before constraint penalties). */
    public abstract double evaluate(EvaluationResult result, Map<String, Object> variableValues);

    public abstract boolean isMinimization();

    /** Objective value plus penalty_weight times total constraint violation. */
    public double evaluateWithPenalty(EvaluationResult result,
                                      Map<String, Object> variableValues,
                                      double penaltyWeight) {
        double obj = evaluate(result, variableValues);
        double penalty = 0.0;
        for (Constraint c : constraints) {
            penalty += c.evaluate(result, variableValues) * penaltyWeight;
        }
        return obj + penalty;
    }
}
