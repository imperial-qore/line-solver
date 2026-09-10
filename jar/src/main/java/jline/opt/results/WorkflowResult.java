package jline.opt.results;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Result from a decomposed workflow optimization: results from each subproblem
 * and overall fixed-point convergence info. Mirrors native-Python
 * {@code line_solver.opt.results.WorkflowResult}.
 */
public class WorkflowResult {

    public double finalObjective = Double.POSITIVE_INFINITY;
    public final Map<String, SubProblemResult> subproblemResults = new LinkedHashMap<String, SubProblemResult>();
    public int cyclesCompleted = 0;
    public boolean converged = false;
    public double totalSolveTime = 0.0;
    public final List<Double> objectiveHistory = new ArrayList<Double>();
    public final Map<String, Object> finalVariableValues = new LinkedHashMap<String, Object>();
    /** LQN diagnostics ({@code solveLayered}): final frozen layer set. */
    public List<String> frozenLayers = new ArrayList<String>();
    /** LQN diagnostics ({@code solveLayered}): total LINE solves performed. */
    public int modelEvaluations = 0;

    public boolean isConverged() {
        return converged;
    }

    public SubProblemResult getSubProblemResult(String name) {
        return subproblemResults.get(name);
    }

    public Object getFinalVariableValue(String name) {
        return finalVariableValues.get(name);
    }

    public String toString() {
        String status = converged ? "converged" : "not converged";
        return String.format("WorkflowResult(obj=%.4f, %d subproblems, %d cycles, %s, time=%.2fs)",
                finalObjective, subproblemResults.size(), cyclesCompleted, status, totalSolveTime);
    }
}
