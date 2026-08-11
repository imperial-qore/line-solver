package jline.opt.pareto;

import jline.opt.results.OptimizationResult;

/**
 * One point of a cost-performance tradeoff curve produced by {@link ParetoSweep}.
 * Mirrors native-Python {@code line_solver.opt.pareto.ParetoPoint}.
 */
public class ParetoPoint {

    public double epsilon;
    public double objectiveValue;
    public boolean feasible;
    public OptimizationResult result;

    public ParetoPoint(double epsilon, double objectiveValue, boolean feasible,
                       OptimizationResult result) {
        this.epsilon = epsilon;
        this.objectiveValue = objectiveValue;
        this.feasible = feasible;
        this.result = result;
    }
}
