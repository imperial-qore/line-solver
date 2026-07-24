package jline.opt.solver;

import java.util.ArrayList;
import java.util.List;

/**
 * Configuration options for {@link LineOptSolver}, mirroring the option dict of
 * native-Python {@code LineOptSolver.defaultOptions()}. Chained setters return
 * {@code this}.
 */
public class LineOptSolverOptions {

    public String strategy = "best1bin";
    public int popsize = 15;
    public double mutationLow = 0.5;
    public double mutationHigh = 1.0;
    public double recombination = 0.7;
    public double tol = 0.01;
    public int maxIterations = 100;
    public double timeLimit = 300.0;
    public Long seed = null;
    public boolean verbose = false;
    public double penaltyWeight = 1e6;
    /** 'worst' (minimax) or 'mean' (weighted) aggregation across scenarios. */
    public String scenarioAggregation = "worst";
    /** 'evolution', 'gradient', or 'auto'. */
    public String optimizer = "evolution";
    public double fdStep = 1e-6;
    public int gradientRestarts = 4;
    /**
     * LayeredNetwork (LQN) gradient source, used only when the model is a
     * LayeredNetwork and the gradient path is taken: 'fd' (whole-model finite
     * difference, robust default), 'partial_sens' (SolverLN per-layer partial
     * derivatives, cheap/biased), or 'partial_plus_fd' (partial with a periodic
     * full-FD correction every {@link #fdRefresh} gradient evaluations).
     */
    public String lqnGradient = "fd";
    /** partial_plus_fd full finite-difference correction period. */
    public int fdRefresh = 5;
    /**
     * Explicit layer freezing (LQN only): host/task layer names whose variables
     * are held at the model's current value instead of being optimized. Null (or
     * empty) freezes nothing.
     */
    public List<String> frozenLayers = null;

    public LineOptSolverOptions setStrategy(String v) {
        this.strategy = v;
        return this;
    }

    public LineOptSolverOptions setPopsize(int v) {
        this.popsize = v;
        return this;
    }

    public LineOptSolverOptions setMutation(double low, double high) {
        this.mutationLow = low;
        this.mutationHigh = high;
        return this;
    }

    public LineOptSolverOptions setRecombination(double v) {
        this.recombination = v;
        return this;
    }

    public LineOptSolverOptions setTolerance(double v) {
        this.tol = v;
        return this;
    }

    public LineOptSolverOptions setMaxIterations(int v) {
        this.maxIterations = v;
        return this;
    }

    public LineOptSolverOptions setTimeLimit(double v) {
        this.timeLimit = v;
        return this;
    }

    public LineOptSolverOptions setSeed(long v) {
        this.seed = v;
        return this;
    }

    public LineOptSolverOptions setVerbose(boolean v) {
        this.verbose = v;
        return this;
    }

    public LineOptSolverOptions setPenaltyWeight(double v) {
        this.penaltyWeight = v;
        return this;
    }

    public LineOptSolverOptions setScenarioAggregation(String v) {
        this.scenarioAggregation = v;
        return this;
    }

    public LineOptSolverOptions setOptimizer(String v) {
        this.optimizer = v;
        return this;
    }

    public LineOptSolverOptions setFdStep(double v) {
        this.fdStep = v;
        return this;
    }

    public LineOptSolverOptions setGradientRestarts(int v) {
        this.gradientRestarts = v;
        return this;
    }

    /**
     * Set the LQN gradient source: 'fd', 'partial_sens', or 'partial_plus_fd'.
     */
    public LineOptSolverOptions setLqnGradient(String v) {
        if (!"fd".equals(v) && !"partial_sens".equals(v) && !"partial_plus_fd".equals(v)) {
            throw new IllegalArgumentException(
                    "lqn_gradient must be 'fd', 'partial_sens', or 'partial_plus_fd'");
        }
        this.lqnGradient = v;
        return this;
    }

    /** Set the partial_plus_fd full finite-difference correction period. */
    public LineOptSolverOptions setFdRefresh(int v) {
        this.fdRefresh = v;
        return this;
    }

    /** Freeze the given LQN layers (list of host/task layer names). */
    public LineOptSolverOptions setFrozenLayers(List<String> layers) {
        this.frozenLayers = (layers != null && !layers.isEmpty())
                ? new ArrayList<String>(layers) : null;
        return this;
    }
}
