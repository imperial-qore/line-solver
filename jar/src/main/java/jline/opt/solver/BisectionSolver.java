package jline.opt.solver;

import jline.lang.Network;
import jline.opt.LineEvaluator;
import jline.opt.OptimizationProblem;
import jline.opt.objectives.Constraint;
import jline.opt.objectives.Objective;
import jline.opt.results.EvaluationResult;
import jline.opt.results.OptimizationResult;
import jline.opt.variables.DecisionVariable;
import jline.util.Pair;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Exact O(log n) solver for a single integer decision variable with monotone
 * feasibility (the standard sizing pattern), an alternative to differential
 * evolution. {@code direction='min_feasible'} finds the smallest feasible value
 * (server sizing); {@code 'max_feasible'} the largest (population sizing).
 * Mirrors native-Python {@code line_solver.opt.sizing.BisectionSolver}.
 */
public class BisectionSolver {

    private final OptimizationProblem problem;
    private final String direction;
    private final DecisionVariable variable;
    private final int lo;
    private final int hi;
    private final Map<String, Object> fixedValueDict = new LinkedHashMap<String, Object>();
    private final List<LineEvaluator> evaluators = new ArrayList<LineEvaluator>();
    private final Map<Integer, ProbeOutcome> probeCache = new LinkedHashMap<Integer, ProbeOutcome>();
    private EvaluationResult baseResultAt = null;
    private Map<String, Double> violationsAt = new LinkedHashMap<String, Double>();

    private static class ProbeOutcome {
        final boolean feasible;
        final Map<String, Double> violations;
        final EvaluationResult baseResult;

        ProbeOutcome(boolean feasible, Map<String, Double> violations, EvaluationResult baseResult) {
            this.feasible = feasible;
            this.violations = violations;
            this.baseResult = baseResult;
        }
    }

    public BisectionSolver(OptimizationProblem problem) {
        this(problem, "min_feasible");
    }

    public BisectionSolver(OptimizationProblem problem, String direction) {
        this.problem = problem;
        this.direction = direction;
        List<DecisionVariable> variables = problem.getVariables();
        if (variables.size() != 1 || variables.get(0).getDimension() != 1) {
            throw new IllegalArgumentException(
                    "BisectionSolver requires exactly one dimension-1 decision variable");
        }
        this.variable = variables.get(0);
        Object loObj = variable.decode(new double[]{0.0});
        Object hiObj = variable.decode(new double[]{1.0});
        if (!(loObj instanceof Integer) || !(hiObj instanceof Integer)) {
            throw new IllegalArgumentException(
                    "BisectionSolver requires an integer-valued variable");
        }
        this.lo = (Integer) loObj;
        this.hi = (Integer) hiObj;

        List<Pair<DecisionVariable, Object>> fixed = problem.getFixedVariables();
        for (Pair<DecisionVariable, Object> fv : fixed) {
            fixedValueDict.put(fv.getLeft().getName(), fv.getRight());
        }
        evaluators.add(new LineEvaluator(problem.getModel(), variables, fixed));
        for (Pair<Network, Double> s : problem.getScenarios()) {
            evaluators.add(new LineEvaluator(s.getLeft(), variables, fixed));
        }
    }

    private List<Constraint> allConstraints() {
        Objective objective = problem.getObjective();
        List<Constraint> constraints = new ArrayList<Constraint>();
        if (objective != null) {
            constraints.addAll(objective.getConstraints());
        }
        constraints.addAll(problem.getConstraints());
        return constraints;
    }

    /** Probe feasibility at an integer value; records base result + violations. */
    private boolean probe(int value) {
        ProbeOutcome cached = probeCache.get(value);
        if (cached != null) {
            this.baseResultAt = cached.baseResult;
            this.violationsAt = cached.violations;
            return cached.feasible;
        }
        Map<String, Object> values = new LinkedHashMap<String, Object>();
        values.put(variable.getName(), value);
        Map<String, Object> allValues = new LinkedHashMap<String, Object>(fixedValueDict);
        allValues.putAll(values);
        List<Constraint> constraints = allConstraints();

        boolean feasible = true;
        Map<String, Double> violations = new LinkedHashMap<String, Double>();
        EvaluationResult baseResult = null;
        for (LineEvaluator evaluator : evaluators) {
            EvaluationResult res = evaluator.evaluateValues(values);
            if (baseResult == null) {
                baseResult = res;
            }
            if (!res.feasible) {
                feasible = false;
                continue;
            }
            for (Constraint c : constraints) {
                double v = c.evaluate(res, allValues);
                if (v > 0) {
                    feasible = false;
                    String name = c.getName();
                    double prev = violations.containsKey(name) ? violations.get(name) : 0.0;
                    violations.put(name, Math.max(v, prev));
                }
            }
        }
        this.baseResultAt = baseResult;
        this.violationsAt = violations;
        probeCache.put(value, new ProbeOutcome(feasible, violations, baseResult));
        return feasible;
    }

    public OptimizationResult solve() {
        long start = System.nanoTime();
        int loV = lo;
        int hiV = hi;
        int iterations = 0;

        if ("min_feasible".equals(direction)) {
            while (loV < hiV) {
                int mid = (loV + hiV) / 2;
                boolean feasible = probe(mid);
                iterations++;
                if (feasible) {
                    hiV = mid;
                } else {
                    loV = mid + 1;
                }
            }
        } else if ("max_feasible".equals(direction)) {
            while (loV < hiV) {
                int mid = (loV + hiV + 1) / 2;
                boolean feasible = probe(mid);
                iterations++;
                if (feasible) {
                    loV = mid;
                } else {
                    hiV = mid - 1;
                }
            }
        } else {
            throw new IllegalArgumentException("Unknown direction: " + direction);
        }

        int chosen = loV;
        boolean feasible = probe(chosen);

        OptimizationResult result = new OptimizationResult();
        result.variableValues.put(variable.getName(), chosen);
        result.feasible = feasible;
        result.constraintViolations.putAll(violationsAt);
        result.iterations = iterations;
        int evals = 0;
        for (LineEvaluator e : evaluators) {
            evals += e.getEvaluationCount();
        }
        result.modelEvaluations = evals;
        result.terminatedBy = "bisection";

        Objective objective = problem.getObjective();
        if (objective != null && baseResultAt != null) {
            Map<String, Object> allValues = new LinkedHashMap<String, Object>(fixedValueDict);
            allValues.putAll(result.variableValues);
            result.objectiveValue = objective.evaluate(baseResultAt, allValues);
        }
        result.solveTime = (System.nanoTime() - start) / 1e9;
        return result;
    }
}
