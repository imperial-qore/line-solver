package jline.opt.decomposition;

import jline.lang.Network;
import jline.lang.layered.LayeredNetwork;
import jline.opt.LineEvaluator;
import jline.opt.OptimizationProblem;
import jline.opt.objectives.Constraint;
import jline.opt.objectives.Objective;
import jline.opt.results.EvaluationResult;
import jline.opt.results.OptimizationResult;
import jline.opt.results.SubProblemResult;
import jline.opt.results.WorkflowResult;
import jline.opt.solver.LineOptSolver;
import jline.opt.solver.LineOptSolverOptions;
import jline.opt.variables.DecisionVariable;
import jline.util.Pair;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Decomposes a joint optimization problem into per-variable-type subproblems
 * solved via Gauss-Seidel cycling with fixed-value propagation. An internal
 * topological sort (no external dependency) orders subproblems when explicit
 * dependencies are set. Mirrors native-Python
 * {@code line_solver.opt.decomposition.DecompositionWorkflow}.
 */
public class DecompositionWorkflow {

    // Flat-network variable types first, then LayeredNetwork (LQN) types; only
    // the types present in a given problem produce subproblems.
    private static final String[] DEFAULT_ORDER = {
            "server_allocation", "station_replicas", "service_rate",
            "job_population", "class_priority", "routing", "class_mapping",
            // LayeredNetwork (LQN) variable types
            "processor_multiplicity", "task_multiplicity", "task_replication",
            "host_demand", "think_time"
    };

    private final OptimizationProblem problem;
    private List<SubProblem> subproblems = new ArrayList<SubProblem>();
    private final Map<String, Set<String>> dependencyGraph = new LinkedHashMap<String, Set<String>>();
    private LineOptSolverOptions solverOptions = new LineOptSolverOptions();

    public DecompositionWorkflow(OptimizationProblem problem) {
        this.problem = problem;
    }

    public OptimizationProblem getProblem() {
        return problem;
    }

    public List<SubProblem> getSubProblems() {
        return subproblems;
    }

    public DecompositionWorkflow setSolverOptions(LineOptSolverOptions options) {
        this.solverOptions = options;
        return this;
    }

    public DecompositionWorkflow autoDecompose() {
        Map<String, List<DecisionVariable>> byType = new LinkedHashMap<String, List<DecisionVariable>>();
        for (DecisionVariable var : problem.getVariables()) {
            String t = var.getVariableType();
            if (!byType.containsKey(t)) {
                byType.put(t, new ArrayList<DecisionVariable>());
            }
            byType.get(t).add(var);
        }
        subproblems = new ArrayList<SubProblem>();
        Set<String> used = new HashSet<String>();
        for (String t : DEFAULT_ORDER) {
            if (byType.containsKey(t)) {
                subproblems.add(new SubProblem(t, t, byType.get(t)));
                used.add(t);
            }
        }
        for (Map.Entry<String, List<DecisionVariable>> e : byType.entrySet()) {
            if (!used.contains(e.getKey())) {
                subproblems.add(new SubProblem(e.getKey(), e.getKey(), e.getValue()));
            }
        }
        return this;
    }

    public DecompositionWorkflow setDependency(String fromProblem, String toProblem) {
        if (!dependencyGraph.containsKey(toProblem)) {
            dependencyGraph.put(toProblem, new LinkedHashSet<String>());
        }
        dependencyGraph.get(toProblem).add(fromProblem);
        return this;
    }

    public DecompositionWorkflow addSubProblem(String name, List<DecisionVariable> variables,
                                               List<String> after) {
        String varType = variables.isEmpty() ? "custom" : variables.get(0).getVariableType();
        subproblems.add(new SubProblem(name, varType, variables));
        if (after != null) {
            for (String dep : after) {
                setDependency(dep, name);
            }
        }
        return this;
    }

    /** Topological order respecting dependencies; insertion order when acyclic-free. */
    private List<SubProblem> getExecutionOrder() {
        if (dependencyGraph.isEmpty()) {
            return new ArrayList<SubProblem>(subproblems);
        }
        Map<String, SubProblem> byName = new LinkedHashMap<String, SubProblem>();
        for (SubProblem sp : subproblems) {
            byName.put(sp.name, sp);
        }
        // Kahn's algorithm
        Map<String, Integer> indeg = new LinkedHashMap<String, Integer>();
        Map<String, List<String>> adj = new LinkedHashMap<String, List<String>>();
        for (SubProblem sp : subproblems) {
            indeg.put(sp.name, 0);
            adj.put(sp.name, new ArrayList<String>());
        }
        for (Map.Entry<String, Set<String>> e : dependencyGraph.entrySet()) {
            String to = e.getKey();
            for (String from : e.getValue()) {
                if (adj.containsKey(from) && indeg.containsKey(to)) {
                    adj.get(from).add(to);
                    indeg.put(to, indeg.get(to) + 1);
                }
            }
        }
        List<String> queue = new ArrayList<String>();
        for (SubProblem sp : subproblems) {
            if (indeg.get(sp.name) == 0) {
                queue.add(sp.name);
            }
        }
        List<SubProblem> order = new ArrayList<SubProblem>();
        int head = 0;
        while (head < queue.size()) {
            String n = queue.get(head++);
            if (byName.containsKey(n)) {
                order.add(byName.get(n));
            }
            for (String m : adj.get(n)) {
                indeg.put(m, indeg.get(m) - 1);
                if (indeg.get(m) == 0) {
                    queue.add(m);
                }
            }
        }
        if (order.size() != subproblems.size()) {
            return new ArrayList<SubProblem>(subproblems);   // cycle: fall back
        }
        return order;
    }

    public WorkflowResult solveSequential(int maxCycles, double tolerance) {
        long start = System.nanoTime();
        WorkflowResult result = new WorkflowResult();
        if (subproblems.isEmpty()) {
            result.converged = true;
            return result;
        }
        List<SubProblem> ordered = getExecutionOrder();
        Map<String, Object> fixedValues = new LinkedHashMap<String, Object>();
        double prevObjective = Double.POSITIVE_INFINITY;

        for (int cycle = 0; cycle < maxCycles; cycle++) {
            for (SubProblem sp : ordered) {
                OptimizationProblem partial = createPartialProblem(sp, fixedValues);
                OptimizationResult spResult = new LineOptSolver(partial, solverOptions).solve();
                SubProblemResult spr = new SubProblemResult(sp.name, spResult);
                spr.variablesFixed.putAll(fixedValues);
                result.subproblemResults.put(sp.name, spr);
                for (Map.Entry<String, Object> e : spResult.variableValues.entrySet()) {
                    fixedValues.put(e.getKey(), e.getValue());
                }
            }
            double currentObjective = evaluateFullObjective(fixedValues);
            result.objectiveHistory.add(currentObjective);
            if (Math.abs(currentObjective - prevObjective) < tolerance) {
                result.converged = true;
                break;
            }
            prevObjective = currentObjective;
            result.cyclesCompleted = cycle + 1;
        }

        result.finalObjective = result.objectiveHistory.isEmpty()
                ? Double.POSITIVE_INFINITY
                : result.objectiveHistory.get(result.objectiveHistory.size() - 1);
        result.finalVariableValues.putAll(fixedValues);
        result.totalSolveTime = (System.nanoTime() - start) / 1e9;
        return result;
    }

    public WorkflowResult solveHierarchical() {
        return solveSequential(1, 0.0);
    }

    // ---- LayeredNetwork (LQN) layer-wise decomposition ------------------

    /** LQN layer-wise decomposition with default settings. */
    public WorkflowResult solveLayered() {
        return solveLayered(10, 0.01, true, 1e-3, null);
    }

    /**
     * Solve an LQN by layer, optionally freezing converged layers. Groups
     * decision variables by the LQN layer they perturb (host or task layer) and
     * cycles Gauss-Seidel over the layer groups, fixing every other layer's
     * variables at their current values while one layer is optimized. This is the
     * LQN analogue of {@link #solveSequential}, but the subproblems are LAYERS.
     *
     * <p>Freezing has two composable sources: {@code frozenLayersSeed}, an
     * explicit set held fixed throughout; and {@code autoFreeze}, adaptive --
     * after each cycle a layer whose representative node metrics moved less than
     * {@code freezeTol} (relative) is frozen and skipped, and unfrozen again if
     * any still-active layer later moves by more than {@code freezeTol}.</p>
     *
     * <p>Convergence is on the full penalized objective delta between cycles, or
     * when every layer is frozen. Falls back to {@link #solveSequential} for a
     * flat network. The returned {@link WorkflowResult} carries the extra fields
     * {@code frozenLayers} (final frozen set) and {@code modelEvaluations} (total
     * LINE solves). Mirrors native-Python {@code solveLayered}.</p>
     *
     * @param maxCycles        maximum Gauss-Seidel cycles
     * @param tolerance        objective-delta convergence tolerance
     * @param autoFreeze       enable adaptive layer freezing
     * @param freezeTol        relative-movement threshold for freezing a layer
     * @param frozenLayersSeed explicit layers held fixed throughout (may be null)
     * @return the workflow result with LQN freezing diagnostics
     */
    public WorkflowResult solveLayered(int maxCycles, double tolerance, boolean autoFreeze,
                                       double freezeTol, List<String> frozenLayersSeed) {
        long start = System.nanoTime();
        if (!problem.isLayered()) {
            return solveSequential(maxCycles, tolerance);
        }
        LayeredNetwork model = problem.getLqnModel();

        // Group variables by their primary (first) layer.
        Map<String, List<DecisionVariable>> groups = new LinkedHashMap<String, List<DecisionVariable>>();
        for (DecisionVariable var : problem.getVariables()) {
            List<String> layers = var.getLayer(model);
            String key = (layers == null || layers.isEmpty()) ? "_nolayer" : layers.get(0);
            if (!groups.containsKey(key)) {
                groups.put(key, new ArrayList<DecisionVariable>());
            }
            groups.get(key).add(var);
        }

        Objective objective = problem.getObjective();
        double penaltyWeight = solverOptions.penaltyWeight;
        LineEvaluator evaluator = new LineEvaluator(model, problem.getVariables(), null);

        Set<String> frozen = new LinkedHashSet<String>();
        if (frozenLayersSeed != null) {
            frozen.addAll(frozenLayersSeed);
        }
        Map<String, Object> fixedValues = new LinkedHashMap<String, Object>();
        Map<String, double[]> prevSig = new LinkedHashMap<String, double[]>();
        double prevObjective = Double.POSITIVE_INFINITY;
        int modelEvaluations = 0;

        WorkflowResult result = new WorkflowResult();
        for (int cycle = 0; cycle < maxCycles; cycle++) {
            for (Map.Entry<String, List<DecisionVariable>> ge : groups.entrySet()) {
                String layer = ge.getKey();
                if (frozen.contains(layer)) {
                    continue;
                }
                List<DecisionVariable> layerVars = ge.getValue();
                OptimizationProblem partial = new OptimizationProblem(model);
                for (DecisionVariable var : layerVars) {
                    partial.addVariable(var);
                }
                // Fix every other layer's variables at their current values.
                Set<String> subNames = new HashSet<String>();
                for (DecisionVariable v : layerVars) {
                    subNames.add(v.getName());
                }
                List<Pair<DecisionVariable, Object>> fixedPairs =
                        new ArrayList<Pair<DecisionVariable, Object>>();
                for (DecisionVariable var : problem.getVariables()) {
                    String name = var.getName();
                    if (!subNames.contains(name) && fixedValues.containsKey(name)) {
                        fixedPairs.add(new Pair<DecisionVariable, Object>(var, fixedValues.get(name)));
                    }
                }
                partial.setFixedVariables(fixedPairs);
                partial.setObjective(objective);
                for (Constraint c : problem.getConstraints()) {
                    partial.addConstraint(c);
                }
                OptimizationResult spResult = new LineOptSolver(partial, solverOptions).solve();
                modelEvaluations += spResult.modelEvaluations;

                SubProblemResult spr = new SubProblemResult(layer, spResult);
                spr.variablesFixed.putAll(fixedValues);
                result.subproblemResults.put(layer, spr);
                for (Map.Entry<String, Object> e : spResult.variableValues.entrySet()) {
                    fixedValues.put(e.getKey(), e.getValue());
                }
            }

            // Full evaluation: objective + per-layer signatures for freezing.
            EvaluationResult evalResult = evaluator.evaluateValues(fixedValues);
            modelEvaluations += 1;
            double currentObjective;
            Map<String, double[]> sig;
            if (!evalResult.feasible) {
                currentObjective = Double.POSITIVE_INFINITY;
                sig = new LinkedHashMap<String, double[]>();
            } else {
                currentObjective = objective.evaluateWithPenalty(evalResult, fixedValues, penaltyWeight);
                for (Constraint c : problem.getConstraints()) {
                    currentObjective += c.evaluate(evalResult, fixedValues) * penaltyWeight;
                }
                sig = layerSignatures(evalResult, groups.keySet());
            }

            if (autoFreeze && !prevSig.isEmpty()) {
                Map<String, Double> moved = new LinkedHashMap<String, Double>();
                for (String L : groups.keySet()) {
                    moved.put(L, sigDelta(prevSig.get(L), sig.get(L)));
                }
                boolean activeMoved = false;
                for (String a : groups.keySet()) {
                    if (!frozen.contains(a) && moved.get(a) > freezeTol) {
                        activeMoved = true;
                        break;
                    }
                }
                for (String L : new ArrayList<String>(groups.keySet())) {
                    if (frozen.contains(L)) {
                        // Unfreeze if a still-active layer has moved.
                        if (activeMoved) {
                            frozen.remove(L);
                        }
                    } else if (moved.get(L) < freezeTol) {
                        frozen.add(L);
                    }
                }
            }

            prevSig = sig;
            result.objectiveHistory.add(currentObjective);
            if (Math.abs(currentObjective - prevObjective) < tolerance) {
                result.converged = true;
                break;
            }
            prevObjective = currentObjective;
            result.cyclesCompleted = cycle + 1;
            // All layers frozen: nothing left to optimize.
            if (frozen.size() >= groups.size()) {
                result.converged = true;
                break;
            }
        }

        result.finalObjective = result.objectiveHistory.isEmpty()
                ? Double.POSITIVE_INFINITY
                : result.objectiveHistory.get(result.objectiveHistory.size() - 1);
        result.finalVariableValues.putAll(fixedValues);
        result.totalSolveTime = (System.nanoTime() - start) / 1e9;
        List<String> frozenSorted = new ArrayList<String>(frozen);
        Collections.sort(frozenSorted);
        result.frozenLayers = frozenSorted;
        result.modelEvaluations = modelEvaluations;
        return result;
    }

    /**
     * Representative (Util, QLen, Tput, RespT) per layer, keyed by its node. A
     * layer named after a processor or task has a same-named node in the LQN
     * average table; its metrics are the layer's convergence signature. Layers
     * without a matching node get a null signature so they never auto-freeze.
     */
    private static Map<String, double[]> layerSignatures(EvaluationResult evalResult, Set<String> layers) {
        Map<String, double[]> sig = new LinkedHashMap<String, double[]>();
        for (String layer : layers) {
            double util = evalResult.getUtilization(layer);
            double qlen = evalResult.getQueueLength(layer);
            double tput = evalResult.getThroughput(layer);
            double respt = evalResult.getResponseTime(layer);
            boolean informative = false;
            double[] triad = {util, qlen, tput};
            for (double v : triad) {
                if (v != 0.0 && !Double.isInfinite(v)) {
                    informative = true;
                    break;
                }
            }
            sig.put(layer, informative ? new double[]{util, qlen, tput, respt} : null);
        }
        return sig;
    }

    /** Max relative change between two layer signatures (inf if unknown). */
    private static double sigDelta(double[] a, double[] b) {
        if (a == null || b == null) {
            return Double.POSITIVE_INFINITY;
        }
        double eps = 1e-12;
        double delta = 0.0;
        int n = Math.min(a.length, b.length);
        for (int i = 0; i < n; i++) {
            double ai = a[i];
            double bi = b[i];
            if (Double.isNaN(ai) || Double.isInfinite(ai) || Double.isNaN(bi) || Double.isInfinite(bi)) {
                continue;
            }
            delta = Math.max(delta, Math.abs(bi - ai) / (Math.abs(ai) + eps));
        }
        return delta;
    }

    /** New empty problem over the same model (flat or LQN). */
    private OptimizationProblem newPartialProblem() {
        return problem.isLayered()
                ? new OptimizationProblem(problem.getLqnModel())
                : new OptimizationProblem(problem.getModel());
    }

    private OptimizationProblem createPartialProblem(SubProblem subproblem,
                                                     Map<String, Object> fixedValues) {
        OptimizationProblem partial = newPartialProblem();
        for (DecisionVariable var : subproblem.variables) {
            partial.addVariable(var);
        }
        Set<String> subNames = new HashSet<String>(subproblem.getVariableNames());
        List<Pair<DecisionVariable, Object>> fixedPairs =
                new ArrayList<Pair<DecisionVariable, Object>>();
        for (DecisionVariable var : problem.getVariables()) {
            String name = var.getName();
            if (!subNames.contains(name) && fixedValues.containsKey(name)) {
                fixedPairs.add(new Pair<DecisionVariable, Object>(var, fixedValues.get(name)));
            }
        }
        partial.setFixedVariables(fixedPairs);
        partial.setObjective(problem.getObjective());
        for (Constraint c : problem.getConstraints()) {
            partial.addConstraint(c);
        }
        for (Pair<Network, Double> s : problem.getScenarios()) {
            partial.addScenario(s.getLeft(), s.getRight());
        }
        subproblem.fixedValues = new LinkedHashMap<String, Object>(fixedValues);
        return partial;
    }

    private double evaluateFullObjective(Map<String, Object> variableValues) {
        double penaltyWeight = solverOptions.penaltyWeight;
        LineEvaluator evaluator = problem.isLayered()
                ? new LineEvaluator(problem.getLqnModel(), problem.getVariables(), null)
                : new LineEvaluator(problem.getModel(), problem.getVariables(), null);
        EvaluationResult res = evaluator.evaluateValues(variableValues);
        if (!res.feasible) {
            return Double.POSITIVE_INFINITY;
        }
        Objective objective = problem.getObjective();
        double value = objective.evaluateWithPenalty(res, variableValues, penaltyWeight);
        for (Constraint c : problem.getConstraints()) {
            value += c.evaluate(res, variableValues) * penaltyWeight;
        }
        return value;
    }
}
