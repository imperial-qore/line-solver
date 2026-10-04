package jline.opt.solver;

import jline.lang.layered.LayeredNetwork;
import jline.opt.LineEvaluator;
import jline.opt.OptimizationProblem;
import jline.opt.objectives.Constraint;
import jline.opt.objectives.Objective;
import jline.opt.results.EvaluationResult;
import jline.opt.results.OptimizationResult;
import jline.opt.results.SensitivityData;
import jline.opt.solver.de.DifferentialEvolution;
import jline.opt.variables.DecisionVariable;
import jline.opt.variables.HostDemand;
import jline.opt.variables.ServiceRate;
import jline.util.Pair;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Main line-opt solver. Minimizes a single penalized scalar objective (with
 * constraints folded in as penalties) over the decision variables, aggregating
 * across workload scenarios. Uses the self-contained {@link DifferentialEvolution}
 * engine (numpy-exact RNG) for the evolutionary path, or an analytic/FD-driven
 * projected-gradient path for continuous problems. Mirrors native-Python
 * {@code line_solver.opt.solver.LineOptSolver}.
 */
public class LineOptSolver {

    // Continuous (differentiable) types, including the continuous LQN knobs.
    private static final Set<String> CONTINUOUS_TYPES = new HashSet<String>();

    static {
        CONTINUOUS_TYPES.add("service_rate");
        CONTINUOUS_TYPES.add("routing");
        CONTINUOUS_TYPES.add("host_demand");
        CONTINUOUS_TYPES.add("think_time");
    }

    private static class TimeLimitReached extends RuntimeException {
    }

    private final OptimizationProblem problem;
    private final LineOptSolverOptions opt;
    private final Map<String, Object> fixedValueMap = new LinkedHashMap<String, Object>();
    /** Free (post-freeze) decision variables the optimizer searches over. */
    private final List<DecisionVariable> freeVariables;
    private final List<LineEvaluator> evaluators = new ArrayList<LineEvaluator>();
    private final double[] scenarioWeights;
    private List<Map<String, EvaluationResult>> caches;

    private int iterations = 0;
    private double bestValue = Double.POSITIVE_INFINITY;
    private double[] bestX = null;
    private final List<Double> convergenceHistory = new ArrayList<Double>();
    private long startTime = 0;
    private long deadlineNanos = Long.MAX_VALUE;

    // LQN gradient bookkeeping: a per-configuration sensitivity cache and a
    // gradient-evaluation counter driving the partial_plus_fd refresh.
    private final Map<String, Map<String, Map<String, Double>>> lqnSensCache =
            new LinkedHashMap<String, Map<String, Map<String, Double>>>();
    private int gradCalls = 0;

    public LineOptSolver(OptimizationProblem problem, LineOptSolverOptions options) {
        this.problem = problem;
        this.opt = options;

        List<Pair<DecisionVariable, Object>> fixed =
                new ArrayList<Pair<DecisionVariable, Object>>(problem.getFixedVariables());
        for (Pair<DecisionVariable, Object> fv : fixed) {
            fixedValueMap.put(fv.getLeft().getName(), fv.getRight());
        }

        // Explicit layer freezing (LQN only): variables whose layer is in
        // frozen_layers move from free to fixed, held at the model's current
        // value.
        this.freeVariables = freezeLayers(problem.getVariables(), fixed);

        if (problem.isLayered()) {
            evaluators.add(new LineEvaluator(problem.getLqnModel(), freeVariables, fixed));
            // Workload scenarios are supported for flat models only.
            this.scenarioWeights = new double[]{1.0};
        } else {
            evaluators.add(new LineEvaluator(problem.getModel(), freeVariables, fixed));
            List<Pair<jline.lang.Network, Double>> scenarios = problem.getScenarios();
            double[] weights = new double[1 + scenarios.size()];
            weights[0] = 1.0;
            for (int i = 0; i < scenarios.size(); i++) {
                evaluators.add(new LineEvaluator(scenarios.get(i).getLeft(),
                        freeVariables, fixed));
                weights[i + 1] = scenarios.get(i).getRight();
            }
            this.scenarioWeights = weights;
        }
    }

    /**
     * Partition variables by the {@code frozenLayers} option (LQN only). A
     * variable whose layer set ({@link DecisionVariable#getLayer}) intersects the
     * frozen set is moved from the free list into {@code fixed}, held at its
     * current model parameter value ({@link DecisionVariable#currentValue}). If
     * the current value cannot be read the variable is dropped from the
     * optimization, leaving the model's built-in value. A no-op when the model is
     * flat or no layers are frozen. Mirrors native-Python {@code _freezeLayers}.
     */
    private List<DecisionVariable> freezeLayers(List<DecisionVariable> freeVars,
                                                List<Pair<DecisionVariable, Object>> fixed) {
        List<String> frozen = opt.frozenLayers;
        if (frozen == null || frozen.isEmpty() || !problem.isLayered()) {
            return freeVars;
        }
        LayeredNetwork model = problem.getLqnModel();
        Set<String> frozenSet = new HashSet<String>(frozen);
        Set<String> already = new HashSet<String>();
        for (Pair<DecisionVariable, Object> fv : fixed) {
            already.add(fv.getLeft().getName());
        }
        List<DecisionVariable> kept = new ArrayList<DecisionVariable>();
        for (DecisionVariable var : freeVars) {
            List<String> layers = var.getLayer(model);
            if (layers == null) {
                layers = Collections.emptyList();
            }
            if (intersects(frozenSet, layers)) {
                if (already.contains(var.getName())) {
                    continue;
                }
                Object value = var.currentValue(model);
                if (value != null) {
                    fixed.add(new Pair<DecisionVariable, Object>(var, value));
                    fixedValueMap.put(var.getName(), value);
                }
                // else: drop the variable; the model keeps its built-in value
            } else {
                kept.add(var);
            }
        }
        return kept;
    }

    private static boolean intersects(Set<String> a, List<String> b) {
        for (String s : b) {
            if (a.contains(s)) {
                return true;
            }
        }
        return false;
    }

    public OptimizationResult solve() {
        startTime = System.nanoTime();
        iterations = 0;
        bestValue = Double.POSITIVE_INFINITY;
        bestX = null;
        convergenceHistory.clear();
        caches = new ArrayList<Map<String, EvaluationResult>>();
        for (int i = 0; i < evaluators.size(); i++) {
            caches.add(new LinkedHashMap<String, EvaluationResult>());
        }
        lqnSensCache.clear();
        gradCalls = 0;
        deadlineNanos = startTime + (long) (opt.timeLimit * 1e9);

        double[][] bounds = evaluators.get(0).getBounds();
        if (bounds.length == 0) {
            return buildEmptyResult();
        }

        if (shouldUseGradient()) {
            return solveGradient(bounds);
        }
        return solveEvolution(bounds);
    }

    // ---- evolutionary path -------------------------------------------------

    private OptimizationResult solveEvolution(double[][] bounds) {
        int dim = bounds.length;
        double[] low = new double[dim];
        double[] high = new double[dim];
        for (int i = 0; i < dim; i++) {
            low[i] = bounds[i][0];
            high[i] = bounds[i][1];
        }
        long seed = (opt.seed != null) ? opt.seed : System.nanoTime();

        DifferentialEvolution de = new DifferentialEvolution(
                new DifferentialEvolution.Objective() {
                    public double evaluate(double[] x) {
                        return objectiveFunction(x);
                    }
                },
                low, high, opt.strategy, opt.popsize, opt.maxIterations,
                opt.mutationLow, opt.mutationHigh, opt.recombination, opt.tol, seed);

        de.setCallback(new DifferentialEvolution.Callback() {
            public boolean call(double[] bestXk, int nit) {
                iterations = nit;
                convergenceHistory.add(bestValue);
                return System.nanoTime() >= deadlineNanos;
            }
        });

        boolean timedOut = false;
        double[] resX;
        double resFun;
        try {
            DifferentialEvolution.Result r = de.solve();
            resX = null;
            resFun = r.fun;
            // recover the decision vector from the engine's best
            resX = engineBestVector(r);
        } catch (TimeLimitReached e) {
            timedOut = true;
            if (bestX == null) {
                return buildEmptyResult();
            }
            resX = bestX;
            resFun = bestValue;
        }

        double solveTime = (System.nanoTime() - startTime) / 1e9;
        OptimizationResult result = buildResult(resX, resFun, solveTime);
        if (timedOut) {
            result.terminatedBy = "time_limit";
        }
        return result;
    }

    private double[] engineBestVector(DifferentialEvolution.Result r) {
        // r.x is the scaled best parameter vector; with [0,1] bounds this is the
        // decision vector directly. Prefer the tracked best-so-far when lower.
        if (bestX != null && bestValue <= r.fun) {
            return bestX;
        }
        return r.x;
    }

    private double objectiveFunction(double[] x) {
        if (System.nanoTime() >= deadlineNanos) {
            throw new TimeLimitReached();
        }
        Map<String, Object> values = evaluators.get(0).decodeVariables(x);
        Map<String, Object> allValues = new LinkedHashMap<String, Object>(fixedValueMap);
        allValues.putAll(values);

        Objective objective = problem.getObjective();
        double pw = opt.penaltyWeight;
        List<Double> scenarioValues = new ArrayList<Double>();
        for (int i = 0; i < evaluators.size(); i++) {
            EvaluationResult res = evaluators.get(i).evaluateValuesWithCache(values, caches.get(i));
            if (!res.feasible) {
                return Double.POSITIVE_INFINITY;
            }
            double value = objective.evaluateWithPenalty(res, allValues, pw);
            for (Constraint c : problem.getConstraints()) {
                value += c.evaluate(res, allValues) * pw;
            }
            scenarioValues.add(value);
        }
        double total = aggregateScenarios(scenarioValues);
        if (total < bestValue) {
            bestValue = total;
            bestX = x.clone();
        }
        return total;
    }

    private double aggregateScenarios(List<Double> values) {
        if (values.size() == 1) {
            return values.get(0);
        }
        if ("mean".equals(opt.scenarioAggregation)) {
            double totalWeight = 0.0;
            double sum = 0.0;
            for (int i = 0; i < values.size(); i++) {
                sum += scenarioWeights[i] * values.get(i);
                totalWeight += scenarioWeights[i];
            }
            return sum / totalWeight;
        }
        double max = Double.NEGATIVE_INFINITY;
        for (double v : values) {
            if (v > max) {
                max = v;
            }
        }
        return max;
    }

    // ---- gradient path -----------------------------------------------------

    private boolean shouldUseGradient() {
        if ("gradient".equals(opt.optimizer)) {
            return true;
        }
        if ("evolution".equals(opt.optimizer)) {
            return false;
        }
        return allContinuous();   // 'auto'
    }

    private boolean allContinuous() {
        List<DecisionVariable> vars = freeVariables;
        if (vars.isEmpty()) {
            return false;
        }
        for (DecisionVariable v : vars) {
            if (!CONTINUOUS_TYPES.contains(v.getVariableType())) {
                return false;
            }
        }
        return true;
    }

    private OptimizationResult solveGradient(double[][] bounds) {
        int dim = bounds.length;
        int nStart = Math.max(1, opt.gradientRestarts);
        java.util.Random rng = new java.util.Random(opt.seed != null ? opt.seed : 0L);
        List<double[]> starts = new ArrayList<double[]>();
        double[] mid = new double[dim];
        for (int i = 0; i < dim; i++) {
            mid[i] = 0.5;
        }
        starts.add(mid);
        for (int s = 1; s < nStart; s++) {
            double[] x0 = new double[dim];
            for (int i = 0; i < dim; i++) {
                x0[i] = rng.nextDouble();
            }
            starts.add(x0);
        }

        double[] best = null;
        double bestFun = Double.POSITIVE_INFINITY;
        boolean timedOut = false;
        for (double[] x0 : starts) {
            if (System.nanoTime() >= deadlineNanos) {
                timedOut = true;
                break;
            }
            try {
                double[] xx = projectedGradientDescent(x0);
                double fx = objectiveFunction(xx);
                if (isFinite(fx) && fx < bestFun) {
                    bestFun = fx;
                    best = xx.clone();
                }
            } catch (TimeLimitReached e) {
                timedOut = true;
                break;
            }
        }

        if (best == null || (bestX != null && bestValue < bestFun)) {
            if (bestX != null) {
                best = bestX;
                bestFun = bestValue;
            } else if (best == null) {
                return buildEmptyResult();
            }
        }

        double solveTime = (System.nanoTime() - startTime) / 1e9;
        OptimizationResult result = buildResult(best, bestFun, solveTime);
        if (timedOut) {
            result.terminatedBy = "time_limit";
        }
        return result;
    }

    /** Projected gradient descent with backtracking line search on [0,1]^dim. */
    private double[] projectedGradientDescent(double[] x0) {
        int dim = x0.length;
        double[] x = clip(x0.clone());
        double f = objectiveFunction(x);
        for (int iter = 0; iter < opt.maxIterations; iter++) {
            double[] g = objectiveGradient(x);
            double gnorm = 0.0;
            for (double gi : g) {
                gnorm += gi * gi;
            }
            if (Math.sqrt(gnorm) < 1e-9) {
                break;
            }
            double step = 1.0;
            boolean improved = false;
            for (int ls = 0; ls < 30; ls++) {
                double[] xn = new double[dim];
                for (int i = 0; i < dim; i++) {
                    xn[i] = x[i] - step * g[i];
                }
                xn = clip(xn);
                double fn = objectiveFunction(xn);
                if (isFinite(fn) && fn < f - 1e-12) {
                    x = xn;
                    f = fn;
                    improved = true;
                    break;
                }
                step *= 0.5;
            }
            if (!improved) {
                break;
            }
        }
        return x;
    }

    private static double[] clip(double[] x) {
        for (int i = 0; i < x.length; i++) {
            if (x[i] < 0.0) {
                x[i] = 0.0;
            }
            if (x[i] > 1.0) {
                x[i] = 1.0;
            }
        }
        return x;
    }

    /**
     * Gradient of the penalized objective in encoded space. For a LayeredNetwork
     * the source is selected by the {@code lqnGradient} option: 'fd' (whole-model
     * finite difference, robust default), 'partial_sens' (SolverLN per-layer
     * partial derivatives, cheap/biased), or 'partial_plus_fd' (partial with a
     * periodic full-FD correction). For a flat network it prefers the analytic
     * product-form gradient. All paths fall back to central finite differences.
     */
    private double[] objectiveGradient(double[] x) {
        if (evaluators.get(0).isLayered()) {
            String mode = opt.lqnGradient;
            if ("partial_sens".equals(mode) || "partial_plus_fd".equals(mode)) {
                gradCalls++;
                int refresh = Math.max(1, opt.fdRefresh);
                // partial_plus_fd periodically replaces the biased partial
                // direction with the correct whole-model finite difference.
                if (!("partial_plus_fd".equals(mode) && gradCalls % refresh == 0)) {
                    double[] g = lqnAnalyticGradient(x);
                    if (g != null) {
                        return g;
                    }
                }
            }
            return finiteDifferenceGradient(x);
        }
        double[] analytic = analyticGradient(x);
        if (analytic != null) {
            return analytic;
        }
        return finiteDifferenceGradient(x);
    }

    /**
     * Central finite-difference gradient of the penalized scalar objective.
     * Works for any model (for a LayeredNetwork each perturbed evaluation
     * re-solves the whole ensemble, giving the correct total derivative).
     * One-sided differences are used near an infeasible boundary.
     */
    private double[] finiteDifferenceGradient(double[] x) {
        // A layered evaluation re-solves an iterative fixed point, so the step
        // must clear its noise floor (see fdStepLayered).
        double h = evaluators.get(0).isLayered() ? opt.fdStepLayered : opt.fdStep;
        int dim = x.length;
        double[] g = new double[dim];
        Double f0 = null;
        for (int i = 0; i < dim; i++) {
            double[] xp = x.clone();
            double[] xm = x.clone();
            xp[i] = Math.min(1.0, x[i] + h);
            xm[i] = Math.max(0.0, x[i] - h);
            double fp = objectiveFunction(xp);
            double fm = objectiveFunction(xm);
            if (isFinite(fp) && isFinite(fm) && xp[i] > xm[i]) {
                g[i] = (fp - fm) / (xp[i] - xm[i]);
                continue;
            }
            if (f0 == null) {
                f0 = objectiveFunction(x);
            }
            if (isFinite(fp) && isFinite(f0) && xp[i] > x[i]) {
                g[i] = (fp - f0) / (xp[i] - x[i]);
            } else if (isFinite(fm) && isFinite(f0) && x[i] > xm[i]) {
                g[i] = (f0 - fm) / (x[i] - xm[i]);
            } else {
                g[i] = 0.0;
            }
        }
        return g;
    }

    /**
     * Partial-sensitivity gradient for a LayeredNetwork, or {@code null}.
     * Assembles d(objective)/dx from SolverLN's per-layer WITHIN-LAYER
     * service-rate derivatives, WITHOUT re-solving per parameter: for each
     * host-demand variable it reads d(layer metric)/d(rate) from the sensitivity
     * table row, maps each layer metric to the LQN node metric it approximates
     * and takes d(penalized scalar)/d(node metric) by cheap metric-space finite
     * differences, then chains through d(rate)/d(demand) = -1/D^2 and the linear
     * decode. Returns {@code null} (caller finite-differences the whole model)
     * when multi-scenario, any variable is not a host-demand variable, or the
     * sensitivity table is unavailable. BIASED (omits cross-layer coupling);
     * 'fd'/'partial_plus_fd' correct that. Mirrors native-Python
     * {@code _lqnAnalyticGradient}.
     */
    private double[] lqnAnalyticGradient(double[] x) {
        List<DecisionVariable> vars = freeVariables;
        if (evaluators.size() != 1 || vars.isEmpty()) {
            return null;
        }
        for (DecisionVariable v : vars) {
            if (!"host_demand".equals(v.getVariableType()) || !(v instanceof HostDemand)) {
                return null;
            }
        }
        Map<String, Object> values = evaluators.get(0).decodeVariables(x);
        Map<String, Object> allValues = new LinkedHashMap<String, Object>(fixedValueMap);
        allValues.putAll(values);
        EvaluationResult res = evaluators.get(0).evaluateValuesWithCache(values, caches.get(0));
        if (!res.feasible) {
            return null;
        }

        // Per-configuration sensitivity cache (the table is expensive).
        String skey = LineEvaluator.valuesKey(values);
        Map<String, Map<String, Double>> sens;
        if (lqnSensCache.containsKey(skey)) {
            sens = lqnSensCache.get(skey);
        } else {
            sens = evaluators.get(0).evaluateLayeredSensitivities(values);
            lqnSensCache.put(skey, sens);
        }
        if (sens == null || sens.isEmpty()) {
            return null;
        }

        double pw = opt.penaltyWeight;
        double h = opt.fdStep;
        LayeredNetwork model = problem.getLqnModel();

        double[] grad = new double[x.length];
        for (int i = 0; i < vars.size(); i++) {
            HostDemand hd = (HostDemand) vars.get(i);
            String sk = hd.sensKey(model);
            if (sk == null || !sens.containsKey(sk)) {
                // No sensitivity row for this variable: leave 0.
                continue;
            }
            Map<String, Double> row = sens.get(sk);
            Map<String, String> targets = hd.sensMetricTargets(model);
            double dSdRate = 0.0;
            for (Map.Entry<String, Double> e : row.entrySet()) {
                double dMetricDRate = e.getValue();
                if (dMetricDRate == 0.0) {
                    continue;
                }
                double dSdMetric = scalarMetricDerivative(res, allValues, e.getKey(),
                        targets.get(e.getKey()), pw, h);
                dSdRate += dSdMetric * dMetricDRate;
            }
            Object demandObj = allValues.get(hd.getName());
            double rateJac = (demandObj instanceof Number)
                    ? hd.rateJacobian(((Number) demandObj).doubleValue()) : 0.0;
            grad[i] = dSdRate * rateJac * hd.decodeJacobian(new double[]{x[i]});
        }
        return grad;
    }

    /**
     * Analytic gradient of the penalized objective for a single-scenario,
     * all-rate-variable problem, assembled from product-form sensitivities and
     * cheap metric/value-space finite differences (no extra solves). Returns
     * {@code null} otherwise. Mirrors native-Python {@code _analyticGradient}.
     */
    private double[] analyticGradient(double[] x) {
        List<DecisionVariable> vars = freeVariables;
        if (evaluators.size() != 1 || vars.isEmpty()) {
            return null;
        }
        for (DecisionVariable v : vars) {
            if (!"service_rate".equals(v.getVariableType()) || !(v instanceof ServiceRate)) {
                return null;
            }
        }
        Map<String, Object> values = evaluators.get(0).decodeVariables(x);
        Map<String, Object> allValues = new LinkedHashMap<String, Object>(fixedValueMap);
        allValues.putAll(values);
        EvaluationResult res = evaluators.get(0).evaluateValuesWithCache(values, caches.get(0));
        SensitivityData sens = res.sensitivities;
        if (!res.feasible || sens == null) {
            return null;
        }

        double pw = opt.penaltyWeight;
        double h = opt.fdStep;

        // d(scalar)/d(rate param) accumulated over all metric kinds
        Map<String, Double> dSdRate = new LinkedHashMap<String, Double>();
        for (Map.Entry<String, Map<String, Map<String, Double>>> kindEntry : sens.data.entrySet()) {
            String kind = kindEntry.getKey();
            for (Map.Entry<String, Map<String, Double>> metricEntry : kindEntry.getValue().entrySet()) {
                double dSdm = scalarMetricDerivative(res, allValues, kind, metricEntry.getKey(), pw, h);
                if (dSdm == 0.0) {
                    continue;
                }
                for (Map.Entry<String, Double> pe : metricEntry.getValue().entrySet()) {
                    double prev = dSdRate.containsKey(pe.getKey()) ? dSdRate.get(pe.getKey()) : 0.0;
                    dSdRate.put(pe.getKey(), prev + dSdm * pe.getValue());
                }
            }
        }

        double[] grad = new double[x.length];
        for (int i = 0; i < vars.size(); i++) {
            ServiceRate var = (ServiceRate) vars.get(i);
            String name = var.getName();
            String pkey = var.paramKey();
            double direct = 0.0;
            Object baseVal = allValues.get(name);
            if (baseVal instanceof Number) {
                double bv = ((Number) baseVal).doubleValue();
                allValues.put(name, bv + h);
                double fp = scalar(res, allValues, pw);
                allValues.put(name, bv - h);
                double fm = scalar(res, allValues, pw);
                allValues.put(name, bv);
                direct = (fp - fm) / (2.0 * h);
            }
            double dSdvalue = direct + (dSdRate.containsKey(pkey) ? dSdRate.get(pkey) : 0.0);
            grad[i] = dSdvalue * var.decodeJacobian(new double[]{x[i]});
        }
        return grad;
    }

    private double scalar(EvaluationResult res, Map<String, Object> vals, double pw) {
        double v = problem.getObjective().evaluateWithPenalty(res, vals, pw);
        for (Constraint c : problem.getConstraints()) {
            v += c.evaluate(res, vals) * pw;
        }
        return v;
    }

    private double scalarMetricDerivative(EvaluationResult res, Map<String, Object> allValues,
                                          String kind, String mkey, double pw, double h) {
        // central FD of the penalized scalar w.r.t. a stored metric (no solving)
        Double base = readMetric(res, kind, mkey);
        if (base == null) {
            return 0.0;
        }
        writeMetric(res, kind, mkey, base + h);
        double fp = scalar(res, allValues, pw);
        writeMetric(res, kind, mkey, base - h);
        double fm = scalar(res, allValues, pw);
        writeMetric(res, kind, mkey, base);
        return (fp - fm) / (2.0 * h);
    }

    private Double readMetric(EvaluationResult res, String kind, String mkey) {
        if ("Util".equals(kind)) {
            return res.utilizations.get(mkey);
        }
        String[] sc = mkey.split("\\|\\|");
        if (sc.length != 2) {
            return null;
        }
        Map<String, Map<String, Double>> m = metricMap(res, kind);
        if (m == null || !m.containsKey(sc[0]) || !m.get(sc[0]).containsKey(sc[1])) {
            return null;
        }
        return m.get(sc[0]).get(sc[1]);
    }

    private void writeMetric(EvaluationResult res, String kind, String mkey, double value) {
        if ("Util".equals(kind)) {
            res.utilizations.put(mkey, value);
            return;
        }
        String[] sc = mkey.split("\\|\\|");
        if (sc.length != 2) {
            return;
        }
        Map<String, Map<String, Double>> m = metricMap(res, kind);
        if (m != null && m.containsKey(sc[0])) {
            m.get(sc[0]).put(sc[1], value);
        }
    }

    private Map<String, Map<String, Double>> metricMap(EvaluationResult res, String kind) {
        if ("RespT".equals(kind)) {
            return res.responseTimes;
        }
        if ("QLen".equals(kind)) {
            return res.queueLengths;
        }
        if ("Tput".equals(kind)) {
            return res.throughputs;
        }
        return null;
    }

    // ---- result assembly ---------------------------------------------------

    private OptimizationResult buildResult(double[] x, double objectiveValue, double solveTime) {
        OptimizationResult result = new OptimizationResult();
        result.objectiveValue = objectiveValue;
        result.variableValues.putAll(evaluators.get(0).decodeVariables(x));
        result.iterations = iterations;
        result.solveTime = solveTime;
        int evals = 0;
        for (LineEvaluator e : evaluators) {
            evals += e.getEvaluationCount();
        }
        result.modelEvaluations = evals;
        result.convergenceHistory.addAll(convergenceHistory);

        Map<String, Object> allValues = new LinkedHashMap<String, Object>(fixedValueMap);
        allValues.putAll(result.variableValues);
        Objective objective = problem.getObjective();
        List<Constraint> allConstraints = new ArrayList<Constraint>(objective.getConstraints());
        allConstraints.addAll(problem.getConstraints());

        result.feasible = true;
        for (int i = 0; i < evaluators.size(); i++) {
            EvaluationResult er = evaluators.get(i).evaluateValuesWithCache(
                    result.variableValues, caches.get(i));
            if (!er.feasible) {
                result.feasible = false;
                continue;
            }
            for (Constraint c : allConstraints) {
                double violation = c.evaluate(er, allValues);
                if (violation > 0) {
                    String name = c.getName();
                    double prev = result.constraintViolations.containsKey(name)
                            ? result.constraintViolations.get(name) : 0.0;
                    result.constraintViolations.put(name, Math.max(violation, prev));
                    result.feasible = false;
                }
            }
        }

        double elapsed = (System.nanoTime() - startTime) / 1e9;
        if (elapsed >= opt.timeLimit) {
            result.terminatedBy = "time_limit";
        } else if (iterations >= opt.maxIterations) {
            result.terminatedBy = "iterations";
        } else {
            result.terminatedBy = "convergence";
        }
        return result;
    }

    private OptimizationResult buildEmptyResult() {
        OptimizationResult result = new OptimizationResult();
        result.objectiveValue = 0.0;
        result.feasible = true;
        result.terminatedBy = "empty";
        return result;
    }

    private static boolean isFinite(double v) {
        return !Double.isNaN(v) && !Double.isInfinite(v);
    }
}
