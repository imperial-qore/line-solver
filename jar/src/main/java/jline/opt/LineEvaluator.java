package jline.opt;

import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.lang.layered.Entry;
import jline.lang.layered.LayeredNetwork;
import jline.lang.layered.Task;
import jline.opt.results.EvaluationResult;
import jline.opt.sensitivity.SensitivityTable;
import jline.opt.variables.DecisionVariable;
import jline.opt.variables.LqnDecisionVariable;
import jline.solvers.LayeredNetworkAvgTable;
import jline.solvers.NetworkAvgSysTable;
import jline.solvers.NetworkAvgTable;
import jline.solvers.auto.SolverAUTO;
import jline.solvers.ln.SolverLN;
import jline.util.Pair;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Interface between the optimizer and LINE's SolverAUTO. Applies decision
 * variable values to a per-evaluation model copy, solves, and extracts per-
 * (station, class) and system metrics into an {@link EvaluationResult}. Mirrors
 * native-Python {@code line_solver.opt.evaluator.LineEvaluator}.
 *
 * <p>Because a copied model carries a cached NetworkStruct that scalar setters
 * do not invalidate, {@link #evaluateValues} forces {@code refreshStruct} after
 * applying variables so the solve reflects the mutated configuration.</p>
 */
public class LineEvaluator {

    /** Non-null iff flat; null when layered. */
    private final Network baseModel;
    /** Non-null iff layered; null when flat. */
    private final LayeredNetwork lqnBase;
    private final boolean layered;
    private final List<DecisionVariable> variables;
    private final List<Pair<DecisionVariable, Object>> fixedVariables;
    private int evaluationCount = 0;
    private final int totalDimension;
    private final int[] varOffsets;

    public LineEvaluator(Network model, List<DecisionVariable> variables,
                         List<Pair<DecisionVariable, Object>> fixedVariables) {
        this(model, null, false, variables, fixedVariables);
    }

    public LineEvaluator(LayeredNetwork lqnModel, List<DecisionVariable> variables,
                         List<Pair<DecisionVariable, Object>> fixedVariables) {
        this(null, lqnModel, true, variables, fixedVariables);
    }

    private LineEvaluator(Network model, LayeredNetwork lqnModel, boolean layered,
                          List<DecisionVariable> variables,
                          List<Pair<DecisionVariable, Object>> fixedVariables) {
        this.baseModel = model;
        this.lqnBase = lqnModel;
        this.layered = layered;
        this.variables = variables;
        this.fixedVariables = (fixedVariables != null)
                ? fixedVariables : new ArrayList<Pair<DecisionVariable, Object>>();
        int dim = 0;
        this.varOffsets = new int[variables.size()];
        for (int i = 0; i < variables.size(); i++) {
            varOffsets[i] = dim;
            dim += variables.get(i).getDimension();
        }
        this.totalDimension = dim;
    }

    public Network getModel() {
        return baseModel;
    }

    /** True if this evaluator drives a LayeredNetwork (LQN) model. */
    public boolean isLayered() {
        return layered;
    }

    public List<DecisionVariable> getVariables() {
        return variables;
    }

    public List<Pair<DecisionVariable, Object>> getFixedVariables() {
        return fixedVariables;
    }

    public int getTotalDimension() {
        return totalDimension;
    }

    public int getEvaluationCount() {
        return evaluationCount;
    }

    /** Concatenated per-dimension bounds for the full optimization vector. */
    public double[][] getBounds() {
        List<double[]> b = new ArrayList<double[]>();
        for (DecisionVariable var : variables) {
            for (double[] row : var.getBounds()) {
                b.add(row);
            }
        }
        double[][] out = new double[b.size()][2];
        for (int i = 0; i < b.size(); i++) {
            out[i] = b.get(i);
        }
        return out;
    }

    /** Decode the full vector to a name -&gt; decoded-value map. */
    public Map<String, Object> decodeVariables(double[] x) {
        Map<String, Object> values = new LinkedHashMap<String, Object>();
        for (int i = 0; i < variables.size(); i++) {
            DecisionVariable var = variables.get(i);
            int offset = varOffsets[i];
            int dim = var.getDimension();
            double[] slice = new double[dim];
            System.arraycopy(x, offset, slice, 0, dim);
            values.put(var.getName(), var.decode(slice));
        }
        return values;
    }

    private void applyVariables(Object model, Map<String, Object> values) {
        for (DecisionVariable var : variables) {
            Object value = values.get(var.getName());
            if (value != null) {
                applyOne(model, var, value);
            }
        }
    }

    /** Dispatch a single variable application to the flat or LQN setter. */
    private void applyOne(Object model, DecisionVariable var, Object value) {
        if (layered) {
            if (var instanceof LqnDecisionVariable) {
                ((LqnDecisionVariable) var).applyLqn((LayeredNetwork) model, value);
            }
        } else {
            var.apply((Network) model, value);
        }
    }

    /** Deep copy the base model (flat Network or LayeredNetwork). */
    public Object copyModel() {
        if (layered) {
            LayeredNetwork copied = lqnBase.copy();
            return (copied != null) ? copied : lqnBase;
        }
        Network copied = baseModel.copy();
        return (copied != null) ? copied : baseModel;
    }

    public EvaluationResult evaluate(double[] x) {
        return evaluateValues(decodeVariables(x));
    }

    public EvaluationResult evaluateValues(Map<String, Object> values) {
        evaluationCount++;
        long start = System.nanoTime();
        EvaluationResult result = new EvaluationResult();
        try {
            Object model = copyModel();
            for (Pair<DecisionVariable, Object> fv : fixedVariables) {
                applyOne(model, fv.getLeft(), fv.getRight());
            }
            applyVariables(model, values);

            if (layered) {
                // LQN path: solve with SolverLN, read the per-node LQN table.
                LayeredNetwork lqn = (LayeredNetwork) model;
                SolverLN solver = LqnAdapter.makeLqnSolver(lqn);
                LayeredNetworkAvgTable avgTable = (LayeredNetworkAvgTable) solver.getAvgTable();
                result.feasible = true;
                result.solverUsed = "SolverLN";
                extractLayeredMetrics(avgTable, result);
                extractLayeredSystemMetrics(lqn, result);
                // LQN sensitivities are expensive (they re-solve every layer),
                // so they are computed lazily by the gradient path only.
            } else {
                Network net = (Network) model;
                // Invalidate the cached struct so the solve reflects mutations.
                net.refreshStruct(true);

                SolverAUTO solver = new SolverAUTO(net);
                NetworkAvgTable avgTable = solver.getAvgTable();
                result.feasible = true;
                try {
                    result.solverUsed = solver.getSelectedSolverName();
                } catch (RuntimeException e) {
                    result.solverUsed = "";
                }
                extractMetrics(avgTable, result);
                extractSystemMetrics(solver, net, result);
                result.sensitivities = SensitivityTable.compute(net);
            }
        } catch (RuntimeException e) {
            result.feasible = false;
        }
        result.solveTime = (System.nanoTime() - start) / 1e9;
        return result;
    }

    private void extractMetrics(NetworkAvgTable avgTable, EvaluationResult result) {
        if (avgTable == null) {
            return;
        }
        List<String> stations = avgTable.getStationNames();
        List<String> classes = avgTable.getClassNames();
        List<Double> respt = avgTable.getRespT();
        List<Double> tput = avgTable.getTput();
        List<Double> util = avgTable.getUtil();
        List<Double> qlen = avgTable.getQLen();
        int n = stations.size();
        for (int i = 0; i < n; i++) {
            String st = stations.get(i);
            String cl = classes.get(i);
            if (respt != null && i < respt.size()) {
                result.setResponseTime(st, cl, respt.get(i));
            }
            if (tput != null && i < tput.size()) {
                result.setThroughput(st, cl, tput.get(i));
            }
            if (qlen != null && i < qlen.size()) {
                result.setQueueLength(st, cl, qlen.get(i));
            }
            if (util != null && i < util.size()) {
                double prev = result.utilizations.containsKey(st) ? result.utilizations.get(st) : 0.0;
                result.utilizations.put(st, prev + util.get(i));
            }
        }
    }

    private void extractSystemMetrics(SolverAUTO solver, Network model, EvaluationResult result) {
        NetworkAvgSysTable sysTable;
        try {
            sysTable = solver.getAvgSysTable();
        } catch (RuntimeException e) {
            return;
        }
        if (sysTable == null) {
            return;
        }
        List<Double> sysRespT = sysTable.getSysRespT();
        List<Double> sysTput = sysTable.getSysTput();

        // chain index -> member class names + open flag, via the network struct
        List<List<String>> chainClasses = new ArrayList<List<String>>();
        List<Boolean> chainIsOpen = new ArrayList<Boolean>();
        try {
            NetworkStruct sn = model.getStruct();
            int nchains = sn.chains.getNumRows();
            int nclasses = sn.chains.getNumCols();
            List<String> classnames = sn.classnames;
            for (int ci = 0; ci < nchains; ci++) {
                List<String> members = new ArrayList<String>();
                boolean open = false;
                for (int k = 0; k < nclasses; k++) {
                    if (sn.chains.get(ci, k) > 0) {
                        members.add(classnames.get(k));
                        if (Double.isInfinite(sn.njobs.get(k))) {
                            open = true;
                        }
                    }
                }
                chainClasses.add(members);
                chainIsOpen.add(open);
            }
        } catch (RuntimeException e) {
            // fall back to row labels below
        }

        List<String> chainNames = sysTable.getChainNames();
        int rows = (sysTput != null) ? sysTput.size() : 0;
        for (int position = 0; position < rows; position++) {
            List<String> keys;
            if (position < chainClasses.size() && !chainClasses.get(position).isEmpty()) {
                keys = chainClasses.get(position);
            } else if (chainNames != null && position < chainNames.size()) {
                keys = new ArrayList<String>();
                keys.add(chainNames.get(position));
            } else {
                continue;
            }
            double tput = (sysTput != null && position < sysTput.size()) ? sysTput.get(position) : 0.0;
            Double respt = (sysRespT != null && position < sysRespT.size()) ? sysRespT.get(position) : null;

            // open chains: recompute end-to-end sojourn by Little's law
            if (position < chainIsOpen.size() && chainIsOpen.get(position)
                    && tput > 0 && !result.queueLengths.isEmpty()) {
                double jobsInSystem = 0.0;
                for (Map.Entry<String, Map<String, Double>> stEntry : result.queueLengths.entrySet()) {
                    for (Map.Entry<String, Double> clEntry : stEntry.getValue().entrySet()) {
                        if (keys.contains(clEntry.getKey())) {
                            jobsInSystem += clEntry.getValue();
                        }
                    }
                }
                respt = jobsInSystem / tput;
            }
            for (String key : keys) {
                if (respt != null) {
                    result.systemResponseTimes.put(key, respt);
                }
                if (sysTput != null) {
                    result.systemThroughputs.put(key, tput);
                }
            }
        }
    }

    // ---- LayeredNetwork (LQN) metric extraction ---------------------------

    /**
     * Extract per-LQN-node metrics from SolverLN's average table. One row per
     * node (Processor/Task/Entry/Activity); throughput/queue-length/response-time
     * are keyed under (node, node) and utilization under node, matching the flat
     * EvaluationResult convention. Non-finite cells are skipped.
     */
    private void extractLayeredMetrics(LayeredNetworkAvgTable avgTable, EvaluationResult result) {
        if (avgTable == null) {
            return;
        }
        List<String> nodes = avgTable.getNodeNames();
        if (nodes == null) {
            return;
        }
        List<Double> respt = avgTable.getRespT();
        List<Double> tput = avgTable.getTput();
        List<Double> util = avgTable.getUtil();
        List<Double> qlen = avgTable.getQLen();
        int n = nodes.size();
        for (int i = 0; i < n; i++) {
            String node = nodes.get(i);
            Double r = finiteCell(respt, i);
            if (r != null) {
                result.setResponseTime(node, node, r);
            }
            Double t = finiteCell(tput, i);
            if (t != null) {
                result.setThroughput(node, node, t);
            }
            Double q = finiteCell(qlen, i);
            if (q != null) {
                result.setQueueLength(node, node, q);
            }
            Double u = finiteCell(util, i);
            if (u != null) {
                result.utilizations.put(node, u);
            }
        }
    }

    private static Double finiteCell(List<Double> col, int i) {
        if (col == null || i >= col.size()) {
            return null;
        }
        Double v = col.get(i);
        if (v == null || Double.isNaN(v) || Double.isInfinite(v)) {
            return null;
        }
        return v;
    }

    /**
     * Derive end-to-end (system) metrics from the reference task(s). A closed
     * LQN's system throughput is the reference task's throughput; its end-to-end
     * response time is the sum of response times over the reference task's
     * entries. Keyed by the reference task's name.
     */
    private void extractLayeredSystemMetrics(LayeredNetwork model, EvaluationResult result) {
        for (Task task : model.getTasks().values()) {
            if (task.getScheduling() != SchedStrategy.REF) {
                continue;
            }
            String tname = task.getName();
            double tput = result.getThroughput(tname, tname);
            if (tput > 0) {
                result.systemThroughputs.put(tname, tput);
            }
            double totalRt = 0.0;
            boolean haveRt = false;
            for (Entry entry : task.getEntries()) {
                String ename = entry.getName();
                Map<String, Double> row = result.responseTimes.get(ename);
                if (row != null && row.containsKey(ename)) {
                    totalRt += row.get(ename);
                    haveRt = true;
                }
            }
            if (haveRt) {
                result.systemResponseTimes.put(tname, totalRt);
            }
        }
    }

    /**
     * Compute LQN per-layer service-rate partial sensitivities on demand.
     * Rebuilds the configured model copy, solves it with SolverLN, and reshapes
     * SolverLN.getSensitivityTable into {@code {station||jobclass: {metric: d/dRate}}}.
     * Called only by the partial-sensitivity gradient path, never on the hot
     * evaluation loop. Returns null on any failure (the caller then
     * finite-differences).
     */
    public Map<String, Map<String, Double>> evaluateLayeredSensitivities(Map<String, Object> values) {
        if (!layered) {
            return null;
        }
        try {
            Object model = copyModel();
            for (Pair<DecisionVariable, Object> fv : fixedVariables) {
                applyOne(model, fv.getLeft(), fv.getRight());
            }
            applyVariables(model, values);
            SolverLN solver = LqnAdapter.makeLqnSolver((LayeredNetwork) model);
            return LqnAdapter.computeLqnSensitivities(solver);
        } catch (RuntimeException e) {
            return null;
        }
    }

    // ---- caching -----------------------------------------------------------

    public static String valuesKey(Map<String, Object> values) {
        List<String> parts = new ArrayList<String>();
        for (Map.Entry<String, Object> e : values.entrySet()) {
            parts.add(e.getKey() + "=" + canonical(e.getValue()));
        }
        java.util.Collections.sort(parts);
        StringBuilder sb = new StringBuilder();
        for (String p : parts) {
            sb.append(p).append(';');
        }
        return sb.toString();
    }

    private static String canonical(Object v) {
        if (v instanceof double[]) {
            StringBuilder sb = new StringBuilder("[");
            for (double d : (double[]) v) {
                sb.append(Math.round(d * 1e9) / 1e9).append(',');
            }
            return sb.append(']').toString();
        }
        if (v instanceof int[]) {
            StringBuilder sb = new StringBuilder("[");
            for (int d : (int[]) v) {
                sb.append(d).append(',');
            }
            return sb.append(']').toString();
        }
        if (v instanceof Double) {
            return Double.toString(Math.round((Double) v * 1e9) / 1e9);
        }
        return String.valueOf(v);
    }

    public EvaluationResult evaluateValuesWithCache(Map<String, Object> values,
                                                    Map<String, EvaluationResult> cache) {
        String key = valuesKey(values);
        if (cache.containsKey(key)) {
            return cache.get(key);
        }
        EvaluationResult result = evaluateValues(values);
        cache.put(key, result);
        return result;
    }

    public EvaluationResult evaluateWithCache(double[] x, Map<String, EvaluationResult> cache) {
        return evaluateValuesWithCache(decodeVariables(x), cache);
    }
}
