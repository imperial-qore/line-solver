package jline.opt;

import jline.VerboseLevel;
import jline.lang.constant.SolverType;
import jline.lang.layered.Activity;
import jline.lang.layered.Entry;
import jline.lang.layered.Host;
import jline.lang.layered.LayeredNetwork;
import jline.lang.layered.Processor;
import jline.lang.layered.Task;
import jline.solvers.LayeredNetworkSensitivityTable;
import jline.solvers.SolverOptions;
import jline.solvers.ln.SolverLN;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * LayeredNetwork (LQN) support for line-opt. Mirrors native-Python
 * {@code line_solver.opt.layered}.
 *
 * <p>line-opt's core ({@link OptimizationProblem}/{@link LineEvaluator}/
 * decision variables/{@link jline.opt.solver.LineOptSolver}) was written for the
 * flat {@link jline.lang.Network}: it copies the model, solves it with
 * SolverAUTO, and reads a per-(station, class) average table. A
 * {@link LayeredNetwork} cannot go through that path -- it is an ensemble of
 * per-layer submodels, solved by {@link SolverLN} into a per-LQN-node average
 * table (Node, NodeType, QLen, Util, RespT, ResidT, ArvR, Tput). This class is
 * the adapter that lets the same optimizer drive an LQN: element resolution by
 * name in a per-evaluation model copy, a quiet SolverLN factory, and a reshape
 * of SolverLN's per-layer service-rate sensitivity table.</p>
 *
 * <p>IMPORTANT (see {@link LayeredNetworkSensitivityTable}): the per-layer table
 * holds WITHIN-LAYER PARTIAL derivatives, taken with the fixed-point layer
 * parameters held constant. It omits the cross-layer coupling term, so it is a
 * biased estimate of the total derivative of the solved layered model. The
 * optimizer's {@code lqn_gradient='fd'} mode finite-differences the whole
 * LayeredNetwork instead (correct total derivative); {@code partial_sens} uses
 * this table directly (cheap, biased); {@code partial_plus_fd} corrects it with
 * a periodic full-model finite difference.</p>
 */
public final class LqnAdapter {

    private LqnAdapter() {
    }

    /** True if the model is a LayeredNetwork (LQN), false for a flat Network. */
    public static boolean isLayered(Object model) {
        return model instanceof LayeredNetwork;
    }

    /** Resolve an Activity by name inside a (possibly copied) LQN model. */
    public static Activity resolveActivity(LayeredNetwork model, String name) {
        if (model == null) {
            return null;
        }
        for (Activity a : model.getActivities().values()) {
            if (a.getName().equals(name)) {
                return a;
            }
        }
        return null;
    }

    /** Resolve a Task by name inside a (possibly copied) LQN model. */
    public static Task resolveTask(LayeredNetwork model, String name) {
        if (model == null) {
            return null;
        }
        for (Task t : model.getTasks().values()) {
            if (t.getName().equals(name)) {
                return t;
            }
        }
        return null;
    }

    /**
     * Resolve a Processor by name inside a (possibly copied) LQN model. The LQN
     * has no getProcessors(); processors are the {@link Processor} entries among
     * the hosts, so this matches by name over {@code getHosts()}.
     */
    public static Host resolveProcessor(LayeredNetwork model, String name) {
        if (model == null) {
            return null;
        }
        for (Host h : model.getHosts().values()) {
            if (h.getName().equals(name)) {
                return h;
            }
        }
        return null;
    }

    /**
     * Name of the processor an activity ultimately runs on
     * (Activity -&gt; parent Task -&gt; deployed Processor), or null if the chain
     * is incomplete. Used both to tag a HostDemand variable with the host layer
     * it perturbs and to key its row in the per-layer sensitivity table, whose
     * host-layer rows are (Layer=processor, Station=processor, JobClass=activity).
     */
    public static String activityProcessorName(LayeredNetwork model, String activityName) {
        Activity act = resolveActivity(model, activityName);
        if (act == null) {
            return null;
        }
        Task task = act.getParent();
        if (task == null) {
            return null;
        }
        Processor proc = task.getProcessor();
        return (proc != null) ? proc.getName() : null;
    }

    /** Name of the processor a task is deployed on (null if undeployed). */
    public static String taskProcessorName(LayeredNetwork model, String taskName) {
        Task task = resolveTask(model, taskName);
        if (task == null) {
            return null;
        }
        Processor proc = task.getProcessor();
        return (proc != null) ? proc.getName() : null;
    }

    /**
     * Construct a quiet SolverLN for a per-evaluation LQN model copy. Silences
     * the layer/avg-table printing so an optimization run of thousands of solves
     * does not flood stdout.
     */
    public static SolverLN makeLqnSolver(LayeredNetwork model) {
        SolverOptions options = new SolverOptions(SolverType.LN);
        options.verbose = VerboseLevel.SILENT;
        return new SolverLN(model, options);
    }

    /**
     * Reshape SolverLN's per-layer service-rate sensitivity table into
     * {@code {station||jobclass: {Tput,RespT,QLen,Util}}}, each value the
     * derivative of that layer row's mean measure with respect to that
     * station-class service RATE. Returns null if the table is unavailable. See
     * the class docstring for the partial-vs-total caveat.
     */
    public static Map<String, Map<String, Double>> computeLqnSensitivities(SolverLN solver) {
        LayeredNetworkSensitivityTable table;
        try {
            table = solver.getSensitivityTable();
        } catch (RuntimeException e) {
            return null;
        }
        if (table == null) {
            return null;
        }
        List<String> stations = table.getStationNames();
        List<String> classes = table.getClassNames();
        if (stations == null || classes == null) {
            return null;
        }
        List<Double> dTput = table.getDTput();
        List<Double> dRespT = table.getDRespT();
        List<Double> dQLen = table.getDQLen();
        List<Double> dUtil = table.getDUtil();

        Map<String, Map<String, Double>> out = new LinkedHashMap<String, Map<String, Double>>();
        int n = stations.size();
        for (int i = 0; i < n; i++) {
            String key = stations.get(i) + "||" + classes.get(i);
            Map<String, Double> row = out.get(key);
            if (row == null) {
                row = new LinkedHashMap<String, Double>();
                out.put(key, row);
            }
            row.put("Tput", cell(dTput, i));
            row.put("RespT", cell(dRespT, i));
            row.put("QLen", cell(dQLen, i));
            row.put("Util", cell(dUtil, i));
        }
        return out;
    }

    private static double cell(List<Double> col, int i) {
        if (col == null || i >= col.size()) {
            return 0.0;
        }
        Double v = col.get(i);
        return (v == null || Double.isNaN(v)) ? 0.0 : v;
    }
}
