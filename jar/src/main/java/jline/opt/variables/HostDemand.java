package jline.opt.variables;

import jline.lang.layered.Activity;
import jline.lang.layered.LayeredNetwork;
import jline.opt.LqnAdapter;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Optimize the mean host demand of an LQN Activity (continuous). The host demand
 * D is the mean service requirement the activity places on its processor; the
 * processor-layer service rate is mu = 1/D. This is the primary LQN tuning knob
 * (analogous to {@link ServiceRate} for a flat station). It exposes the hooks
 * the partial-sensitivity gradient path needs ({@link #sensKey},
 * {@link #sensMetricTargets}, {@link #rateJacobian}, {@link #decodeJacobian}).
 * Mirrors native-Python {@code HostDemand}.
 */
public class HostDemand extends LqnDecisionVariable {

    private final String activity;
    private final double minDemand;
    private final double maxDemand;

    public HostDemand(Activity activity, double minDemand, double maxDemand) {
        this(activity.getName(), minDemand, maxDemand, activity.getName() + "_hostdemand");
    }

    public HostDemand(Activity activity, double minDemand, double maxDemand, String name) {
        this(activity.getName(), minDemand, maxDemand, name);
    }

    public HostDemand(String activityName, double minDemand, double maxDemand) {
        this(activityName, minDemand, maxDemand, activityName + "_hostdemand");
    }

    public HostDemand(String activityName, double minDemand, double maxDemand, String name) {
        super(name);
        this.activity = activityName;
        this.minDemand = minDemand;
        this.maxDemand = maxDemand;
    }

    public String getActivity() {
        return activity;
    }

    public double[][] getBounds() {
        return unitBounds(1);
    }

    public Object decode(double[] x) {
        return minDemand + x[0] * (maxDemand - minDemand);
    }

    public void applyLqn(LayeredNetwork model, Object value) {
        Activity act = LqnAdapter.resolveActivity(model, activity);
        if (act != null) {
            act.setHostDemand(((Number) value).doubleValue());
        }
    }

    public String getVariableType() {
        return "host_demand";
    }

    public List<String> getLayer(Object model) {
        String proc = LqnAdapter.activityProcessorName((LayeredNetwork) model, activity);
        if (proc == null) {
            return null;
        }
        List<String> layers = new ArrayList<String>();
        layers.add(proc);
        return layers;
    }

    public Object currentValue(Object model) {
        Activity act = LqnAdapter.resolveActivity((LayeredNetwork) model, activity);
        return (act != null) ? Double.valueOf(act.getHostDemandMean()) : null;
    }

    // ---- partial-sensitivity gradient hooks --------------------------------

    /**
     * Row key {@code station||jobclass} in the per-layer sensitivity table. The
     * host-layer rows are (Layer=processor, Station=processor,
     * JobClass=activity); the sensitivity value is d(metric)/d(service rate).
     * Returns null if the activity is not deployed on a processor.
     */
    public String sensKey(LayeredNetwork model) {
        String proc = LqnAdapter.activityProcessorName(model, activity);
        if (proc == null) {
            return null;
        }
        return proc + "||" + activity;
    }

    /**
     * Map each layer-row metric to the EvaluationResult key it approximates. The
     * host-layer row's utilization tracks the processor node's utilization; its
     * throughput/queue-length/response-time track the activity node's. Util is
     * keyed by node name, the others by {@code node||node}.
     */
    public Map<String, String> sensMetricTargets(LayeredNetwork model) {
        String proc = LqnAdapter.activityProcessorName(model, activity);
        Map<String, String> targets = new LinkedHashMap<String, String>();
        targets.put("Util", proc);
        targets.put("Tput", activity + "||" + activity);
        targets.put("QLen", activity + "||" + activity);
        targets.put("RespT", activity + "||" + activity);
        return targets;
    }

    /** d(service rate)/d(demand) = d(1/D)/dD = -1/D^2 at demand D=value. */
    public double rateJacobian(double value) {
        if (value <= 0) {
            return 0.0;
        }
        return -1.0 / (value * value);
    }

    /** d(decoded demand)/d(encoded x): constant slope of the linear map. */
    public double decodeJacobian(double[] x) {
        return maxDemand - minDemand;
    }
}
