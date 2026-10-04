package jline.opt.variables;

import jline.lang.layered.Activity;
import jline.lang.layered.LayeredNetwork;
import jline.lang.layered.Task;
import jline.opt.LqnAdapter;

import java.util.ArrayList;
import java.util.List;

/**
 * Optimize the activity-level think time of an LQN Activity (continuous).
 * Mirrors native-Python {@code ActivityThinkTime}.
 */
public class ActivityThinkTime extends LqnDecisionVariable {

    private final String activity;
    private final double minValue;
    private final double maxValue;

    public ActivityThinkTime(Activity activity, double minValue, double maxValue) {
        this(activity.getName(), minValue, maxValue, activity.getName() + "_thinktime");
    }

    public ActivityThinkTime(String activityName, double minValue, double maxValue) {
        this(activityName, minValue, maxValue, activityName + "_thinktime");
    }

    public ActivityThinkTime(String activityName, double minValue, double maxValue, String name) {
        super(name);
        this.activity = activityName;
        this.minValue = minValue;
        this.maxValue = maxValue;
    }

    public double[][] getBounds() {
        return unitBounds(1);
    }

    public Object decode(double[] x) {
        return minValue + x[0] * (maxValue - minValue);
    }

    public void applyLqn(LayeredNetwork model, Object value) {
        Activity act = LqnAdapter.resolveActivity(model, activity);
        if (act != null) {
            act.setThinkTime(((Number) value).doubleValue());
        }
    }

    public String getVariableType() {
        return "think_time";
    }

    public List<String> getLayer(Object model) {
        LayeredNetwork lqn = (LayeredNetwork) model;
        Activity act = LqnAdapter.resolveActivity(lqn, activity);
        Task task = (act != null) ? act.getParent() : null;
        String proc = LqnAdapter.activityProcessorName(lqn, activity);
        List<String> layers = new ArrayList<String>();
        if (task != null) {
            layers.add(task.getName());
        }
        if (proc != null) {
            layers.add(proc);
        }
        return layers.isEmpty() ? null : layers;
    }

    public Object currentValue(Object model) {
        Activity act = LqnAdapter.resolveActivity((LayeredNetwork) model, activity);
        return (act != null) ? Double.valueOf(act.getThinkTimeMean()) : null;
    }
}
