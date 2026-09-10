package jline.opt.variables;

import jline.lang.layered.LayeredNetwork;
import jline.lang.layered.Task;
import jline.opt.LqnAdapter;

import java.util.ArrayList;
import java.util.List;

/**
 * Optimize the think time of an LQN Task (continuous). Mirrors native-Python
 * {@code TaskThinkTime}.
 */
public class TaskThinkTime extends LqnDecisionVariable {

    private final String task;
    private final double minValue;
    private final double maxValue;

    public TaskThinkTime(Task task, double minValue, double maxValue) {
        this(task.getName(), minValue, maxValue, task.getName() + "_thinktime");
    }

    public TaskThinkTime(String taskName, double minValue, double maxValue) {
        this(taskName, minValue, maxValue, taskName + "_thinktime");
    }

    public TaskThinkTime(String taskName, double minValue, double maxValue, String name) {
        super(name);
        this.task = taskName;
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
        Task t = LqnAdapter.resolveTask(model, task);
        if (t != null) {
            t.setThinkTime(((Number) value).doubleValue());
        }
    }

    public String getVariableType() {
        return "think_time";
    }

    public List<String> getLayer(Object model) {
        List<String> layers = new ArrayList<String>();
        layers.add(task);
        String proc = LqnAdapter.taskProcessorName((LayeredNetwork) model, task);
        if (proc != null) {
            layers.add(proc);
        }
        return layers;
    }

    public Object currentValue(Object model) {
        Task t = LqnAdapter.resolveTask((LayeredNetwork) model, task);
        return (t != null) ? Double.valueOf(t.getThinkTimeMean()) : null;
    }
}
