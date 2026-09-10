package jline.opt.variables;

import jline.lang.layered.LayeredNetwork;
import jline.lang.layered.Task;
import jline.opt.LqnAdapter;

import java.util.ArrayList;
import java.util.List;

/**
 * Optimize the multiplicity (thread/instance count) of an LQN Task (integer).
 * Mirrors native-Python {@code TaskMultiplicity}.
 */
public class TaskMultiplicity extends LqnDecisionVariable {

    private final String task;
    private final int minValue;
    private final int maxValue;

    public TaskMultiplicity(Task task, int minValue, int maxValue) {
        this(task.getName(), minValue, maxValue, task.getName() + "_multiplicity");
    }

    public TaskMultiplicity(String taskName, int minValue, int maxValue) {
        this(taskName, minValue, maxValue, taskName + "_multiplicity");
    }

    public TaskMultiplicity(String taskName, int minValue, int maxValue, String name) {
        super(name);
        this.task = taskName;
        this.minValue = minValue;
        this.maxValue = maxValue;
    }

    public double[][] getBounds() {
        return unitBounds(1);
    }

    public Object decode(double[] x) {
        double continuous = minValue + x[0] * (maxValue - minValue);
        long rounded = Math.round(continuous);
        if (rounded < minValue) {
            rounded = minValue;
        }
        if (rounded > maxValue) {
            rounded = maxValue;
        }
        return (int) rounded;
    }

    public void applyLqn(LayeredNetwork model, Object value) {
        Task t = LqnAdapter.resolveTask(model, task);
        if (t != null) {
            t.setMultiplicity(((Number) value).intValue());
        }
    }

    public String getVariableType() {
        return "task_multiplicity";
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
        return (t != null) ? Integer.valueOf(t.getMultiplicity()) : null;
    }
}
