package jline.opt.variables;

import jline.lang.layered.Host;
import jline.lang.layered.LayeredNetwork;
import jline.lang.layered.Processor;
import jline.opt.LqnAdapter;

import java.util.ArrayList;
import java.util.List;

/**
 * Optimize the multiplicity (core count) of an LQN Processor (integer). Mirrors
 * native-Python {@code ProcessorMultiplicity}.
 */
public class ProcessorMultiplicity extends LqnDecisionVariable {

    private final String processor;
    private final int minValue;
    private final int maxValue;

    public ProcessorMultiplicity(Processor processor, int minValue, int maxValue) {
        this(processor.getName(), minValue, maxValue, processor.getName() + "_multiplicity");
    }

    public ProcessorMultiplicity(String processorName, int minValue, int maxValue) {
        this(processorName, minValue, maxValue, processorName + "_multiplicity");
    }

    public ProcessorMultiplicity(String processorName, int minValue, int maxValue, String name) {
        super(name);
        this.processor = processorName;
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
        Host proc = LqnAdapter.resolveProcessor(model, processor);
        if (proc != null) {
            proc.setMultiplicity(((Number) value).intValue());
        }
    }

    public String getVariableType() {
        return "processor_multiplicity";
    }

    public List<String> getLayer(Object model) {
        // A processor owns its own host layer, named after the processor.
        List<String> layers = new ArrayList<String>();
        layers.add(processor);
        return layers;
    }

    public Object currentValue(Object model) {
        Host proc = LqnAdapter.resolveProcessor((LayeredNetwork) model, processor);
        return (proc != null) ? Integer.valueOf(proc.getMultiplicity()) : null;
    }
}
