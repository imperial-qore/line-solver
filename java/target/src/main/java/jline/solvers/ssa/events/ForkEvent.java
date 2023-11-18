package jline.solvers.ssa.events;

import jline.lang.JobClass;
import jline.lang.OutputStrategy;
import jline.lang.nodes.Node;
import jline.lang.sections.OutputSection;

import java.util.List;

public class ForkEvent extends SynchedEvent {
    public List<OutputStrategy> outputStrategies;
    public ForkEvent(OutputSection outputSection, Node targetNode, JobClass jobClass) {
        super(outputSection, targetNode, jobClass);
        this.outputStrategies = outputSection.getOutputStrategies();
    }

}
