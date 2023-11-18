package jline.solvers.ssa.events;

import jline.lang.JobClass;
import jline.lang.OutputStrategy;
import jline.lang.nodes.Node;
import jline.lang.sections.OutputSection;


import java.util.*;

public class JoinEvent extends SynchedEvent {
    public int nSources;
    public Map<OutputSection, Integer> waitingJobs;
    public Set<OutputSection> seenSet;
    public JoinEvent(OutputSection outputSection, Node targetNode, JobClass jobClass) {
        super(outputSection, targetNode, jobClass);

        Node corresNode = null;

        for (Node node : targetNode.getModel().getNodes()) {
            if (node.getOutput() == outputSection) {
                corresNode = node;
            }
        }

        if (corresNode == null) {
            return;
        }

        this.waitingJobs = new HashMap<OutputSection, Integer>();
        this.seenSet = new HashSet<OutputSection>();

        for (Node node : targetNode.getModel().getNodes()) {
            for (OutputStrategy os : node.getOutputStrategies()) {
                if (os.getDestination() == corresNode) {
                    this.nSources += 1;
                    this.waitingJobs.put(node.getOutput(), 0);
                    break;
                }
            }
        }
    }

}