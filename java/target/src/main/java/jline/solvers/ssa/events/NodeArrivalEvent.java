package jline.solvers.ssa.events;

import jline.lang.JobClass;
import jline.lang.nodes.Node;

public class NodeArrivalEvent extends ArrivalEvent {
    public NodeArrivalEvent(Node node, JobClass jobClass) {
        super(node, jobClass);
    }

}
