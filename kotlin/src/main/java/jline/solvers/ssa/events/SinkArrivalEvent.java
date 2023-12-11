package jline.solvers.ssa.events;

import jline.lang.nodes.Node;

public class SinkArrivalEvent extends ArrivalEvent{
    public SinkArrivalEvent(Node node) {
        super(node, null);
    }


}
