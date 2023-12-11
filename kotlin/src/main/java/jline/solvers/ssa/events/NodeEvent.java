package jline.solvers.ssa.events;
import jline.lang.nodes.Node;

public interface NodeEvent {
    Node getNode();
    int getNodeStatefulIdx();
    int getClassIdx();
    boolean isStateful();
}
