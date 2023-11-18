package jline.solvers.ssa.events;

import jline.lang.JobClass;
import jline.lang.nodes.Node;
import jline.lang.nodes.StatefulNode;
import jline.lang.processes.MAP;
import jline.solvers.ssa.state.SSAStateMatrix;

public class MAPPhaseEvent extends PhaseEvent {
    public MAP MAP;

    public final int statefulIndex;
    public final int classIndex;

    public MAPPhaseEvent(Node node, JobClass jobClass, MAP MAP) {
        super();
        this.MAP = MAP;

        if (node instanceof StatefulNode) {
            this.statefulIndex = ((StatefulNode)node).getStatefulIndex();
        } else {
            this.statefulIndex = -1;
        }
        this.classIndex = node.getModel().getJobClassIndex(jobClass);
    }

    @Override
    public long getNPhases() {
        return MAP.getNumberOfPhases();
    }

    @Override
    public double getRate(SSAStateMatrix networkState) {
        return (long)MAP.getTotalPhaseRate(networkState.getGlobalPhase(this.statefulIndex, this.classIndex));
    }

    public boolean updateGlobalPhase(int classIdx, int newPhase) {
        throw new RuntimeException("Not implemented");
    }
}
