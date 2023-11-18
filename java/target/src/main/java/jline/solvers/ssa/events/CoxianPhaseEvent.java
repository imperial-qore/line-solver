package jline.solvers.ssa.events;

import jline.lang.HasSchedStrategy;
import jline.lang.JobClass;
import jline.lang.constant.SchedStrategy;
import jline.lang.distributions.Distribution;
import jline.lang.distributions.Coxian;
import jline.lang.nodes.Node;
import jline.lang.nodes.Source;
import jline.lang.nodes.StatefulNode;
import jline.solvers.ssa.state.SSAStateMatrix;
import jline.util.Matrix;

public class CoxianPhaseEvent extends PhaseEvent implements NodeEvent {
    public final int statefulIndex;
    public final int classIndex;
    public final SchedStrategy schedStrategy;
    public final boolean isSource;

    public Node node;
    public final JobClass jobClass;
    public boolean isProcessorSharing;
    public Matrix subGenerator;

    public final DepartureEvent departureEvent;


    public final Coxian serviceProcess;
    @SuppressWarnings("unchecked")
    public CoxianPhaseEvent(Node node, JobClass jobClass, DepartureEvent departureEvent) {
        super();
        this.node = node;
        this.jobClass = jobClass;

        if (node instanceof StatefulNode) {
            this.statefulIndex = ((StatefulNode)this.node).getStatefulIndex();
        } else {
            this.statefulIndex = -1;
        }
        this.classIndex = this.node.getModel().getJobClassIndex(this.jobClass);

        this.isSource = node instanceof Source;
        if (!(node instanceof HasSchedStrategy)) {
            throw new RuntimeException("Scheduling strategy required");
        }

        this.schedStrategy = ((HasSchedStrategy)node).getSchedStrategy();

        Distribution distServiceProcess = ((HasSchedStrategy)node).getServiceProcess(this.jobClass);
        if (!(distServiceProcess instanceof Coxian)) {
            throw new RuntimeException("Coxian distribution required");
        }
        this.serviceProcess = (Coxian)distServiceProcess;


        this.departureEvent = departureEvent;
        this.isProcessorSharing = this.schedStrategy == SchedStrategy.PS;

        this.subGenerator = this.serviceProcess.getRepres().get(0);
    }

    @Override
    public long getNPhases() {
        return this.subGenerator.length();
    }

    @Override
    public double getRate(SSAStateMatrix networkState) {
        if (this.isProcessorSharing) {
            double totalRate = 0;
            for (int i = 0; i < this.subGenerator.length(); i++) {
                int inPhase = networkState.getInPhase(this.statefulIndex, this.classIndex, i);
                totalRate += inPhase * this.serviceProcess.getTotalPhaseRate(i);
            }

            double serviceRatio = (double) networkState.getState(this.statefulIndex, this.classIndex)/(double) networkState.totalStateAtNode(this.statefulIndex);
            serviceRatio *= networkState.psTotalCapacity(this.statefulIndex);
            return totalRate * serviceRatio;
        }

        int activeServers = 1;
        if (this.node instanceof StatefulNode) {
            activeServers = networkState.inProcess(this.statefulIndex, this.classIndex);
            if (this.node instanceof Source) {
                // NOTE: Pay active attention to this part
                activeServers = 1; //stateMatrix.getPhaseListSize(this)+1;
            } else if (activeServers == 0) {
                return Double.NaN;
            }
        }

        double totalRate = 0;

        for (int i = 0; i < this.subGenerator.length(); i++) {
            int inPhase = networkState.getInPhase(this.statefulIndex, this.classIndex, i);
            totalRate += inPhase * this.serviceProcess.getTotalPhaseRate(i);
        }

        return totalRate;
    }


    public Node getNode() {
        return this.node;
    }

    public int getNodeStatefulIdx() {
        return this.statefulIndex;
    }
    public int getClassIdx() {
        return this.classIndex;
    }

    public boolean isStateful() {
        return this.statefulIndex != -1;
    }
}
