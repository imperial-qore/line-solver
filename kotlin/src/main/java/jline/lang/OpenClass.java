package jline.lang;

import java.io.Serializable;

import jline.lang.constant.JobClassType;
import jline.lang.constant.JoinStrategy;
import jline.lang.constant.RoutingStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Join;
import jline.lang.nodes.Node;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;

/**
 * A class of jobs that arrives from the external world to the Network and, after completion, leaves it
 */
public class OpenClass extends JobClass  implements Serializable {
    protected int classIndex;
    protected Network model;

    public OpenClass(Network model, String name, int priority) {
        super(JobClassType.Open, name);
        this.index = model.getNumberOfClasses()+1;
        model.addJobClass(this);
        try {
            setReference(model.getSource());
        } catch (Exception e){
            System.err.println("The model requires a Source prior to instantiating open classes.");
        }
        this.classIndex = -1;

        for (int i = 0; i < model.getNumberOfNodes(); i++) {
            Node currentNode = model.getNodes().get(i);
            if (currentNode != null && !(currentNode instanceof Source) && !(currentNode instanceof Sink) && !(currentNode instanceof Join)){
                model.setNodeScheduling(i, this, SchedStrategy.FCFS);
            }
            if (currentNode instanceof Join){
                model.setJoinNodeStrategy(i, this, JoinStrategy.STD);
                model.setJoinNodeRequired(i, this, -1);
            }
//            if (currentNode != null){
//                model.setNodeRouting(i, this, RoutingStrategy.RAND);
//            }
        }
        this.model = model;
        this.setPriority(priority);
    }
    public OpenClass(Network model, String name) {
        this(model, name,0);
    }

    @Override
    public void setReference(Node source) throws Exception {
        if (!(source instanceof Source)) {
            throw new Exception("The reference station for an open class must be a jline.Source.");
        }
        super.setReference(source);
    }

    @Override
    public void printSummary() {
        System.out.format("Open class: %s\n", this.getName());
    }

    @Override
    public int getJobClassIdx() {
        if (this.classIndex == -1) {
            this.classIndex = this.model.getJobClassIndex(this);
        }
        return this.classIndex;
    }
}
