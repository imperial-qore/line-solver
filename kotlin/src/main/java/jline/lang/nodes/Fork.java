package jline.lang.nodes;

import jline.lang.JobClass;
import jline.lang.Network;
import jline.lang.constant.SchedStrategy;
import jline.lang.sections.Buffer;
import jline.lang.sections.Forker;
import jline.lang.sections.ServiceTunnel;

import java.io.Serializable;
import java.util.List;

/**
 * A node that forks an incoming job into a set of sibling tasks
 */
public class Fork extends Node implements Serializable {
    private final double cap;
    private final SchedStrategy schedStrategy;

    public Fork(Network model){
        this(model, "Fork");
    }

    public Fork(Network model, String name) {
        super(name);
        List<JobClass> classes = model.getClasses();
        this.cap = Double.POSITIVE_INFINITY;
        this.input = new Buffer(classes);
        this.schedStrategy = SchedStrategy.FORK;
        this.server = new ServiceTunnel();
        this.output = new Forker(classes);
        this.setModel(model);
        model.addNode(this);
    }

    /**
     * Sets the number of tasks sent out on each outgoing link
     * @param nTasks - the number of tasks
     */
    public void setTasksPerLink(int nTasks){
        Forker f = (Forker) this.output;
        f.taskPerLink = nTasks;
    }

    @Override
    public Network getModel() {
        return this.model;
    }
}
