package jline.lang;

import jline.lang.constant.JobClassType;
import jline.lang.constant.JoinStrategy;
import jline.lang.nodes.Join;
import jline.lang.nodes.Node;
import jline.lang.nodes.Station;

import java.io.Serializable;

/**
 * Class where jobs perpetually loop without arriving or leaving (Closed class)
 */
public class ClosedClass extends JobClass implements Serializable {
    protected double population;
    protected int classIndex;
    protected Network model;
    public ClosedClass(Network model, String name, double njobs, Station refstat, int priority) {
        super(JobClassType.CLOSED, name);
        this.index = model.getNumberOfClasses()+1;
        model.addJobClass(this);
        this.population = njobs;
        for (int i = 0; i < model.getNumberOfNodes(); i++) {
            Node currentNode = model.getNodes().get(i);
//            model.setNodeRouting(i, this, RoutingStrategy.RAND);
            if (currentNode instanceof Join){
                model.setJoinNodeStrategy(i, this, JoinStrategy.STD);
                model.setJoinNodeRequired(i, this, -1);
            }
        }
        this.model = model;
        this.refstat = refstat;
        this.classIndex = -1;
        this.setPriority(priority);
    }
    public ClosedClass(Network model, String name, double njobs, Station refstat) {
        this(model, name,njobs, refstat,0);
    }

    @Override
    public void printSummary() {
        System.out.format("Closed class: %s\n", this.getName());
    }

    @Override
    public double getNumberOfJobs() {
        return population;
    }

    public double getPopulation() {
        return this.population;
    }

    @Override
    public int getJobClassIdx() {
        if (this.classIndex == -1) {
            this.classIndex = this.model.getJobClassIndex(this);
        }
        return this.classIndex;
    }

    public void setPopulation(double pop){
        this.population = pop;
    }
}
