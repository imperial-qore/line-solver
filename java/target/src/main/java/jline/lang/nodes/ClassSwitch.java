package jline.lang.nodes;

import jline.lang.JobClass;
import jline.lang.Network;
import jline.lang.constant.RoutingStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SchedStrategyType;
import jline.lang.sections.Buffer;
import jline.lang.sections.ClassSwitchOutputSection;
import jline.lang.sections.StatelessClassSwitcher;
import jline.solvers.ssa.events.ArrivalEvent;
import jline.solvers.ssa.events.ClassSwitchArrivalEvent;
import jline.util.Matrix;

import java.io.Serializable;
import java.util.List;

/**
 * A node that switches the class of an incoming job based on a probability table
 */
public class ClassSwitch extends Node implements Serializable {
    protected SchedStrategyType schedPolicy;
    protected SchedStrategy schedStrategy;
    public boolean autoAdded;

    public ClassSwitch(Network model, String name) {
        super(name);

        List<JobClass> jobClasses = model.getClasses();
        int m = jobClasses.size();

        this.input = new Buffer(jobClasses);
        this.output = new ClassSwitchOutputSection(jobClasses);

        Matrix csFun = new Matrix(m , m, m*m);

        this.setModel(model);
        this.model.addNode(this);

        this.schedPolicy = SchedStrategyType.NP;
        this.schedStrategy = SchedStrategy.FCFS;
        this.server = new StatelessClassSwitcher(jobClasses, csFun);
        this.autoAdded = false;
    }

    public ClassSwitch(Network model, String name, Matrix csFun) {
        super(name);

        List<JobClass> jobClasses = model.getClasses();
        this.input = new Buffer(jobClasses);
        this.output = new ClassSwitchOutputSection(jobClasses);


        this.setModel(model);
        this.model.addNode(this);

        this.schedPolicy = SchedStrategyType.NP;
        this.schedStrategy = SchedStrategy.FCFS;
        this.server = new StatelessClassSwitcher(jobClasses, csFun);
    }

    public Matrix initClassSwitchMatrix() {
        int m = this.model.getNumberOfClasses();
        return new Matrix(m,m,m*m);
    }

    public void setClassSwitchingMatrix(Matrix csMatrix) {
        ((StatelessClassSwitcher)this.server).updateClasses(this.model.getClasses());
        ((StatelessClassSwitcher)this.server).updateClassSwitch(csMatrix);
    }

    public void setProbRouting(JobClass jobClass, Node destination, double probability) {
        this.setRouting(jobClass, RoutingStrategy.PROB, destination, probability);
    }

    @Override
    public ArrivalEvent getArrivalEvent(JobClass jobClass) {
        if (!this.arrivalEvents.containsKey(jobClass)) {
            this.arrivalEvents.put(jobClass, new ClassSwitchArrivalEvent(this, jobClass, (StatelessClassSwitcher) this.server));
        }
        return this.arrivalEvents.get(jobClass);
    }
}