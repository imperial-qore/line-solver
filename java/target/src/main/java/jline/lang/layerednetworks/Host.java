package jline.lang.layerednetworks;

import jline.lang.constant.SchedStrategy;

import java.util.ArrayList;
import java.util.List;

/**
 * A processor that can run Tasks
 */
public class Host extends LayeredNetworkElement {
    protected int multiplicity;
    protected int replication;
    protected SchedStrategy scheduling;
    protected double quantum;
    protected double speedFactor;
    protected List<Task> tasks;
    private int ID;

    public Host(LayeredNetwork model, String name, int multiplicity, SchedStrategy scheduling, double quantum, double speedFactor) {
        super(name);
        this.multiplicity = multiplicity;
        this.replication = 1;
        this.scheduling = scheduling;
        this.quantum = quantum;
        this.speedFactor = speedFactor;
        this.model = model;
        this.tasks = new ArrayList<>();
        model.hosts.put(model.hosts.size(),this);//TODO
    }

    public Host(LayeredNetwork model, String name, int multiplicity, SchedStrategy scheduling, double quantum) {
        super(name);
        this.multiplicity = multiplicity;
        this.replication = 1;
        this.scheduling = scheduling;
        this.quantum = quantum;
        this.speedFactor = 1;
        this.model = model;
        this.tasks = new ArrayList<>();
        model.hosts.put(model.hosts.size(),this);
    }

    public Host(LayeredNetwork model, String name, int multiplicity, SchedStrategy scheduling) {
        super(name);
        this.multiplicity = multiplicity;
        this.replication = 1;
        this.scheduling = scheduling;
        this.quantum = 0.01;
        this.speedFactor = 1;
        this.model = model;
        this.tasks = new ArrayList<>();
        model.hosts.put(model.hosts.size(),this);
    }

    public Host(LayeredNetwork model, String name, int multiplicity) {
        super(name);
        this.multiplicity = multiplicity;
        this.replication = 1;
        this.scheduling = SchedStrategy.PS;
        this.quantum = 0.01;
        this.speedFactor = 1;
        this.model = model;
        this.tasks = new ArrayList<>();
        model.hosts.put(model.hosts.size(),this);
    }

    public Host(LayeredNetwork model, String name) {
        super(name);
        this.multiplicity = 1;
        this.replication = 1;
        this.scheduling = SchedStrategy.PS;
        this.quantum = 0.01;
        this.speedFactor = 1;
        this.model = model;
        this.tasks = new ArrayList<>();
        model.hosts.put(model.hosts.size(),this);
    }

    public void setReplication(int replication) {
        this.replication = replication;
    }

    public void addTask(Task newTask){
        tasks.add(newTask);
    }
}
