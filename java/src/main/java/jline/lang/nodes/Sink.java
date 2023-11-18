package jline.lang.nodes;

import jline.lang.sections.ServiceSection;
import jline.solvers.ssa.events.ArrivalEvent;
import jline.solvers.ssa.events.SinkArrivalEvent;

import java.io.Serializable;
import jline.lang.*;
import jline.lang.constant.SchedStrategy;

/**
 * An abstraction of the external world jobs in open classes depart to
 */
public class Sink extends Node implements Serializable {
    protected SchedStrategy schedStrategy;
    protected jline.solvers.ssa.events.ArrivalEvent arrivalEvent;

    public Sink(Network model, String name) {
        super(name);


        if (model != null) {
            this.setModel(model);

            this.server = new ServiceSection("JobSink");
            this.model.addNode(this);
            this.setModel(model);
            this.schedStrategy = SchedStrategy.EXT;
        }

        this.arrivalEvent = new SinkArrivalEvent(this);
    }

    @Override
    public void printSummary() {
        System.out.format("jline.Sink: %s\n", this.getName());
    }


    @Override
    public ArrivalEvent getArrivalEvent(JobClass jobClass) {
        return this.arrivalEvent;
    }
}
