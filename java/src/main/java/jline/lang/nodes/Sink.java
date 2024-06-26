package jline.lang.nodes;

import jline.lang.sections.ServiceSection;

import java.io.Serializable;
import jline.lang.*;
import jline.lang.constant.SchedStrategy;

/**
 * An abstraction of the external world jobs in open classes depart to
 */
public class Sink extends Node implements Serializable {
    protected SchedStrategy schedStrategy;

    public Sink(Network model, String name) {
        super(name);


        if (model != null) {
            this.setModel(model);

            this.server = new ServiceSection("JobSink");
            this.model.addNode(this);
            this.setModel(model);
            this.schedStrategy = SchedStrategy.EXT;
        }
    }

    @Override
    public void printSummary() {
        System.out.format("jline.Sink: %s\n", this.getName());
    }

}
