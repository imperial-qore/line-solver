/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.nodes;

import jline.lang.Network;
import jline.lang.constant.RoutingStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.sections.ServiceSection;

import java.io.Serializable;

/**
 * An abstraction of the external world jobs in open classes depart to
 */
public class Sink extends Node implements Serializable {
    protected SchedStrategy schedStrategy;

    /**
     * Creates a new sink node where jobs in open classes depart from the network.
     * Configures disabled routing for all job classes.
     * 
     * @param model the network model to add this sink to (can be null for standalone creation)
     * @param name the name for this sink node
     */
    public Sink(Network model, String name) {
        super(name);


        if (model != null) {
            this.setModel(model);

            this.server = new ServiceSection("JobSink");
            this.setModel(model);
            this.model.addNode(this);
            this.schedStrategy = SchedStrategy.EXT;
            if (this.model.hasClasses()) {
                for (int r = 0; r < this.model.getNumberOfClasses(); r++) {
                    this.setRouting(this.model.getClassByIndex(r), RoutingStrategy.DISABLED);
                }
            }
        }
    }

    @Override
    public void printSummary() {
        System.out.format("jline.Sink: %s\n", this.getName());
    }

}
