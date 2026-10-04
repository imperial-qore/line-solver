/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang;

import jline.lang.constant.JobClassType;
import jline.lang.constant.JoinStrategy;
import jline.lang.constant.RoutingStrategy;
import jline.lang.nodes.Join;
import jline.lang.nodes.Node;
import jline.lang.nodes.Station;

import java.io.Serializable;

/**
 * Class of jobs that perpetually loop at a given station
 */
public class DisabledClass extends JobClass implements Serializable {

    /**
     * Creates a new disabled job class that perpetually loops at the reference station.
     * Configures all nodes in the network to use disabled routing for this class.
     * 
     * @param model the network model to add this class to
     * @param name the name for this disabled class
     * @param refstat the reference station where jobs of this class loop
     * @throws Exception if there's an error setting up the disabled class
     */
    public DisabledClass(Network model, String name, Station refstat) throws Exception {
        super(JobClassType.DISABLED, name);
        this.type = JobClassType.DISABLED;
        this.priority = 0;
        model.addJobClass(this);
        this.setReferenceStation(refstat);
        for (int i = 0; i < model.getNumberOfNodes(); i++) {
            Node node_i = model.getNodeByIndex(i);
            node_i.setRouting(this, RoutingStrategy.DISABLED);
            if (node_i instanceof Join) {
                ((Join) node_i).setStrategy(this, JoinStrategy.STD);
                ((Join) node_i).setRequired(this, -1);
            }
        }
    }

}
