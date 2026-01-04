/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.nodes;

import jline.lang.Network;
import jline.lang.constant.SchedStrategy;

import java.io.Serializable;

/**
 * An infinite server station, i.e. a node imposing a delay without queueing to an incoming job
 */
public class Delay extends Queue implements Serializable {
    /**
     * Creates a new delay station with infinite servers.
     * Jobs experience service delay without queueing (infinite server scheduling).
     * 
     * @param model the network model to add this delay station to
     * @param name the name for this delay station
     */
    public Delay(Network model, String name) {
        super(model, name, SchedStrategy.INF);
        this.numberOfServers = Integer.MAX_VALUE;
    }
}
