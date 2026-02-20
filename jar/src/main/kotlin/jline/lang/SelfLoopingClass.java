/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang;

import jline.lang.nodes.Station;

import java.io.Serializable;

/**
 * Class of jobs that perpetually loop at a given station
 */
public class SelfLoopingClass extends ClosedClass implements Serializable {

    /**
     * Creates a new self-looping job class with specified priority.
     * Jobs in this class perpetually loop at the reference station.
     * 
     * @param model the network model to add this class to
     * @param name the name for this self-looping class
     * @param njobs the number of jobs in this closed class
     * @param refstat the reference station where jobs loop
     * @param priority the priority level for this class
     */
    public SelfLoopingClass(Network model, String name, long njobs, Station refstat, int priority) {
        super(model, name, njobs, refstat, priority);
    }

    /**
     * Creates a new self-looping job class with default priority (0).
     * Jobs in this class perpetually loop at the reference station.
     * 
     * @param model the network model to add this class to
     * @param name the name for this self-looping class
     * @param njobs the number of jobs in this closed class
     * @param refstat the reference station where jobs loop
     */
    public SelfLoopingClass(Network model, String name, long njobs, Station refstat) {
        this(model, name, njobs, refstat, 0);
    }

}
