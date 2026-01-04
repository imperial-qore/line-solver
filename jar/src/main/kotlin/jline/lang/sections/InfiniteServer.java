/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.sections;

import static jline.GlobalConstants.Inf;

import jline.lang.JobClass;

import java.io.Serializable;
import java.util.List;

/**
 * A service section with an infinite number of servers (pure delay).
 * <p>
 * This server models a pure delay station where every arriving job immediately
 * enters service without queueing. Also known as an M/G/∞ queue, it represents
 * scenarios where resources are unlimited, such as think time in interactive systems
 * or self-service operations where capacity is effectively unlimited.
 */
public class InfiniteServer extends Server implements Serializable {
    /**
     * Creates a new infinite server section for the specified job classes.
     * Sets the number of servers to infinity, ensuring no queueing delays.
     * 
     * @param jobClasses the list of job classes this infinite server will handle
     */
    public InfiniteServer(List<JobClass> jobClasses) {
        super(jobClasses);
        this.numberOfServers = Inf;
        this.className = "InfiniteServer";
    }
}
