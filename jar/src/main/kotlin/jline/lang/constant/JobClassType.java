/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.constant;

/**
 * Enumeration of job class types in queueing network models.
 * 
 * <p>Job classes represent different categories of customers or workloads in a queueing
 * network. Each class can have distinct service requirements, routing behavior, and
 * performance characteristics. The type of a job class determines its fundamental
 * behavior in the network.</p>
 * 
 * <p>The three types are:
 * <ul>
 *   <li><b>OPEN:</b> Jobs arrive from external sources and eventually leave the system</li>
 *   <li><b>CLOSED:</b> Fixed population of jobs that circulate indefinitely in the network</li>
 *   <li><b>DISABLED:</b> Job class is inactive and not processed by the network</li>
 * </ul>
 * </p>
 * 
 * <p>Networks can contain a mix of open and closed classes, forming a mixed queueing
 * network. Open classes model external arrivals (e.g., web requests), while closed
 * classes model finite resources (e.g., thread pools).</p>
 * 
 * @see jline.lang.JobClass
 * @see jline.lang.OpenClass
 * @see jline.lang.ClosedClass
 * @since 1.0
 */
public enum JobClassType {
    /** Open class - jobs arrive from external sources and depart the system */
    OPEN,
    
    /** Closed class - fixed population of jobs circulating in the network */
    CLOSED,
    
    /** Disabled class - inactive job class not processed by the system */
    DISABLED
}
