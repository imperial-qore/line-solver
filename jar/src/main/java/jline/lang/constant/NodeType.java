/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.constant;

/**
 * Enumeration of node types available in queueing network models.
 * 
 * <p>Each node type represents a different functional component in a queueing network,
 * with specific behavior for processing, routing, or transforming jobs. Nodes are the
 * fundamental building blocks that are connected to form complete network models.</p>
 * 
 * <p>Node categories include:
 * <ul>
 *   <li><b>Service nodes:</b> Queue, Delay - process jobs with service times</li>
 *   <li><b>Routing nodes:</b> Router, Fork, Join - control job flow</li>
 *   <li><b>Source/Sink nodes:</b> Source, Sink - job creation and destruction</li>
 *   <li><b>Special nodes:</b> Cache, ClassSwitch, Logger - specialized functionality</li>
 *   <li><b>Petri net nodes:</b> Place, Transition - for Petri net models</li>
 * </ul>
 * </p>
 * 
 * @see jline.lang.nodes.Node
 * @since 1.0
 */
public enum NodeType {
    /** Transition node in Petri net models - fires when enabled */
    Transition,
    
    /** Place node in Petri net models - holds tokens */
    Place,
    
    /** Fork node - splits jobs into parallel tasks */
    Fork,
    
    /** Router node - directs jobs without service time */
    Router,
    
    /** Cache node - implements cache replacement with hit/miss routing */
    Cache,
    
    /** Logger node - records job passage for monitoring */
    Logger,
    
    /** Class switch node - changes job class based on probabilities */
    ClassSwitch,
    
    /** Delay station - infinite servers with no queueing */
    Delay,
    
    /** Source node - generates jobs according to arrival process */
    Source,
    
    /** Sink node - removes completed jobs from the network */
    Sink,
    
    /** Join node - synchronizes parallel tasks from fork */
    Join,
    
    /** Queue station - finite servers with queueing */
    Queue,

    /** Finite Capacity Region - virtual node representing a blocking region */
    Region
}
