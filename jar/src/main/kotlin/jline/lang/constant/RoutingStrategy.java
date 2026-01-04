/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.constant;

/**
 * Enumeration of routing strategies that determine how jobs are dispatched to downstream stations.
 * 
 * <p>Routing strategies control the decision-making process when a job completes service
 * at a station and multiple downstream stations are available. These strategies implement
 * different load balancing, fairness, and optimization policies.</p>
 * 
 * <p>Common strategies include:
 * <ul>
 *   <li>RAND - Random selection with equal probability</li>
 *   <li>PROB - Probabilistic routing with specified probabilities</li>
 *   <li>RROBIN - Round-robin cycling through destinations</li>
 *   <li>JSQ - Join the Shortest Queue for load balancing</li>
 * </ul>
 * </p>
 * 
 * @see jline.lang.nodes.Node#setRouting
 * @since 1.0
 */
public enum RoutingStrategy {
    /** Random routing - destination selected uniformly at random */
    RAND,
    
    /** Probabilistic routing - destination selected according to specified probabilities */
    PROB,
    
    /** Round Robin - destinations visited in cyclic order */
    RROBIN,
    
    /** Weighted Round Robin - cyclic routing with weighted visits */
    WRROBIN,
    
    /** Join Shortest Queue - route to destination with fewest waiting jobs */
    JSQ,
    
    /** Firing routing - used for transitions in Petri net models */
    FIRING,
    
    /** Power of K choices - select best among K randomly chosen destinations */
    KCHOICES,
    
    /** Reinforcement Learning - adaptive routing based on learned policies */
    RL,
    
    /** Disabled routing - no routing allowed (jobs are dropped) */
    DISABLED;

    /**
     * Converts a RoutingStrategy enum value to a feature string for analysis.
     * 
     * @param routing the routing strategy to convert
     * @return the feature string representation, or empty string if not mapped
     */
    public static String toFeature(RoutingStrategy routing) {
        switch (routing) {
            case RAND:
                return "RoutingStrategy_RAND";
            case PROB:
                return "RoutingStrategy_PROB";
            case RROBIN:
                return "RoutingStrategy_RROBIN";
            case WRROBIN:
                return "RoutingStrategy_WRROBIN";
            case KCHOICES:
                return "RoutingStrategy_KCHOICES";
            default:
                return "";
        }
    }
}
