/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.jmt;

import jline.lang.constant.SolverType;
import jline.solvers.SolverOptions;

/**
 * Configuration options for Java Modelling Tools (JMT) solver integration.
 * 
 * <p>JMTOptions extends the base solver options to provide configuration
 * for interfacing with the JMT simulation engine. JMT provides discrete-event
 * simulation capabilities for complex queueing networks that may not be
 * analytically tractable.</p>
 * 
 * <p>Key JMT integration features:
 * <ul>
 *   <li>Discrete-event simulation via JMT engine</li>
 *   <li>Supports complex network topologies</li>
 *   <li>Handles general service and interarrival time distributions</li>
 *   <li>Fork-join networks and advanced scheduling policies</li>
 *   <li>Statistical analysis with confidence intervals</li>
 * </ul>
 * </p>
 * 
 * <p><strong>Note:</strong> Requires JMT.jar to be available in the classpath
 * for simulation functionality.</p>
 * 
 * @see SolverJMT
 * @see SolverOptions
 * @since 1.0
 */
public class JMTOptions extends SolverOptions {
    
    /**
     * Constructs a new JMTOptions instance with default JMT solver configuration.
     */
    public JMTOptions() {
        super(SolverType.JMT);
    }
}
