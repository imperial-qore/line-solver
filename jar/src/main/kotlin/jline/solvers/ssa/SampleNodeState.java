/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.ssa;

import jline.lang.Event;
import jline.lang.nodes.Node;
import jline.util.matrix.Matrix;

import java.io.Serializable;
import java.util.List;

/**
 * Container for node state sampling results from SSA solver
 */
public class SampleNodeState implements Serializable {
    
    /** The node that was sampled */
    public Node handle;
    
    /** Time points of the sampling */
    public Matrix t;
    
    /** State information at each time point */
    public Matrix state;
    
    /** List of events that occurred during sampling */
    public List<Event> event;
    
    /** Whether this represents aggregate data */
    public boolean isaggregate;
    
    public SampleNodeState() {
        this.isaggregate = false;
    }
    
    public SampleNodeState(Node handle, Matrix t, Matrix state, List<Event> event, boolean isaggregate) {
        this.handle = handle;
        this.t = t;
        this.state = state;
        this.event = event;
        this.isaggregate = isaggregate;
    }
}