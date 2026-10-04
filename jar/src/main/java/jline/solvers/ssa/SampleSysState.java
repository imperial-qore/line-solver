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
 * Container for system-wide state sampling results from SSA solver
 */
public class SampleSysState implements Serializable {
    
    /** List of stateful nodes that were sampled */
    public List<Node> handle;
    
    /** Time points of the sampling */
    public Matrix t;
    
    /** State information for each stateful node at each time point */
    public List<Matrix> state;
    
    /** List of events that occurred during sampling */
    public List<Event> event;
    
    /** Whether this represents aggregate data */
    public boolean isaggregate;
    
    public SampleSysState() {
        this.isaggregate = false;
    }
    
    public SampleSysState(List<Node> handle, Matrix t, List<Matrix> state, List<Event> event, boolean isaggregate) {
        this.handle = handle;
        this.t = t;
        this.state = state;
        this.event = event;
        this.isaggregate = isaggregate;
    }
}