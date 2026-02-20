/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.nodes;

import jline.util.matrix.Matrix;

import java.io.Serializable;

/**
 * A node that can have a state
 */
public class StatefulNode extends Node implements Serializable {
    private Integer statefulIndex;
    private Matrix state;
    private Matrix space;
    private Matrix statePrior;

    public StatefulNode(String name) {
        super(name);
        statefulIndex = null;
        state = new Matrix(0, 0, 0);
        space = new Matrix(0, 0, 0);
        statePrior = new Matrix(0, 0, 0);
    }

    protected void clearState() {
        this.state.reshape(0, 0, 0);
    }

    public int getNumberOfServers() {
        return 1;
    }

    public Matrix getState() {
        return this.state.copy();
    }

    public void setState(int state) {
        setState(Matrix.singleton(state));
    }

    public void setState(Matrix state) {
        this.state = state.copy();
    }

    public Matrix getStatePrior() {
        return this.statePrior.copy();
    }

    public void setStatePrior(Matrix prior) {
        // Validate dimensions before modifying state
        if (space.getNumRows() != prior.getNumRows()) {
            throw new RuntimeException("The prior probability vector must have the same rows of the station state space vector");
        }
        this.statePrior = prior.copy();
        // Invalidate Network's cached struct so stateprior is re-read
        if (this.model != null) {
            this.model.setHasStruct(false);
        }
    }

    public Matrix getStateSpace() {
        return this.space.copy();
    }

    public void setStateSpace(Matrix space) {
        this.space = space.copy();
    }

    public int getStatefulIndex() {
        if (this.statefulIndex == null) {
            this.statefulIndex = this.model.getStatefulNodeIndex(this);
            if (this.statefulIndex == -1) {
                throw new RuntimeException("StatefulNode '" + this.getName() + "' could not find its index in the network's stateful nodes list. " +
                        "This may indicate the node was not properly added to the network or there is a reference mismatch.");
            }
        }
        return this.statefulIndex;
    }

    public void resetStateSpace() {
        this.space = new Matrix(0, 0, 0);
    }

}
