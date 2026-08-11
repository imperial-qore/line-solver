/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang;

import jline.lang.nodes.Transition;

import java.io.Serializable;

/**
 * Superclass representing a class of jobs
 */
public class Mode extends NetworkElement implements Serializable {

    Transition transition;
    int index;

    public Mode(Transition parent, String name) {
        super(name);
        this.transition = parent;
        this.index = parent.getModes().size() + 1;
    }

    public int getIndex() {
        return index;
    }

    public Transition getTransition() {
        return transition;
    }

    public void printSummary() {
        System.out.format("Mode: %s\n", this.getName());
        System.out.println();
    }
}
