/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.sections;

import jline.lang.NetworkElement;

import java.io.Serializable;

/**
 * A general class modeling a node section
 */
public abstract class Section extends NetworkElement implements Serializable {
    String className;

    public Section(String className) {
        super("Section");
        this.className = className;
    }

    public String getClassName() {
        return className;
    }
}
