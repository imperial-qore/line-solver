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
public abstract class Section extends NetworkElement implements Serializable, Cloneable {
    String className;

    public Section(String className) {
        super("Section");
        this.className = className;
    }

    public String getClassName() {
        return className;
    }

    /**
     * Returns a copy of this section, twin of the MATLAB {@code copyElement} of
     * a Copyable section.
     *
     * The copy is shallow in the same sense as MATLAB's: the section object and
     * its containers are new, while the model elements they refer to (job
     * classes, destination nodes, distributions) stay shared. Subclasses holding
     * containers override this to copy them.
     *
     * @return a copy of this section
     */
    public Section copyElement() {
        try {
            return (Section) super.clone();
        } catch (CloneNotSupportedException e) {
            throw new RuntimeException("Failed to copy section " + this.getClassName(), e);
        }
    }
}
