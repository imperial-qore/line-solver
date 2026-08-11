/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang;

import java.io.Serializable;

/**
 * Superclass for model elements
 */
public class Element implements Copyable {
    protected String name;

    /**
     * Creates a new element with the specified name.
     * 
     * @param name the name to assign to this element
     */
    public Element(String name) {
        this.name = name;
    }

    /**
     * Gets the name of this element.
     * 
     * @return the element's name
     */
    public String getName() {
        return name;
    }

    /**
     * Sets the name of this element.
     * 
     * @param name the new name for this element
     */
    public void setName(String name) {
        this.name = name;
    }
}
