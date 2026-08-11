/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.util;

import java.io.Serializable;

/**
 * A container for storing named parameters with string identifiers and object values.
 * 
 * <p>This utility class provides a simple key-value pair specifically designed
 * for parameter storage in configuration contexts. It's commonly used in
 * distribution parameters, solver options, and other configurable components.</p>
 * 
 * @see jline.lang.processes.Distribution
 * @since 1.0
 */
public class NamedParam implements Serializable {
    /** The parameter name/identifier */
    protected String name;
    
    /** The parameter value */
    protected Object value;

    /**
     * Constructs a new named parameter with the specified name and value.
     * 
     * @param name  the parameter name/identifier
     * @param value the parameter value
     */
    public NamedParam(String name, Object value) {
        this.name = name;
        this.value = value;
    }

    /**
     * Returns the value of this named parameter.
     * 
     * @return the parameter value
     */
    public Object getValue() {
        return this.value;
    }
}
