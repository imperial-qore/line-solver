/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import jline.util.NamedParam;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

/**
 * An abstract class for stochastic processes
 */
public abstract class Process implements Serializable {
    protected String name;
    protected List<NamedParam> params;

    /**
     * Creates a new process with the specified name and number of parameters.
     *
     * @param name     the name of this process type
     * @param numParam the number of parameters required
     */
    public Process(String name, int numParam) {
        this.name = name;
        this.params = new ArrayList<NamedParam>(numParam);
        for (int i = 0; i < numParam; i++) {
            this.params.add(new NamedParam("", -1));
        }
    }

    /**
     * Sample a value from the process
     * @return sampled value
     */
    public abstract Object sample();

    /**
     * Get the name of this process
     * @return process name
     */
    public String getName() {
        return name;
    }

    /**
     * Get the number of parameters
     * @return number of parameters
     */
    public int getNumParams() {
        return params.size();
    }

    /**
     * Set a parameter by id, name and value
     * @param id parameter index (1-based)
     * @param paramName parameter name
     * @param paramValue parameter value
     */
    public void setParam(int id, String paramName, Object paramValue) {
        if (id > params.size()) {
            // Expand params list if necessary
            while (params.size() < id) {
                params.add(new NamedParam("", -1));
            }
        }
        params.set(id - 1, new NamedParam(paramName, paramValue));
    }

    /**
     * Get a parameter by id
     * @param id parameter index (1-based)
     * @return the named parameter
     */
    public NamedParam getParam(int id) {
        if (id <= 0 || id > params.size()) {
            throw new IndexOutOfBoundsException("Parameter id " + id + " is out of bounds");
        }
        return params.get(id - 1);
    }

    // =================== KOTLIN-STYLE PROPERTY ALIASES ===================
    
    /**
     * Kotlin-style property alias for getName()
     */
    public String name() {
        return getName();
    }
    
    /**
     * Kotlin-style property alias for getNumParams()
     */
    public int numParams() {
        return getNumParams();
    }
    
    /**
     * Kotlin-style property alias for getParam(int id)
     */
    public NamedParam param(int id) {
        return getParam(id);
    }
}