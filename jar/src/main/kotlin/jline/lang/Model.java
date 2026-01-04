/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang;

import jline.GlobalConstants;

import java.io.Serializable;
import java.util.Locale;

/**
 * Class representing a model supported by the library
 */
public class Model implements Serializable {
    private String network_name;
    private String lineVersion;

    /**
     * Creates a new model with the specified name.
     * Sets the locale to US and initializes the LINE version from GlobalConstants.
     * 
     * @param modelName the name to assign to this model
     */
    public Model(String modelName) {
        Locale.setDefault(Locale.US);
        this.network_name = modelName;
        this.setVersion(GlobalConstants.Version);
    }

    /**
     * Gets the name of this model.
     * 
     * @return the model name
     */
    public String getName() {
        return this.network_name;
    }

    /**
     * Sets the name of this model.
     * 
     * @param setName the new name for this model
     */
    public void setName(String setName) {
        this.network_name = setName;
    }

    /**
     * Gets the LINE version associated with this model.
     * 
     * @return the LINE version string
     */
    public String getVersion() {
        return this.lineVersion;
    }

    /**
     * Sets the LINE version for this model.
     * 
     * @param version the LINE version string
     */
    public void setVersion(String version) {
        this.lineVersion = version;
    }

    /*
     * ===================================================================================
     * MISSING METHODS FROM MATLAB MODEL IMPLEMENTATION - NOT YET MIGRATED
     * ===================================================================================
     *
     * Based on analysis of /matlab/src/lang/Model.m
     */

    // =================== MISSING PROPERTIES ===================
    // private ModelAttribute attribute;  // Model attribute object for additional metadata

    // =================== ENHANCED INITIALIZATION ===================
    // Enhanced constructor that includes:
    // - Automatic LINE initialization check
    // - Version detection from GlobalConstants
    // - Splash screen management

    // =================== COPYABLE INTERFACE SUPPORT ===================
    // public Model copy()          // Create deep copy of model
    // public Model copyElement()   // Protected method for deep copying

}