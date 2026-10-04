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
public class Model implements Copyable {
    private String network_name;
    private String lineVersion;
    protected ModelAttribute attribute;

    /**
     * Creates a new model with the specified name.
     * Sets the locale to US and initializes the LINE version from GlobalConstants.
     *
     * Mirrors the constructor of {@code matlab/src/lang/Model.m}: it first makes
     * sure LINE is initialized (the MATLAB {@code lineStart} guard on an empty
     * {@code GlobalConstants.Verbose}, here the singleton), then records the
     * trimmed version and the name.
     *
     * @param modelName the name to assign to this model
     */
    public Model(String modelName) {
        Locale.setDefault(Locale.US);
        GlobalConstants.getInstance();
        this.attribute = new ModelAttribute();
        this.setVersion(GlobalConstants.Version.trim());
        this.setName(modelName);
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

    /**
     * Gets the metadata container of this model, twin of the MATLAB
     * {@code Model.attribute} property.
     *
     * @return the model attribute container
     */
    public ModelAttribute getAttribute() {
        return this.attribute;
    }

    /**
     * Sets the metadata container of this model.
     *
     * @param attribute the container to install
     */
    public void setAttribute(ModelAttribute attribute) {
        this.attribute = attribute;
    }

}