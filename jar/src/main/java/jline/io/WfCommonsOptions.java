/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.io;

/**
 * Options for loading WfCommons workflow files.
 */
public class WfCommonsOptions {

    /**
     * Distribution type to use for task service times.
     */
    public enum DistributionType {
        /** Exponential distribution (default) */
        EXP,
        /** Deterministic (fixed) service time */
        DET,
        /** Acyclic Phase-type distribution */
        APH,
        /** Hyper-exponential distribution */
        HYPEREXP
    }

    private DistributionType distributionType = DistributionType.EXP;
    private double defaultSCV = 1.0;
    private double defaultRuntime = 1.0;
    private boolean useExecutionData = true;
    private boolean storeMetadata = true;

    /**
     * Create options with default values.
     */
    public WfCommonsOptions() {
    }

    /**
     * Get the distribution type.
     * @return Distribution type
     */
    public DistributionType getDistributionType() {
        return distributionType;
    }

    /**
     * Set the distribution type for task service times.
     * @param distributionType Distribution type
     * @return this for chaining
     */
    public WfCommonsOptions setDistributionType(DistributionType distributionType) {
        this.distributionType = distributionType;
        return this;
    }

    /**
     * Get the default SCV for APH/HyperExp distributions.
     * @return Default SCV
     */
    public double getDefaultSCV() {
        return defaultSCV;
    }

    /**
     * Set the default SCV for APH/HyperExp distributions.
     * @param defaultSCV Default SCV value
     * @return this for chaining
     */
    public WfCommonsOptions setDefaultSCV(double defaultSCV) {
        this.defaultSCV = defaultSCV;
        return this;
    }

    /**
     * Get the default runtime when execution data is missing.
     * @return Default runtime in seconds
     */
    public double getDefaultRuntime() {
        return defaultRuntime;
    }

    /**
     * Set the default runtime when execution data is missing.
     * @param defaultRuntime Default runtime in seconds
     * @return this for chaining
     */
    public WfCommonsOptions setDefaultRuntime(double defaultRuntime) {
        this.defaultRuntime = defaultRuntime;
        return this;
    }

    /**
     * Check if execution data should be used.
     * @return true if execution data is used
     */
    public boolean isUseExecutionData() {
        return useExecutionData;
    }

    /**
     * Set whether to use execution data if available.
     * @param useExecutionData true to use execution data
     * @return this for chaining
     */
    public WfCommonsOptions setUseExecutionData(boolean useExecutionData) {
        this.useExecutionData = useExecutionData;
        return this;
    }

    /**
     * Check if metadata should be stored.
     * @return true if metadata is stored
     */
    public boolean isStoreMetadata() {
        return storeMetadata;
    }

    /**
     * Set whether to store WfCommons metadata in activities.
     * @param storeMetadata true to store metadata
     * @return this for chaining
     */
    public WfCommonsOptions setStoreMetadata(boolean storeMetadata) {
        this.storeMetadata = storeMetadata;
        return this;
    }

    /**
     * Create options with exponential distribution.
     * @return Options configured for exponential distribution
     */
    public static WfCommonsOptions exponential() {
        return new WfCommonsOptions().setDistributionType(DistributionType.EXP);
    }

    /**
     * Create options with deterministic service times.
     * @return Options configured for deterministic distribution
     */
    public static WfCommonsOptions deterministic() {
        return new WfCommonsOptions().setDistributionType(DistributionType.DET);
    }
}
