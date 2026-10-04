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
    private boolean resolveChildrenByName = true;

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
     * Check if child and parent references fall back to task names when no task id matches.
     * @return true if children are resolved by id, then by unique task name
     */
    public boolean isResolveChildrenByName() {
        return resolveChildrenByName;
    }

    /**
     * Set whether a child or parent reference that matches no task id is looked up among the task names.
     *
     * <p>The Pegasus traces of schema 1.4 (e.g. Montage) give each task an id such as
     * {@code ID0000001} but list children and parents by name ({@code mBackground_ID0000013}), so an
     * id-only lookup drops every edge. A name shared by several tasks is ambiguous and does not
     * resolve. A reference that resolves to no task is reported by a warning and dropped.
     * Default true; false restores the id-only lookup.</p>
     *
     * @param resolveChildrenByName true to resolve by id, then by unique task name
     * @return this for chaining
     */
    public WfCommonsOptions setResolveChildrenByName(boolean resolveChildrenByName) {
        this.resolveChildrenByName = resolveChildrenByName;
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
