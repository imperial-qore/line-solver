/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.rest.jobs;

/**
 * Status of an asynchronous job.
 */
public enum JobStatus {
    /**
     * Job has been created but not yet started.
     */
    PENDING("pending"),

    /**
     * Job is currently running.
     */
    RUNNING("running"),

    /**
     * Job completed successfully.
     */
    COMPLETED("completed"),

    /**
     * Job failed with an error.
     */
    FAILED("failed"),

    /**
     * Job was cancelled by user request.
     */
    CANCELLED("cancelled");

    private final String value;

    JobStatus(String value) {
        this.value = value;
    }

    public String getValue() {
        return value;
    }

    @Override
    public String toString() {
        return value;
    }
}
