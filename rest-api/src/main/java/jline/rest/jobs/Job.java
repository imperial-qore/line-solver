/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.rest.jobs;

import jline.rest.model.SolveRequest;
import jline.rest.model.SolveResponse;

import java.time.Instant;
import java.util.UUID;
import java.util.concurrent.Future;

/**
 * Represents an asynchronous solve job.
 */
public class Job {

    private final String id;
    private final SolveRequest request;
    private final Instant createdAt;
    private volatile JobStatus status;
    private volatile double progress;
    private volatile String progressMessage;
    private volatile Instant startedAt;
    private volatile Instant completedAt;
    private volatile SolveResponse response;
    private volatile String error;
    private volatile Future<?> future;

    /**
     * Create a new job with the given request.
     *
     * @param request The solve request
     */
    public Job(SolveRequest request) {
        this.id = UUID.randomUUID().toString();
        this.request = request;
        this.createdAt = Instant.now();
        this.status = JobStatus.PENDING;
        this.progress = 0.0;
        this.progressMessage = "Pending";
    }

    /**
     * Create a new job with a specific ID.
     *
     * @param id      The job ID
     * @param request The solve request
     */
    public Job(String id, SolveRequest request) {
        this.id = id;
        this.request = request;
        this.createdAt = Instant.now();
        this.status = JobStatus.PENDING;
        this.progress = 0.0;
        this.progressMessage = "Pending";
    }

    public String getId() {
        return id;
    }

    public SolveRequest getRequest() {
        return request;
    }

    public Instant getCreatedAt() {
        return createdAt;
    }

    public JobStatus getStatus() {
        return status;
    }

    public void setStatus(JobStatus status) {
        this.status = status;
        if (status == JobStatus.RUNNING && startedAt == null) {
            this.startedAt = Instant.now();
        } else if (status == JobStatus.COMPLETED || status == JobStatus.FAILED || status == JobStatus.CANCELLED) {
            this.completedAt = Instant.now();
        }
    }

    public double getProgress() {
        return progress;
    }

    public void setProgress(double progress) {
        this.progress = Math.max(0.0, Math.min(1.0, progress));
    }

    public String getProgressMessage() {
        return progressMessage;
    }

    public void setProgressMessage(String progressMessage) {
        this.progressMessage = progressMessage;
    }

    public Instant getStartedAt() {
        return startedAt;
    }

    public Instant getCompletedAt() {
        return completedAt;
    }

    public SolveResponse getResponse() {
        return response;
    }

    public void setResponse(SolveResponse response) {
        this.response = response;
    }

    public String getError() {
        return error;
    }

    public void setError(String error) {
        this.error = error;
    }

    public Future<?> getFuture() {
        return future;
    }

    public void setFuture(Future<?> future) {
        this.future = future;
    }

    /**
     * Check if the job is in a terminal state.
     *
     * @return true if completed, failed, or cancelled
     */
    public boolean isTerminal() {
        return status == JobStatus.COMPLETED ||
               status == JobStatus.FAILED ||
               status == JobStatus.CANCELLED;
    }

    /**
     * Get the elapsed time in seconds since job creation.
     *
     * @return Elapsed time in seconds
     */
    public double getElapsedSeconds() {
        Instant end = completedAt != null ? completedAt : Instant.now();
        Instant start = startedAt != null ? startedAt : createdAt;
        return (end.toEpochMilli() - start.toEpochMilli()) / 1000.0;
    }

    /**
     * Cancel this job if it's still running.
     *
     * @return true if cancellation was successful
     */
    public boolean cancel() {
        if (isTerminal()) {
            return false;
        }
        if (future != null && !future.isDone()) {
            future.cancel(true);
        }
        setStatus(JobStatus.CANCELLED);
        setProgressMessage("Cancelled by user");
        return true;
    }

    /**
     * Mark the job as started.
     */
    public void markStarted() {
        setStatus(JobStatus.RUNNING);
        setProgress(0.0);
        setProgressMessage("Running");
    }

    /**
     * Mark the job as completed with a response.
     *
     * @param response The solve response
     */
    public void markCompleted(SolveResponse response) {
        setResponse(response);
        setStatus(JobStatus.COMPLETED);
        setProgress(1.0);
        setProgressMessage("Completed");
    }

    /**
     * Mark the job as failed with an error.
     *
     * @param error The error message
     */
    public void markFailed(String error) {
        setError(error);
        setStatus(JobStatus.FAILED);
        setProgressMessage("Failed: " + error);
    }
}
