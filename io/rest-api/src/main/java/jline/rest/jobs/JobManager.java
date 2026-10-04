/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.rest.jobs;

import jline.rest.handlers.SolveHandler;
import jline.rest.model.SolveRequest;
import jline.rest.model.SolveResponse;

import java.time.Duration;
import java.time.Instant;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Manages asynchronous solve jobs.
 * <p>
 * This class maintains an in-memory store of jobs and provides methods
 * for submitting, querying, and cancelling jobs. Jobs are executed
 * asynchronously using a thread pool.
 */
public class JobManager {

    private static final int DEFAULT_THREAD_POOL_SIZE = Runtime.getRuntime().availableProcessors();
    private static final Duration DEFAULT_JOB_RETENTION = Duration.ofHours(1);
    private static final Duration DEFAULT_CLEANUP_INTERVAL = Duration.ofMinutes(5);

    private final ConcurrentHashMap<String, Job> jobs;
    private final ExecutorService executor;
    private final ScheduledExecutorService cleanupExecutor;
    private final SolveHandler solveHandler;
    private final Duration jobRetention;
    private final AtomicInteger activeJobCount;

    /**
     * Create a new JobManager with default settings.
     */
    public JobManager() {
        this(DEFAULT_THREAD_POOL_SIZE, DEFAULT_JOB_RETENTION);
    }

    /**
     * Create a new JobManager with specified settings.
     *
     * @param threadPoolSize Number of threads for job execution
     * @param jobRetention   How long to keep completed jobs
     */
    public JobManager(int threadPoolSize, Duration jobRetention) {
        this.jobs = new ConcurrentHashMap<String, Job>();
        this.executor = Executors.newFixedThreadPool(threadPoolSize, new ThreadFactory() {
            private final AtomicInteger counter = new AtomicInteger(0);

            @Override
            public Thread newThread(Runnable r) {
                Thread t = new Thread(r, "line-solver-" + counter.incrementAndGet());
                t.setDaemon(true);
                return t;
            }
        });
        this.cleanupExecutor = Executors.newSingleThreadScheduledExecutor(new ThreadFactory() {
            @Override
            public Thread newThread(Runnable r) {
                Thread t = new Thread(r, "line-job-cleanup");
                t.setDaemon(true);
                return t;
            }
        });
        this.solveHandler = new SolveHandler();
        this.jobRetention = jobRetention;
        this.activeJobCount = new AtomicInteger(0);

        // Schedule cleanup task
        cleanupExecutor.scheduleAtFixedRate(
            this::cleanupExpiredJobs,
            DEFAULT_CLEANUP_INTERVAL.toMillis(),
            DEFAULT_CLEANUP_INTERVAL.toMillis(),
            TimeUnit.MILLISECONDS
        );
    }

    /**
     * Submit a new job for asynchronous execution.
     *
     * @param request The solve request
     * @return The created job
     */
    public Job submitJob(SolveRequest request) {
        Job job = new Job(request);
        jobs.put(job.getId(), job);

        Future<?> future = executor.submit(() -> executeJob(job));
        job.setFuture(future);

        return job;
    }

    /**
     * Execute a job.
     */
    private void executeJob(Job job) {
        try {
            activeJobCount.incrementAndGet();
            job.markStarted();

            // Execute the solve
            SolveResponse response = solveHandler.solve(job.getRequest());

            if (job.getStatus() == JobStatus.CANCELLED) {
                // Job was cancelled while running
                return;
            }

            if ("completed".equals(response.getStatus())) {
                job.markCompleted(response);
            } else {
                job.markFailed(response.getError());
            }
        } catch (Exception e) {
            if (job.getStatus() != JobStatus.CANCELLED) {
                job.markFailed(e.getMessage());
            }
        } finally {
            activeJobCount.decrementAndGet();
        }
    }

    /**
     * Get a job by ID.
     *
     * @param jobId The job ID
     * @return The job, or null if not found
     */
    public Job getJob(String jobId) {
        return jobs.get(jobId);
    }

    /**
     * Cancel a job.
     *
     * @param jobId The job ID
     * @return true if the job was cancelled, false if not found or already terminal
     */
    public boolean cancelJob(String jobId) {
        Job job = jobs.get(jobId);
        if (job == null) {
            return false;
        }
        return job.cancel();
    }

    /**
     * Get all jobs.
     *
     * @return Collection of all jobs
     */
    public Collection<Job> getAllJobs() {
        return jobs.values();
    }

    /**
     * Get jobs by status.
     *
     * @param status The status to filter by
     * @return List of jobs with the given status
     */
    public List<Job> getJobsByStatus(JobStatus status) {
        List<Job> result = new ArrayList<Job>();
        for (Job job : jobs.values()) {
            if (job.getStatus() == status) {
                result.add(job);
            }
        }
        return result;
    }

    /**
     * Get the count of active (running) jobs.
     *
     * @return Number of running jobs
     */
    public int getActiveJobCount() {
        return activeJobCount.get();
    }

    /**
     * Get the total number of jobs (including completed).
     *
     * @return Total number of jobs
     */
    public int getTotalJobCount() {
        return jobs.size();
    }

    /**
     * Clean up expired jobs.
     */
    private void cleanupExpiredJobs() {
        Instant cutoff = Instant.now().minus(jobRetention);
        List<String> toRemove = new ArrayList<String>();

        for (Job job : jobs.values()) {
            if (job.isTerminal() && job.getCompletedAt() != null) {
                if (job.getCompletedAt().isBefore(cutoff)) {
                    toRemove.add(job.getId());
                }
            }
        }

        for (String jobId : toRemove) {
            jobs.remove(jobId);
        }
    }

    /**
     * Remove a specific job.
     *
     * @param jobId The job ID to remove
     * @return true if removed, false if not found
     */
    public boolean removeJob(String jobId) {
        Job job = jobs.get(jobId);
        if (job == null) {
            return false;
        }
        if (!job.isTerminal()) {
            job.cancel();
        }
        jobs.remove(jobId);
        return true;
    }

    /**
     * Shutdown the job manager.
     */
    public void shutdown() {
        // Cancel all running jobs
        for (Job job : jobs.values()) {
            if (!job.isTerminal()) {
                job.cancel();
            }
        }

        // Shutdown executors
        executor.shutdown();
        cleanupExecutor.shutdown();

        try {
            if (!executor.awaitTermination(5, TimeUnit.SECONDS)) {
                executor.shutdownNow();
            }
            if (!cleanupExecutor.awaitTermination(1, TimeUnit.SECONDS)) {
                cleanupExecutor.shutdownNow();
            }
        } catch (InterruptedException e) {
            executor.shutdownNow();
            cleanupExecutor.shutdownNow();
            Thread.currentThread().interrupt();
        }
    }
}
