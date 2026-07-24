/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.rest.routes;

import jline.rest.jobs.Job;
import jline.rest.jobs.JobManager;
import jline.rest.jobs.JobStatus;
import jline.rest.model.ErrorResponse;
import jline.rest.model.SolveRequest;
import jline.rest.util.JsonTransformer;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static spark.Spark.*;

/**
 * Job management endpoints for asynchronous solving.
 * Provides endpoints for submitting, monitoring, and cancelling jobs.
 */
public class JobRoutes {

    private static JobManager jobManager;

    /**
     * Register job routes.
     *
     * @param basePath The API base path (e.g., "/api/v1")
     * @param json     The JSON transformer
     * @param manager  The job manager instance
     */
    public static void register(String basePath, JsonTransformer json, JobManager manager) {
        jobManager = manager;

        // Submit async job
        post(basePath + "/models/solve/async", (req, res) -> {
            try {
                SolveRequest request = json.fromJson(req.body(), SolveRequest.class);
                if (request == null) {
                    res.status(400);
                    return new ErrorResponse("INVALID_REQUEST", "Request body is required");
                }

                // Validate request
                String validationError = request.validate();
                if (validationError != null) {
                    res.status(400);
                    return new ErrorResponse("VALIDATION_ERROR", validationError);
                }

                // Submit job
                Job job = jobManager.submitJob(request);

                res.status(202); // Accepted
                return createJobResponse(job);
            } catch (Exception e) {
                res.status(500);
                return new ErrorResponse("SUBMIT_ERROR", "Failed to submit job: " + e.getMessage());
            }
        }, json);

        // Get job status
        get(basePath + "/jobs/:jobId", (req, res) -> {
            String jobId = req.params(":jobId");
            Job job = jobManager.getJob(jobId);

            if (job == null) {
                res.status(404);
                return new ErrorResponse("NOT_FOUND", "Job not found: " + jobId);
            }

            return createJobResponse(job);
        }, json);

        // List all jobs
        get(basePath + "/jobs", (req, res) -> {
            String statusFilter = req.queryParams("status");

            List<Map<String, Object>> jobList = new ArrayList<Map<String, Object>>();

            for (Job job : jobManager.getAllJobs()) {
                if (statusFilter == null || job.getStatus().getValue().equals(statusFilter)) {
                    jobList.add(createJobSummary(job));
                }
            }

            Map<String, Object> response = new HashMap<String, Object>();
            response.put("jobs", jobList);
            response.put("count", jobList.size());
            response.put("activeCount", jobManager.getActiveJobCount());

            return response;
        }, json);

        // Cancel job
        delete(basePath + "/jobs/:jobId", (req, res) -> {
            String jobId = req.params(":jobId");
            Job job = jobManager.getJob(jobId);

            if (job == null) {
                res.status(404);
                return new ErrorResponse("NOT_FOUND", "Job not found: " + jobId);
            }

            if (job.isTerminal()) {
                res.status(400);
                return new ErrorResponse("ALREADY_TERMINAL",
                    "Job is already in terminal state: " + job.getStatus().getValue());
            }

            boolean cancelled = jobManager.cancelJob(jobId);
            if (cancelled) {
                Map<String, Object> response = new HashMap<String, Object>();
                response.put("jobId", jobId);
                response.put("status", "cancelled");
                response.put("message", "Job cancelled successfully");
                return response;
            } else {
                res.status(500);
                return new ErrorResponse("CANCEL_FAILED", "Failed to cancel job");
            }
        }, json);

        // Stream job progress (Server-Sent Events)
        get(basePath + "/jobs/:jobId/stream", (req, res) -> {
            String jobId = req.params(":jobId");
            Job job = jobManager.getJob(jobId);

            if (job == null) {
                res.status(404);
                res.type("application/json");
                return json.toJson(new ErrorResponse("NOT_FOUND", "Job not found: " + jobId));
            }

            // Set SSE headers
            res.type("text/event-stream");
            res.header("Cache-Control", "no-cache");
            res.header("Connection", "keep-alive");
            res.header("X-Accel-Buffering", "no");

            try {
                PrintWriter writer = res.raw().getWriter();

                // Send initial status
                sendSSEEvent(writer, "status", json.toJson(createJobResponse(job)));

                // Poll for updates until job is terminal or client disconnects
                JobStatus lastStatus = job.getStatus();
                double lastProgress = job.getProgress();

                while (!job.isTerminal()) {
                    Thread.sleep(500); // Poll every 500ms

                    // Check if client disconnected
                    if (res.raw().getOutputStream() == null) {
                        break;
                    }

                    // Send update if status or progress changed
                    if (job.getStatus() != lastStatus || job.getProgress() != lastProgress) {
                        lastStatus = job.getStatus();
                        lastProgress = job.getProgress();

                        Map<String, Object> update = new HashMap<String, Object>();
                        update.put("jobId", job.getId());
                        update.put("status", job.getStatus().getValue());
                        update.put("progress", job.getProgress());
                        update.put("message", job.getProgressMessage());

                        sendSSEEvent(writer, "progress", json.toJson(update));
                    }

                    // Send heartbeat
                    writer.write(": heartbeat\n\n");
                    writer.flush();

                    if (writer.checkError()) {
                        break; // Client disconnected
                    }
                }

                // Send final result
                sendSSEEvent(writer, "complete", json.toJson(createJobResponse(job)));
                writer.flush();

            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
            } catch (IOException e) {
                // Client disconnected, ignore
            }

            return "";
        });
    }

    /**
     * Send an SSE event.
     */
    private static void sendSSEEvent(PrintWriter writer, String event, String data) {
        writer.write("event: " + event + "\n");
        writer.write("data: " + data + "\n\n");
        writer.flush();
    }

    /**
     * Create a full job response.
     */
    private static Map<String, Object> createJobResponse(Job job) {
        Map<String, Object> response = new HashMap<String, Object>();
        response.put("jobId", job.getId());
        response.put("status", job.getStatus().getValue());
        response.put("progress", job.getProgress());
        response.put("progressMessage", job.getProgressMessage());
        response.put("createdAt", job.getCreatedAt().toString());

        if (job.getStartedAt() != null) {
            response.put("startedAt", job.getStartedAt().toString());
        }
        if (job.getCompletedAt() != null) {
            response.put("completedAt", job.getCompletedAt().toString());
        }

        response.put("elapsedSeconds", job.getElapsedSeconds());

        if (job.getStatus() == JobStatus.COMPLETED && job.getResponse() != null) {
            response.put("results", job.getResponse().getResults());
            response.put("solver", job.getResponse().getSolver());
            response.put("runtime", job.getResponse().getRuntime());
        }

        if (job.getStatus() == JobStatus.FAILED && job.getError() != null) {
            response.put("error", job.getError());
        }

        // Include solver from request
        if (job.getRequest() != null) {
            response.put("solver", job.getRequest().getSolver());
        }

        return response;
    }

    /**
     * Create a job summary (for listing).
     */
    private static Map<String, Object> createJobSummary(Job job) {
        Map<String, Object> summary = new HashMap<String, Object>();
        summary.put("jobId", job.getId());
        summary.put("status", job.getStatus().getValue());
        summary.put("progress", job.getProgress());
        summary.put("createdAt", job.getCreatedAt().toString());
        summary.put("elapsedSeconds", job.getElapsedSeconds());

        if (job.getRequest() != null) {
            summary.put("solver", job.getRequest().getSolver());
        }

        return summary;
    }

    /**
     * Get the job manager instance.
     *
     * @return The job manager
     */
    public static JobManager getJobManager() {
        return jobManager;
    }
}
