/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.rest.routes;

import jline.rest.jobs.JobManager;
import jline.rest.util.JsonTransformer;

import java.lang.management.ManagementFactory;
import java.lang.management.MemoryMXBean;
import java.lang.management.MemoryUsage;
import java.util.concurrent.atomic.AtomicLong;

import static spark.Spark.*;

/**
 * Prometheus-compatible metrics endpoint for SRE/DevOps integration.
 */
public class MetricsRoutes {

    private static final AtomicLong totalRequests = new AtomicLong(0);
    private static final AtomicLong totalSolveRequests = new AtomicLong(0);
    private static final AtomicLong totalAnalysisRequests = new AtomicLong(0);
    private static final AtomicLong failedRequests = new AtomicLong(0);
    private static final AtomicLong totalSolveTimeMs = new AtomicLong(0);

    private static JobManager jobManager;

    /**
     * Register metrics routes.
     *
     * @param basePath The API base path (e.g., "/api/v1")
     * @param json     The JSON transformer
     * @param manager  The job manager for job metrics
     */
    public static void register(String basePath, JsonTransformer json, JobManager manager) {
        jobManager = manager;

        // Prometheus metrics endpoint (text format)
        get(basePath + "/metrics", (req, res) -> {
            res.type("text/plain; version=0.0.4; charset=utf-8");
            return generatePrometheusMetrics();
        });

        // JSON metrics endpoint (for non-Prometheus consumers)
        get(basePath + "/metrics/json", (req, res) -> {
            return generateJsonMetrics();
        }, json);
    }

    /**
     * Generate Prometheus text format metrics.
     */
    private static String generatePrometheusMetrics() {
        StringBuilder sb = new StringBuilder();
        long timestamp = System.currentTimeMillis();

        // Request counters
        sb.append("# HELP line_requests_total Total number of API requests\n");
        sb.append("# TYPE line_requests_total counter\n");
        sb.append("line_requests_total ").append(totalRequests.get()).append("\n\n");

        sb.append("# HELP line_solve_requests_total Total number of solve requests\n");
        sb.append("# TYPE line_solve_requests_total counter\n");
        sb.append("line_solve_requests_total ").append(totalSolveRequests.get()).append("\n\n");

        sb.append("# HELP line_analysis_requests_total Total number of analysis requests\n");
        sb.append("# TYPE line_analysis_requests_total counter\n");
        sb.append("line_analysis_requests_total ").append(totalAnalysisRequests.get()).append("\n\n");

        sb.append("# HELP line_failed_requests_total Total number of failed requests\n");
        sb.append("# TYPE line_failed_requests_total counter\n");
        sb.append("line_failed_requests_total ").append(failedRequests.get()).append("\n\n");

        // Solve time
        sb.append("# HELP line_solve_time_seconds_total Total time spent solving models\n");
        sb.append("# TYPE line_solve_time_seconds_total counter\n");
        sb.append("line_solve_time_seconds_total ").append(totalSolveTimeMs.get() / 1000.0).append("\n\n");

        // Average solve time
        long solveCount = totalSolveRequests.get();
        double avgSolveTime = solveCount > 0 ? (totalSolveTimeMs.get() / 1000.0) / solveCount : 0;
        sb.append("# HELP line_solve_time_seconds_avg Average solve time in seconds\n");
        sb.append("# TYPE line_solve_time_seconds_avg gauge\n");
        sb.append("line_solve_time_seconds_avg ").append(avgSolveTime).append("\n\n");

        // Job metrics
        if (jobManager != null) {
            sb.append("# HELP line_jobs_active Number of currently running jobs\n");
            sb.append("# TYPE line_jobs_active gauge\n");
            sb.append("line_jobs_active ").append(jobManager.getActiveJobCount()).append("\n\n");

            sb.append("# HELP line_jobs_total Total number of jobs (including completed)\n");
            sb.append("# TYPE line_jobs_total gauge\n");
            sb.append("line_jobs_total ").append(jobManager.getTotalJobCount()).append("\n\n");
        }

        // JVM metrics
        MemoryMXBean memoryBean = ManagementFactory.getMemoryMXBean();
        MemoryUsage heapUsage = memoryBean.getHeapMemoryUsage();
        MemoryUsage nonHeapUsage = memoryBean.getNonHeapMemoryUsage();

        sb.append("# HELP jvm_memory_heap_bytes_used JVM heap memory used\n");
        sb.append("# TYPE jvm_memory_heap_bytes_used gauge\n");
        sb.append("jvm_memory_heap_bytes_used ").append(heapUsage.getUsed()).append("\n\n");

        sb.append("# HELP jvm_memory_heap_bytes_max JVM heap memory max\n");
        sb.append("# TYPE jvm_memory_heap_bytes_max gauge\n");
        sb.append("jvm_memory_heap_bytes_max ").append(heapUsage.getMax()).append("\n\n");

        sb.append("# HELP jvm_memory_nonheap_bytes_used JVM non-heap memory used\n");
        sb.append("# TYPE jvm_memory_nonheap_bytes_used gauge\n");
        sb.append("jvm_memory_nonheap_bytes_used ").append(nonHeapUsage.getUsed()).append("\n\n");

        // Thread count
        sb.append("# HELP jvm_threads_current Current number of JVM threads\n");
        sb.append("# TYPE jvm_threads_current gauge\n");
        sb.append("jvm_threads_current ").append(Thread.activeCount()).append("\n\n");

        // Uptime
        long uptimeMs = ManagementFactory.getRuntimeMXBean().getUptime();
        sb.append("# HELP line_uptime_seconds Server uptime in seconds\n");
        sb.append("# TYPE line_uptime_seconds counter\n");
        sb.append("line_uptime_seconds ").append(uptimeMs / 1000.0).append("\n\n");

        // Available processors
        sb.append("# HELP line_processors_available Number of available processors\n");
        sb.append("# TYPE line_processors_available gauge\n");
        sb.append("line_processors_available ").append(Runtime.getRuntime().availableProcessors()).append("\n");

        return sb.toString();
    }

    /**
     * Generate JSON format metrics.
     */
    private static java.util.Map<String, Object> generateJsonMetrics() {
        java.util.Map<String, Object> metrics = new java.util.LinkedHashMap<String, Object>();

        // Request metrics
        java.util.Map<String, Object> requests = new java.util.LinkedHashMap<String, Object>();
        requests.put("total", totalRequests.get());
        requests.put("solve", totalSolveRequests.get());
        requests.put("analysis", totalAnalysisRequests.get());
        requests.put("failed", failedRequests.get());
        metrics.put("requests", requests);

        // Timing metrics
        java.util.Map<String, Object> timing = new java.util.LinkedHashMap<String, Object>();
        timing.put("totalSolveTimeSeconds", totalSolveTimeMs.get() / 1000.0);
        long solveCount = totalSolveRequests.get();
        timing.put("avgSolveTimeSeconds", solveCount > 0 ? (totalSolveTimeMs.get() / 1000.0) / solveCount : 0);
        metrics.put("timing", timing);

        // Job metrics
        if (jobManager != null) {
            java.util.Map<String, Object> jobs = new java.util.LinkedHashMap<String, Object>();
            jobs.put("active", jobManager.getActiveJobCount());
            jobs.put("total", jobManager.getTotalJobCount());
            metrics.put("jobs", jobs);
        }

        // JVM metrics
        MemoryMXBean memoryBean = ManagementFactory.getMemoryMXBean();
        MemoryUsage heapUsage = memoryBean.getHeapMemoryUsage();

        java.util.Map<String, Object> jvm = new java.util.LinkedHashMap<String, Object>();
        jvm.put("heapUsedBytes", heapUsage.getUsed());
        jvm.put("heapMaxBytes", heapUsage.getMax());
        jvm.put("heapUtilization", heapUsage.getMax() > 0 ? (double) heapUsage.getUsed() / heapUsage.getMax() : 0);
        jvm.put("threads", Thread.activeCount());
        jvm.put("processors", Runtime.getRuntime().availableProcessors());
        jvm.put("uptimeSeconds", ManagementFactory.getRuntimeMXBean().getUptime() / 1000.0);
        metrics.put("jvm", jvm);

        return metrics;
    }

    // Methods for incrementing counters (called by other routes)

    public static void incrementRequests() {
        totalRequests.incrementAndGet();
    }

    public static void incrementSolveRequests() {
        totalSolveRequests.incrementAndGet();
    }

    public static void incrementAnalysisRequests() {
        totalAnalysisRequests.incrementAndGet();
    }

    public static void incrementFailedRequests() {
        failedRequests.incrementAndGet();
    }

    public static void addSolveTime(long milliseconds) {
        totalSolveTimeMs.addAndGet(milliseconds);
    }

    /**
     * Reset all counters (for testing).
     */
    public static void resetCounters() {
        totalRequests.set(0);
        totalSolveRequests.set(0);
        totalAnalysisRequests.set(0);
        failedRequests.set(0);
        totalSolveTimeMs.set(0);
    }
}
