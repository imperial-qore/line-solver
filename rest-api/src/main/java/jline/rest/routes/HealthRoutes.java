/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.rest.routes;

import jline.lang.Model;
import jline.rest.util.JsonTransformer;

import java.lang.management.ManagementFactory;
import java.lang.management.MemoryMXBean;
import java.lang.management.RuntimeMXBean;
import java.util.HashMap;
import java.util.Map;

import static spark.Spark.get;

/**
 * Health and info endpoints for the REST API.
 * Provides Kubernetes-compatible health checks and server information.
 */
public class HealthRoutes {

    /**
     * Register health routes.
     *
     * @param basePath The API base path (e.g., "/api/v1")
     * @param json The JSON transformer
     */
    public static void register(String basePath, JsonTransformer json) {
        // Health check - simple liveness probe
        get(basePath + "/health", (req, res) -> {
            Map<String, Object> health = new HashMap<String, Object>();
            health.put("status", "healthy");
            health.put("timestamp", System.currentTimeMillis());
            return health;
        }, json);

        // Readiness check - more detailed
        get(basePath + "/ready", (req, res) -> {
            Map<String, Object> ready = new HashMap<String, Object>();
            ready.put("status", "ready");
            ready.put("timestamp", System.currentTimeMillis());

            // Check memory availability
            Runtime runtime = Runtime.getRuntime();
            long freeMemory = runtime.freeMemory();
            long maxMemory = runtime.maxMemory();
            double memoryUsage = 1.0 - ((double) freeMemory / maxMemory);

            ready.put("memoryUsagePercent", Math.round(memoryUsage * 100));

            if (memoryUsage > 0.95) {
                res.status(503); // Service Unavailable
                ready.put("status", "not_ready");
                ready.put("reason", "High memory usage");
            }

            return ready;
        }, json);

        // Server info
        get(basePath + "/info", (req, res) -> {
            Map<String, Object> info = new HashMap<String, Object>();

            // Version info
            String version = new Model("").getVersion();
            info.put("name", "LINE Solver");
            info.put("version", version);
            info.put("apiVersion", "v1");

            // Runtime info
            RuntimeMXBean runtimeBean = ManagementFactory.getRuntimeMXBean();
            info.put("javaVersion", System.getProperty("java.version"));
            info.put("javaVendor", System.getProperty("java.vendor"));
            info.put("uptimeMillis", runtimeBean.getUptime());

            // Memory info
            MemoryMXBean memoryBean = ManagementFactory.getMemoryMXBean();
            Map<String, Object> memory = new HashMap<String, Object>();
            memory.put("heapUsed", memoryBean.getHeapMemoryUsage().getUsed());
            memory.put("heapMax", memoryBean.getHeapMemoryUsage().getMax());
            memory.put("heapCommitted", memoryBean.getHeapMemoryUsage().getCommitted());
            info.put("memory", memory);

            // Processor info
            info.put("availableProcessors", Runtime.getRuntime().availableProcessors());

            // Capabilities
            Map<String, Object> capabilities = new HashMap<String, Object>();
            capabilities.put("inputFormats", new String[]{"jsim", "jsimg", "jsimw", "lqnx", "xml"});
            capabilities.put("outputFormats", new String[]{"json"});
            capabilities.put("solvers", new String[]{"mva", "ctmc", "fluid", "jmt", "nc", "ssa", "ln", "lqns", "des"});
            capabilities.put("analysisTypes", new String[]{"all", "avg", "sys"});
            info.put("capabilities", capabilities);

            return info;
        }, json);
    }
}
