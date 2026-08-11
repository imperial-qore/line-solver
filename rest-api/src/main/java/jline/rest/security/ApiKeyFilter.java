/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.rest.security;

import spark.Filter;
import spark.Request;
import spark.Response;

import java.util.HashSet;
import java.util.Set;

import static spark.Spark.halt;

/**
 * API key authentication filter.
 * Validates requests against configured API keys.
 */
public class ApiKeyFilter implements Filter {

    private static final String API_KEY_HEADER = "X-API-Key";
    private static final String API_KEY_PARAM = "api_key";

    private final Set<String> validApiKeys;
    private final boolean enabled;
    private final Set<String> excludedPaths;

    /**
     * Create a new API key filter.
     *
     * @param apiKeys Comma-separated list of valid API keys (null or empty to disable)
     */
    public ApiKeyFilter(String apiKeys) {
        this.validApiKeys = new HashSet<String>();
        this.excludedPaths = new HashSet<String>();

        // Always exclude health endpoints for Kubernetes probes
        excludedPaths.add("/api/v1/health");
        excludedPaths.add("/api/v1/ready");
        excludedPaths.add("/api/v1/metrics");

        if (apiKeys != null && !apiKeys.trim().isEmpty()) {
            this.enabled = true;
            for (String key : apiKeys.split(",")) {
                String trimmed = key.trim();
                if (!trimmed.isEmpty()) {
                    validApiKeys.add(trimmed);
                }
            }
        } else {
            this.enabled = false;
        }
    }

    /**
     * Create a disabled API key filter.
     */
    public ApiKeyFilter() {
        this(null);
    }

    @Override
    public void handle(Request request, Response response) throws Exception {
        if (!enabled) {
            return;
        }

        // Skip excluded paths
        String path = request.pathInfo();
        for (String excluded : excludedPaths) {
            if (path.startsWith(excluded)) {
                return;
            }
        }

        // Skip OPTIONS requests (CORS preflight)
        if ("OPTIONS".equals(request.requestMethod())) {
            return;
        }

        // Check for API key in header or query parameter
        String apiKey = request.headers(API_KEY_HEADER);
        if (apiKey == null || apiKey.isEmpty()) {
            apiKey = request.queryParams(API_KEY_PARAM);
        }

        if (apiKey == null || apiKey.isEmpty()) {
            halt(401, "{\"error\": \"API key required\", \"code\": \"UNAUTHORIZED\"}");
            return;
        }

        if (!validApiKeys.contains(apiKey)) {
            halt(403, "{\"error\": \"Invalid API key\", \"code\": \"FORBIDDEN\"}");
        }
    }

    /**
     * Add a path to exclude from authentication.
     *
     * @param path Path prefix to exclude
     */
    public void addExcludedPath(String path) {
        excludedPaths.add(path);
    }

    /**
     * Check if the filter is enabled.
     *
     * @return true if API key validation is enabled
     */
    public boolean isEnabled() {
        return enabled;
    }

    /**
     * Get the number of configured API keys.
     *
     * @return Number of valid API keys
     */
    public int getKeyCount() {
        return validApiKeys.size();
    }
}
