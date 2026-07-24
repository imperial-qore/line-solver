/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.rest.security;

import spark.Filter;
import spark.Request;
import spark.Response;

import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLong;

import static spark.Spark.halt;

/**
 * Rate limiting filter using token bucket algorithm.
 * Limits requests per IP address or API key.
 */
public class RateLimitFilter implements Filter {

    private static final String RATE_LIMIT_HEADER = "X-RateLimit-Limit";
    private static final String RATE_LIMIT_REMAINING_HEADER = "X-RateLimit-Remaining";
    private static final String RATE_LIMIT_RESET_HEADER = "X-RateLimit-Reset";
    private static final String RETRY_AFTER_HEADER = "Retry-After";

    private final int maxRequests;
    private final long windowMs;
    private final boolean enabled;
    private final Map<String, RateLimitBucket> buckets;

    /**
     * Create a new rate limit filter.
     *
     * @param maxRequests Maximum requests per window (0 to disable)
     * @param windowSeconds Time window in seconds
     */
    public RateLimitFilter(int maxRequests, int windowSeconds) {
        this.maxRequests = maxRequests;
        this.windowMs = windowSeconds * 1000L;
        this.enabled = maxRequests > 0;
        this.buckets = new ConcurrentHashMap<String, RateLimitBucket>();

        // Start cleanup thread
        if (enabled) {
            Thread cleanupThread = new Thread(new Runnable() {
                @Override
                public void run() {
                    while (true) {
                        try {
                            Thread.sleep(60000); // Clean up every minute
                            cleanupExpiredBuckets();
                        } catch (InterruptedException e) {
                            Thread.currentThread().interrupt();
                            break;
                        }
                    }
                }
            }, "rate-limit-cleanup");
            cleanupThread.setDaemon(true);
            cleanupThread.start();
        }
    }

    /**
     * Create a disabled rate limit filter.
     */
    public RateLimitFilter() {
        this(0, 60);
    }

    @Override
    public void handle(Request request, Response response) throws Exception {
        if (!enabled) {
            return;
        }

        // Skip OPTIONS requests (CORS preflight)
        if ("OPTIONS".equals(request.requestMethod())) {
            return;
        }

        // Get client identifier (API key or IP)
        String clientId = getClientId(request);
        RateLimitBucket bucket = getBucket(clientId);

        // Try to consume a token
        long now = System.currentTimeMillis();
        boolean allowed = bucket.tryConsume(now);

        // Set rate limit headers
        response.header(RATE_LIMIT_HEADER, String.valueOf(maxRequests));
        response.header(RATE_LIMIT_REMAINING_HEADER, String.valueOf(bucket.getRemaining()));
        response.header(RATE_LIMIT_RESET_HEADER, String.valueOf(bucket.getResetTime() / 1000));

        if (!allowed) {
            long retryAfter = (bucket.getResetTime() - now) / 1000;
            if (retryAfter < 1) retryAfter = 1;
            response.header(RETRY_AFTER_HEADER, String.valueOf(retryAfter));
            halt(429, "{\"error\": \"Rate limit exceeded\", \"code\": \"TOO_MANY_REQUESTS\", " +
                      "\"retryAfter\": " + retryAfter + "}");
        }
    }

    /**
     * Get the client identifier from the request.
     */
    private String getClientId(Request request) {
        // Try API key first
        String apiKey = request.headers("X-API-Key");
        if (apiKey != null && !apiKey.isEmpty()) {
            return "key:" + apiKey;
        }

        // Fall back to IP address
        String ip = request.headers("X-Forwarded-For");
        if (ip != null && !ip.isEmpty()) {
            // Take first IP if multiple
            int comma = ip.indexOf(',');
            if (comma > 0) {
                ip = ip.substring(0, comma).trim();
            }
        } else {
            ip = request.ip();
        }
        return "ip:" + ip;
    }

    /**
     * Get or create a bucket for the client.
     */
    private RateLimitBucket getBucket(String clientId) {
        RateLimitBucket bucket = buckets.get(clientId);
        if (bucket == null) {
            bucket = new RateLimitBucket(maxRequests, windowMs);
            RateLimitBucket existing = buckets.putIfAbsent(clientId, bucket);
            if (existing != null) {
                bucket = existing;
            }
        }
        return bucket;
    }

    /**
     * Clean up expired buckets.
     */
    private void cleanupExpiredBuckets() {
        long now = System.currentTimeMillis();
        for (Map.Entry<String, RateLimitBucket> entry : buckets.entrySet()) {
            if (entry.getValue().isExpired(now)) {
                buckets.remove(entry.getKey());
            }
        }
    }

    /**
     * Check if the filter is enabled.
     */
    public boolean isEnabled() {
        return enabled;
    }

    /**
     * Get the maximum requests per window.
     */
    public int getMaxRequests() {
        return maxRequests;
    }

    /**
     * Get the window size in milliseconds.
     */
    public long getWindowMs() {
        return windowMs;
    }

    /**
     * Token bucket for rate limiting.
     */
    private static class RateLimitBucket {
        private final int maxTokens;
        private final long windowMs;
        private final AtomicInteger tokens;
        private final AtomicLong windowStart;
        private final AtomicLong lastAccess;

        RateLimitBucket(int maxTokens, long windowMs) {
            this.maxTokens = maxTokens;
            this.windowMs = windowMs;
            this.tokens = new AtomicInteger(maxTokens);
            this.windowStart = new AtomicLong(System.currentTimeMillis());
            this.lastAccess = new AtomicLong(System.currentTimeMillis());
        }

        boolean tryConsume(long now) {
            lastAccess.set(now);

            // Check if window has expired and reset
            long start = windowStart.get();
            if (now - start >= windowMs) {
                // Reset the window
                if (windowStart.compareAndSet(start, now)) {
                    tokens.set(maxTokens);
                }
            }

            // Try to consume a token
            while (true) {
                int current = tokens.get();
                if (current <= 0) {
                    return false;
                }
                if (tokens.compareAndSet(current, current - 1)) {
                    return true;
                }
            }
        }

        int getRemaining() {
            return Math.max(0, tokens.get());
        }

        long getResetTime() {
            return windowStart.get() + windowMs;
        }

        boolean isExpired(long now) {
            return now - lastAccess.get() > windowMs * 2;
        }
    }
}
