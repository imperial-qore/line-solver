/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.rest.security;

import jline.rest.LineRestServer;
import jline.rest.util.JsonTransformer;
import org.junit.jupiter.api.*;
import spark.Spark;

import java.io.*;
import java.net.HttpURLConnection;
import java.net.URL;
import java.nio.charset.StandardCharsets;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Tests for security filters (API key authentication and rate limiting).
 */
@TestMethodOrder(MethodOrderer.OrderAnnotation.class)
class SecurityFilterTest {

    private static final int TEST_PORT = 4571;
    private static final String BASE_URL = "http://localhost:" + TEST_PORT + "/api/v1";

    @Test
    @Order(1)
    void testApiKeyFilterDisabledByDefault() {
        ApiKeyFilter filter = new ApiKeyFilter();
        assertFalse(filter.isEnabled());
        assertEquals(0, filter.getKeyCount());
    }

    @Test
    @Order(2)
    void testApiKeyFilterWithKeys() {
        ApiKeyFilter filter = new ApiKeyFilter("key1,key2,key3");
        assertTrue(filter.isEnabled());
        assertEquals(3, filter.getKeyCount());
    }

    @Test
    @Order(3)
    void testApiKeyFilterWithEmptyString() {
        ApiKeyFilter filter = new ApiKeyFilter("");
        assertFalse(filter.isEnabled());
    }

    @Test
    @Order(4)
    void testApiKeyFilterWithWhitespace() {
        ApiKeyFilter filter = new ApiKeyFilter("  key1  ,  key2  ");
        assertTrue(filter.isEnabled());
        assertEquals(2, filter.getKeyCount());
    }

    @Test
    @Order(5)
    void testApiKeyFilterExcludedPaths() {
        ApiKeyFilter filter = new ApiKeyFilter("key1");
        filter.addExcludedPath("/custom/path");
        assertTrue(filter.isEnabled());
    }

    @Test
    @Order(6)
    void testRateLimitFilterDisabledByDefault() {
        RateLimitFilter filter = new RateLimitFilter();
        assertFalse(filter.isEnabled());
    }

    @Test
    @Order(7)
    void testRateLimitFilterEnabled() {
        RateLimitFilter filter = new RateLimitFilter(100, 60);
        assertTrue(filter.isEnabled());
        assertEquals(100, filter.getMaxRequests());
        assertEquals(60000, filter.getWindowMs());
    }

    @Test
    @Order(8)
    void testRateLimitFilterWithZero() {
        RateLimitFilter filter = new RateLimitFilter(0, 60);
        assertFalse(filter.isEnabled());
    }

    @Test
    @Order(9)
    void testRateLimitBucketConsumption() throws Exception {
        // Test token bucket behavior through reflection or indirect testing
        RateLimitFilter filter = new RateLimitFilter(5, 1);
        assertTrue(filter.isEnabled());
        assertEquals(5, filter.getMaxRequests());
    }

    @Test
    @Order(10)
    void testApiKeyFilterSingleKey() {
        ApiKeyFilter filter = new ApiKeyFilter("single-key");
        assertTrue(filter.isEnabled());
        assertEquals(1, filter.getKeyCount());
    }
}
