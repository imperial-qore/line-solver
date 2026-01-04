/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.rest;

import jline.lang.Model;
import jline.rest.jobs.JobManager;
import jline.rest.routes.AnalysisRoutes;
import jline.rest.routes.HealthRoutes;
import jline.rest.routes.JobRoutes;
import jline.rest.routes.MetricsRoutes;
import jline.rest.routes.ModelRoutes;
import jline.rest.routes.SolverRoutes;
import jline.rest.routes.SreRoutes;
import jline.rest.security.ApiKeyFilter;
import jline.rest.security.RateLimitFilter;
import jline.rest.util.JsonTransformer;

import static spark.Spark.*;

/**
 * LINE REST API Server.
 * Provides a REST API for solving queueing network models.
 *
 * Usage:
 *   java -cp jline.jar jline.rest.LineRestServer [--port PORT]
 *
 * Default port: 8080
 */
public class LineRestServer {

    private static final int DEFAULT_PORT = 8080;
    private static final String API_VERSION = "v1";
    private static final String API_BASE = "/api/" + API_VERSION;

    private final int port;
    private final JsonTransformer jsonTransformer;
    private final JobManager jobManager;
    private final ApiKeyFilter apiKeyFilter;
    private final RateLimitFilter rateLimitFilter;

    /**
     * Create a new REST server with the specified port.
     * @param port The port to listen on
     */
    public LineRestServer(int port) {
        this(port, null, 0, 60);
    }

    /**
     * Create a new REST server with security configuration.
     * @param port The port to listen on
     * @param apiKeys Comma-separated API keys (null to disable)
     * @param rateLimit Max requests per window (0 to disable)
     * @param rateLimitWindow Window size in seconds
     */
    public LineRestServer(int port, String apiKeys, int rateLimit, int rateLimitWindow) {
        this.port = port;
        this.jsonTransformer = new JsonTransformer();
        this.jobManager = new JobManager();
        this.apiKeyFilter = new ApiKeyFilter(apiKeys);
        this.rateLimitFilter = new RateLimitFilter(rateLimit, rateLimitWindow);
    }

    /**
     * Create a new REST server with the default port (8080).
     */
    public LineRestServer() {
        this(DEFAULT_PORT);
    }

    /**
     * Start the REST server.
     */
    public void start() {
        // Configure server
        port(port);

        // Enable CORS for all origins (configurable in production)
        enableCORS("*", "GET, POST, PUT, DELETE, OPTIONS", "Content-Type, Authorization, X-API-Key");

        // Apply security filters
        before(apiKeyFilter);
        before(rateLimitFilter);

        // Configure JSON response type
        before((request, response) -> {
            response.type("application/json");
        });

        // Register routes
        registerRoutes();

        // Global exception handling
        configureExceptionHandling();

        // Wait for server to initialize
        awaitInitialization();

        printStartupBanner();
    }

    /**
     * Stop the REST server.
     */
    public void stopServer() {
        jobManager.shutdown();
        AnalysisRoutes.shutdown();
        stop();
        awaitStop();
    }

    /**
     * Register all API routes.
     */
    private void registerRoutes() {
        // Health and info endpoints
        HealthRoutes.register(API_BASE, jsonTransformer);

        // Solver listing endpoint
        SolverRoutes.register(API_BASE, jsonTransformer);

        // Model solving endpoints
        ModelRoutes.register(API_BASE, jsonTransformer);

        // Async job endpoints
        JobRoutes.register(API_BASE, jsonTransformer, jobManager);

        // Analysis endpoints (what-if, sensitivity, bottleneck)
        AnalysisRoutes.register(API_BASE, jsonTransformer);

        // Prometheus metrics endpoint
        MetricsRoutes.register(API_BASE, jsonTransformer, jobManager);

        // SRE/DevOps integration (calibration, trace import)
        SreRoutes.register(API_BASE, jsonTransformer);
    }

    /**
     * Configure global exception handling.
     */
    private void configureExceptionHandling() {
        // Handle 404
        notFound((req, res) -> {
            res.type("application/json");
            return "{\"error\": \"Not found\", \"path\": \"" + req.pathInfo() + "\"}";
        });

        // Handle 500
        internalServerError((req, res) -> {
            res.type("application/json");
            return "{\"error\": \"Internal server error\"}";
        });

        // Handle all other exceptions
        exception(Exception.class, (e, req, res) -> {
            res.type("application/json");
            res.status(500);
            res.body("{\"error\": \"" + escapeJson(e.getMessage()) + "\"}");
        });
    }

    /**
     * Enable CORS for the specified origins.
     */
    private void enableCORS(String origin, String methods, String headers) {
        options("/*", (request, response) -> {
            String accessControlRequestHeaders = request.headers("Access-Control-Request-Headers");
            if (accessControlRequestHeaders != null) {
                response.header("Access-Control-Allow-Headers", accessControlRequestHeaders);
            }
            String accessControlRequestMethod = request.headers("Access-Control-Request-Method");
            if (accessControlRequestMethod != null) {
                response.header("Access-Control-Allow-Methods", accessControlRequestMethod);
            }
            return "OK";
        });

        before((request, response) -> {
            response.header("Access-Control-Allow-Origin", origin);
            response.header("Access-Control-Request-Method", methods);
            response.header("Access-Control-Allow-Headers", headers);
        });
    }

    /**
     * Print startup banner to console.
     */
    private void printStartupBanner() {
        String version = new Model("").getVersion();
        System.out.println("====================================================================");
        System.out.println("LINE Solver - REST API Server");
        System.out.println("Copyright (c) 2012-2026, QORE Lab, Imperial College London");
        System.out.println("Version " + version + ". All rights reserved.");
        System.out.println("====================================================================");
        System.out.println();
        System.out.println("REST API running on http://localhost:" + port + API_BASE);
        System.out.println();
        System.out.println("Endpoints:");
        System.out.println("  GET  " + API_BASE + "/health              - Health check");
        System.out.println("  GET  " + API_BASE + "/info                - Server info");
        System.out.println("  GET  " + API_BASE + "/solvers             - List available solvers");
        System.out.println("  POST " + API_BASE + "/models/solve        - Solve a model (sync)");
        System.out.println("  POST " + API_BASE + "/models/solve/async  - Solve a model (async)");
        System.out.println("  POST " + API_BASE + "/models/validate     - Validate a model");
        System.out.println("  POST " + API_BASE + "/models/convert      - Convert model format");
        System.out.println("  GET  " + API_BASE + "/jobs                - List all jobs");
        System.out.println("  GET  " + API_BASE + "/jobs/:id            - Get job status");
        System.out.println("  GET  " + API_BASE + "/jobs/:id/stream     - Stream job progress (SSE)");
        System.out.println("  DELETE " + API_BASE + "/jobs/:id          - Cancel a job");
        System.out.println("  POST " + API_BASE + "/analysis/whatif     - What-if parameter sweep");
        System.out.println("  POST " + API_BASE + "/analysis/sensitivity - Sensitivity analysis");
        System.out.println("  POST " + API_BASE + "/analysis/bottleneck - Bottleneck detection");
        System.out.println("  GET  " + API_BASE + "/metrics             - Prometheus metrics");
        System.out.println("  GET  " + API_BASE + "/metrics/json        - JSON metrics");
        System.out.println("  POST " + API_BASE + "/models/calibrate    - Calibrate from metrics");
        System.out.println("  POST " + API_BASE + "/traces/import       - Import traces to model");
        System.out.println();
        System.out.println("Press Ctrl+C to stop the server.");
    }

    /**
     * Escape a string for JSON output.
     */
    private static String escapeJson(String str) {
        if (str == null) {
            return "";
        }
        return str.replace("\\", "\\\\")
                  .replace("\"", "\\\"")
                  .replace("\n", "\\n")
                  .replace("\r", "\\r")
                  .replace("\t", "\\t");
    }

    /**
     * Main entry point.
     */
    public static void main(String[] args) {
        int port = DEFAULT_PORT;
        String apiKeys = System.getenv("LINE_API_KEYS");
        int rateLimit = parseEnvInt("LINE_RATE_LIMIT", 0);
        int rateLimitWindow = parseEnvInt("LINE_RATE_WINDOW", 60);

        // Parse command line arguments
        for (int i = 0; i < args.length; i++) {
            if ("--port".equals(args[i]) || "-p".equals(args[i])) {
                if (i + 1 < args.length) {
                    try {
                        port = Integer.parseInt(args[i + 1]);
                        if (port < 1 || port > 65535) {
                            System.err.println("Error: Port must be between 1 and 65535");
                            System.exit(1);
                        }
                        i++;
                    } catch (NumberFormatException e) {
                        System.err.println("Error: Invalid port number: " + args[i + 1]);
                        System.exit(1);
                    }
                } else {
                    System.err.println("Error: --port requires a value");
                    System.exit(1);
                }
            } else if ("--api-keys".equals(args[i])) {
                if (i + 1 < args.length) {
                    apiKeys = args[i + 1];
                    i++;
                } else {
                    System.err.println("Error: --api-keys requires a value");
                    System.exit(1);
                }
            } else if ("--rate-limit".equals(args[i])) {
                if (i + 1 < args.length) {
                    try {
                        rateLimit = Integer.parseInt(args[i + 1]);
                        i++;
                    } catch (NumberFormatException e) {
                        System.err.println("Error: Invalid rate limit: " + args[i + 1]);
                        System.exit(1);
                    }
                } else {
                    System.err.println("Error: --rate-limit requires a value");
                    System.exit(1);
                }
            } else if ("--help".equals(args[i]) || "-h".equals(args[i])) {
                printHelp();
                System.exit(0);
            }
        }

        // Start server
        LineRestServer server = new LineRestServer(port, apiKeys, rateLimit, rateLimitWindow);

        // Add shutdown hook for graceful shutdown
        Runtime.getRuntime().addShutdownHook(new Thread(() -> {
            System.out.println("\nShutting down REST server...");
            server.stopServer();
            System.out.println("Server stopped.");
        }));

        server.start();
    }

    /**
     * Print help message.
     */
    private static void printHelp() {
        System.out.println("LINE Solver REST API Server");
        System.out.println();
        System.out.println("Usage: java -cp jline.jar jline.rest.LineRestServer [OPTIONS]");
        System.out.println();
        System.out.println("Options:");
        System.out.println("  -p, --port PORT       Port to listen on (default: 8080)");
        System.out.println("  --api-keys KEYS       Comma-separated API keys for authentication");
        System.out.println("  --rate-limit N        Max requests per window (0 to disable)");
        System.out.println("  -h, --help            Show this help message");
        System.out.println();
        System.out.println("Environment variables:");
        System.out.println("  LINE_API_KEYS         API keys (alternative to --api-keys)");
        System.out.println("  LINE_RATE_LIMIT       Rate limit (alternative to --rate-limit)");
        System.out.println("  LINE_RATE_WINDOW      Rate limit window in seconds (default: 60)");
        System.out.println();
        System.out.println("Examples:");
        System.out.println("  java -cp jline.jar jline.rest.LineRestServer");
        System.out.println("  java -cp jline.jar jline.rest.LineRestServer --port 9090");
        System.out.println("  java -cp jline.jar jline.rest.LineRestServer --api-keys key1,key2");
        System.out.println("  java -cp jline.jar jline.rest.LineRestServer --rate-limit 100");
    }

    /**
     * Parse an integer from environment variable with default.
     */
    private static int parseEnvInt(String name, int defaultValue) {
        String value = System.getenv(name);
        if (value == null || value.isEmpty()) {
            return defaultValue;
        }
        try {
            return Integer.parseInt(value);
        } catch (NumberFormatException e) {
            return defaultValue;
        }
    }
}
