/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.rest.routes;

import jline.rest.handlers.WhatIfHandler;
import jline.rest.model.*;
import jline.rest.util.JsonTransformer;

import static spark.Spark.*;

/**
 * Analysis endpoints for what-if, sensitivity, and bottleneck analysis.
 */
public class AnalysisRoutes {

    private static WhatIfHandler whatIfHandler;

    /**
     * Register analysis routes.
     *
     * @param basePath The API base path (e.g., "/api/v1")
     * @param json     The JSON transformer
     */
    public static void register(String basePath, JsonTransformer json) {
        whatIfHandler = new WhatIfHandler();

        // What-if parameter sweep
        post(basePath + "/analysis/whatif", (req, res) -> {
            try {
                WhatIfRequest request = json.fromJson(req.body(), WhatIfRequest.class);
                if (request == null) {
                    res.status(400);
                    return new ErrorResponse("INVALID_REQUEST", "Request body is required");
                }

                String validationError = request.validate();
                if (validationError != null) {
                    res.status(400);
                    return new ErrorResponse("VALIDATION_ERROR", validationError);
                }

                WhatIfResponse response = whatIfHandler.whatIf(request);

                if ("failed".equals(response.getStatus()) && response.getResults() == null) {
                    res.status(500);
                }

                return response;
            } catch (Exception e) {
                res.status(500);
                return new ErrorResponse("WHATIF_ERROR", "What-if analysis failed: " + e.getMessage());
            }
        }, json);

        // Sensitivity analysis
        post(basePath + "/analysis/sensitivity", (req, res) -> {
            try {
                SensitivityRequest request = json.fromJson(req.body(), SensitivityRequest.class);
                if (request == null) {
                    res.status(400);
                    return new ErrorResponse("INVALID_REQUEST", "Request body is required");
                }

                String validationError = request.validate();
                if (validationError != null) {
                    res.status(400);
                    return new ErrorResponse("VALIDATION_ERROR", validationError);
                }

                SensitivityResponse response = whatIfHandler.sensitivity(request);

                if ("failed".equals(response.getStatus())) {
                    res.status(500);
                }

                return response;
            } catch (Exception e) {
                res.status(500);
                return new ErrorResponse("SENSITIVITY_ERROR", "Sensitivity analysis failed: " + e.getMessage());
            }
        }, json);

        // Bottleneck detection
        post(basePath + "/analysis/bottleneck", (req, res) -> {
            try {
                BottleneckRequest request = json.fromJson(req.body(), BottleneckRequest.class);
                if (request == null) {
                    res.status(400);
                    return new ErrorResponse("INVALID_REQUEST", "Request body is required");
                }

                String validationError = request.validate();
                if (validationError != null) {
                    res.status(400);
                    return new ErrorResponse("VALIDATION_ERROR", validationError);
                }

                BottleneckResponse response = whatIfHandler.bottleneck(request);

                if ("failed".equals(response.getStatus())) {
                    res.status(500);
                }

                return response;
            } catch (Exception e) {
                res.status(500);
                return new ErrorResponse("BOTTLENECK_ERROR", "Bottleneck analysis failed: " + e.getMessage());
            }
        }, json);
    }

    /**
     * Get the WhatIfHandler instance.
     *
     * @return The handler
     */
    public static WhatIfHandler getHandler() {
        return whatIfHandler;
    }

    /**
     * Shutdown the handler.
     */
    public static void shutdown() {
        if (whatIfHandler != null) {
            whatIfHandler.shutdown();
        }
    }
}
