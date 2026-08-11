/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.rest.routes;

import com.google.gson.GsonBuilder;
import com.google.gson.JsonObject;
import jline.io.LineModelIO;
import jline.io.M2M;
import jline.lang.Model;
import jline.lang.Network;
import jline.lang.layered.LayeredNetwork;
import jline.rest.handlers.SolveHandler;
import jline.rest.model.ErrorResponse;
import jline.rest.model.ModelInput;
import jline.rest.model.SolveRequest;
import jline.rest.model.SolveResponse;
import jline.rest.util.JsonTransformer;
import jline.rest.util.ModelParser;
import jline.rest.util.ModelParser.ModelParseException;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.Map;

import static spark.Spark.post;

/**
 * Model-related endpoints for the REST API.
 * Provides endpoints for solving, validating, and converting models.
 */
public class ModelRoutes {

    private static final SolveHandler solveHandler = new SolveHandler();
    private static final ModelParser modelParser = new ModelParser();

    /**
     * Register model routes.
     *
     * @param basePath The API base path (e.g., "/api/v1")
     * @param json The JSON transformer
     */
    public static void register(String basePath, JsonTransformer json) {
        // Solve a model (synchronous)
        post(basePath + "/models/solve", (req, res) -> {
            try {
                SolveRequest request = json.fromJson(req.body(), SolveRequest.class);
                if (request == null) {
                    res.status(400);
                    return new ErrorResponse("INVALID_REQUEST", "Request body is required");
                }

                SolveResponse response = solveHandler.solve(request);

                if ("failed".equals(response.getStatus())) {
                    res.status(400);
                }

                return response;
            } catch (Exception e) {
                res.status(500);
                return new ErrorResponse("SOLVER_ERROR", "Unexpected error: " + e.getMessage());
            }
        }, json);

        // Validate a model
        post(basePath + "/models/validate", (req, res) -> {
            try {
                ModelInput input = json.fromJson(req.body(), ModelInput.class);
                if (input == null) {
                    res.status(400);
                    return new ErrorResponse("INVALID_REQUEST", "Request body is required");
                }

                // Validate input fields
                String validationError = input.validate();
                if (validationError != null) {
                    res.status(400);
                    Map<String, Object> result = new HashMap<String, Object>();
                    result.put("valid", false);
                    result.put("error", validationError);
                    return result;
                }

                // Try to parse the model
                try {
                    Model model = modelParser.parse(input);

                    Map<String, Object> result = new HashMap<String, Object>();
                    result.put("valid", true);

                    // Add model info
                    Map<String, Object> modelInfo = new HashMap<String, Object>();
                    modelInfo.put("name", model.getName());
                    modelInfo.put("type", model instanceof LayeredNetwork ? "layered" : "queueing");

                    if (model instanceof Network) {
                        Network network = (Network) model;
                        modelInfo.put("stations", network.getNumberOfStations());
                        modelInfo.put("classes", network.getNumberOfClasses());
                        modelInfo.put("chains", network.getNumberOfChains());
                    }

                    result.put("model", modelInfo);
                    return result;
                } catch (ModelParseException e) {
                    Map<String, Object> result = new HashMap<String, Object>();
                    result.put("valid", false);
                    result.put("error", e.getMessage());
                    return result;
                }
            } catch (Exception e) {
                res.status(500);
                return new ErrorResponse("VALIDATION_ERROR", "Unexpected error: " + e.getMessage());
            }
        }, json);

        // Convert between model formats
        post(basePath + "/models/convert", (req, res) -> {
            try {
                Map<String, Object> request = json.fromJson(req.body(), Map.class);
                if (request == null) {
                    res.status(400);
                    return new ErrorResponse("INVALID_REQUEST", "Request body is required");
                }

                // Get input model
                Map<String, Object> modelMap = (Map<String, Object>) request.get("model");
                if (modelMap == null) {
                    res.status(400);
                    return new ErrorResponse("INVALID_REQUEST", "Model is required");
                }

                ModelInput input = new ModelInput();
                input.setFormat((String) modelMap.get("format"));
                input.setContent((String) modelMap.get("content"));
                Object base64Obj = modelMap.get("base64");
                if (base64Obj != null) {
                    input.setBase64((Boolean) base64Obj);
                }

                // Get target format
                String targetFormat = (String) request.get("targetFormat");
                if (targetFormat == null) {
                    res.status(400);
                    return new ErrorResponse("INVALID_REQUEST", "targetFormat is required");
                }

                // Validate target format
                if (!isValidFormat(targetFormat)) {
                    res.status(400);
                    return new ErrorResponse("INVALID_FORMAT",
                        "Invalid target format: " + targetFormat + ". Valid formats: jsimg, lqnx");
                }

                // Parse input model
                Model model;
                try {
                    model = modelParser.parse(input);
                } catch (ModelParseException e) {
                    res.status(400);
                    return new ErrorResponse("PARSE_ERROR", e.getMessage());
                }

                // Convert to target format
                String convertedContent = convertModel(model, targetFormat);
                if (convertedContent == null) {
                    res.status(400);
                    return new ErrorResponse("CONVERSION_ERROR",
                        "Cannot convert model to format: " + targetFormat);
                }

                Map<String, Object> result = new HashMap<String, Object>();
                result.put("format", targetFormat);
                result.put("content", convertedContent);
                return result;

            } catch (Exception e) {
                res.status(500);
                return new ErrorResponse("CONVERSION_ERROR", "Unexpected error: " + e.getMessage());
            }
        }, json);

        // Describe a model (return line-model JSON)
        // Note: no JsonTransformer — we return raw JSON strings to avoid double-encoding
        post(basePath + "/models/describe", (req, res) -> {
            res.type("application/json");
            try {
                ModelInput input = json.fromJson(req.body(), ModelInput.class);
                if (input == null) {
                    res.status(400);
                    return json.toJson(new ErrorResponse("INVALID_REQUEST", "Request body is required"));
                }

                String validationError = input.validate();
                if (validationError != null) {
                    res.status(400);
                    return json.toJson(new ErrorResponse("INVALID_INPUT", validationError));
                }

                Model model;
                try {
                    model = modelParser.parse(input);
                } catch (ModelParseException e) {
                    res.status(400);
                    return json.toJson(new ErrorResponse("PARSE_ERROR", e.getMessage()));
                }

                JsonObject doc;
                if (model instanceof Network) {
                    doc = LineModelIO.toJsonObject((Network) model);
                } else if (model instanceof LayeredNetwork) {
                    doc = LineModelIO.toJsonObject((LayeredNetwork) model);
                } else {
                    res.status(400);
                    return json.toJson(new ErrorResponse("UNSUPPORTED_TYPE", "Unsupported model type: " + model.getClass().getSimpleName()));
                }

                return doc.toString();
            } catch (Exception e) {
                res.status(500);
                return json.toJson(new ErrorResponse("DESCRIBE_ERROR", "Unexpected error: " + e.getMessage()));
            }
        });
    }

    /**
     * Check if the format is valid for conversion.
     */
    private static boolean isValidFormat(String format) {
        return "jsimg".equals(format) || "lqnx".equals(format) || "json".equals(format);
    }

    /**
     * Convert a model to the target format.
     */
    private static String convertModel(Model model, String targetFormat) {
        // JSON format works for any model type
        if ("json".equals(targetFormat)) {
            if (model instanceof Network) {
                JsonObject doc = LineModelIO.toJsonObject((Network) model);
                return new GsonBuilder().setPrettyPrinting().serializeSpecialFloatingPointValues().create().toJson(doc);
            } else if (model instanceof LayeredNetwork) {
                JsonObject doc = LineModelIO.toJsonObject((LayeredNetwork) model);
                return new GsonBuilder().setPrettyPrinting().serializeSpecialFloatingPointValues().create().toJson(doc);
            }
            return null;
        }

        if (!(model instanceof Network)) {
            return null;
        }

        Network network = (Network) model;
        M2M m2m = new M2M();

        Path tempFile = null;
        try {
            if ("jsimg".equals(targetFormat)) {
                tempFile = Files.createTempFile("line_convert_", ".jsimg");
                m2m.LINE2JSIMG(network, tempFile.toString());
                return new String(Files.readAllBytes(tempFile), StandardCharsets.UTF_8);
            }
        } catch (IOException e) {
            return null;
        } finally {
            if (tempFile != null) {
                try {
                    Files.deleteIfExists(tempFile);
                } catch (IOException e) {
                    // Ignore cleanup errors
                }
            }
        }

        return null;
    }
}
