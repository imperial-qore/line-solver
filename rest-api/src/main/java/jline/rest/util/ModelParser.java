/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.rest.util;

import jline.io.M2M;
import jline.lang.Model;
import jline.lang.Network;
import jline.lang.layered.LayeredNetwork;
import jline.rest.model.ModelInput;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Base64;

/**
 * Utility class for parsing model input from REST API requests.
 */
public class ModelParser {

    private final M2M m2m;

    /**
     * Create a new model parser.
     */
    public ModelParser() {
        this.m2m = new M2M();
    }

    /**
     * Parse a ModelInput into a LINE Model.
     *
     * @param input The model input from the API request
     * @return The parsed Model (Network or LayeredNetwork)
     * @throws ModelParseException if parsing fails
     */
    public Model parse(ModelInput input) throws ModelParseException {
        if (input == null) {
            throw new ModelParseException("Model input is null");
        }

        String validationError = input.validate();
        if (validationError != null) {
            throw new ModelParseException(validationError);
        }

        String content = input.getContent();
        String format = input.getFormat();

        // Decode base64 if needed
        if (input.isBase64()) {
            try {
                byte[] decoded = Base64.getDecoder().decode(content);
                content = new String(decoded, StandardCharsets.UTF_8);
            } catch (IllegalArgumentException e) {
                throw new ModelParseException("Invalid base64 encoding: " + e.getMessage());
            }
        }

        // Write content to temp file for parsing
        Path tempFile = null;
        try {
            tempFile = Files.createTempFile("line_rest_", "." + format);
            Files.write(tempFile, content.getBytes(StandardCharsets.UTF_8));

            // Parse based on format
            switch (format) {
                case "jsim":
                case "jsimg":
                case "jsimw":
                    return m2m.JSIM2LINE(tempFile.toString());
                case "lqnx":
                case "xml":
                    return m2m.LQN2LINE(tempFile.toString());
                default:
                    throw new ModelParseException("Unsupported format: " + format);
            }
        } catch (IOException e) {
            throw new ModelParseException("Failed to create temp file: " + e.getMessage());
        } catch (Exception e) {
            throw new ModelParseException("Failed to parse model: " + e.getMessage());
        } finally {
            // Clean up temp file
            if (tempFile != null) {
                try {
                    Files.deleteIfExists(tempFile);
                } catch (IOException e) {
                    // Ignore cleanup errors
                }
            }
        }
    }

    /**
     * Check if the parsed model is a LayeredNetwork.
     *
     * @param model The parsed model
     * @return true if it's a LayeredNetwork, false otherwise
     */
    public static boolean isLayeredNetwork(Model model) {
        return model instanceof LayeredNetwork;
    }

    /**
     * Check if the parsed model is a Network (queueing network).
     *
     * @param model The parsed model
     * @return true if it's a Network, false otherwise
     */
    public static boolean isNetwork(Model model) {
        return model instanceof Network;
    }

    /**
     * Cast the model to Network.
     *
     * @param model The model to cast
     * @return The model as Network
     * @throws ModelParseException if the model is not a Network
     */
    public static Network asNetwork(Model model) throws ModelParseException {
        if (!(model instanceof Network)) {
            throw new ModelParseException("Model is not a queueing network");
        }
        return (Network) model;
    }

    /**
     * Cast the model to LayeredNetwork.
     *
     * @param model The model to cast
     * @return The model as LayeredNetwork
     * @throws ModelParseException if the model is not a LayeredNetwork
     */
    public static LayeredNetwork asLayeredNetwork(Model model) throws ModelParseException {
        if (!(model instanceof LayeredNetwork)) {
            throw new ModelParseException("Model is not a layered network");
        }
        return (LayeredNetwork) model;
    }

    /**
     * Exception thrown when model parsing fails.
     */
    public static class ModelParseException extends Exception {
        public ModelParseException(String message) {
            super(message);
        }

        public ModelParseException(String message, Throwable cause) {
            super(message, cause);
        }
    }
}
