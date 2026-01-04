/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.rest.model;

/**
 * Represents the input model in a solve request.
 */
public class ModelInput {

    /**
     * The format of the model content.
     * Supported values: jsim, jsimg, jsimw, lqnx, xml
     */
    private String format;

    /**
     * The model content, either as raw XML/text or base64-encoded.
     */
    private String content;

    /**
     * Whether the content is base64-encoded.
     */
    private boolean base64;

    /**
     * Default constructor for JSON deserialization.
     */
    public ModelInput() {
    }

    /**
     * Create a new ModelInput with the specified format and content.
     * @param format The model format
     * @param content The model content
     */
    public ModelInput(String format, String content) {
        this.format = format;
        this.content = content;
        this.base64 = false;
    }

    /**
     * Create a new ModelInput with all fields specified.
     * @param format The model format
     * @param content The model content
     * @param base64 Whether the content is base64-encoded
     */
    public ModelInput(String format, String content, boolean base64) {
        this.format = format;
        this.content = content;
        this.base64 = base64;
    }

    public String getFormat() {
        return format;
    }

    public void setFormat(String format) {
        this.format = format;
    }

    public String getContent() {
        return content;
    }

    public void setContent(String content) {
        this.content = content;
    }

    public boolean isBase64() {
        return base64;
    }

    public void setBase64(boolean base64) {
        this.base64 = base64;
    }

    /**
     * Validate the model input.
     * @return null if valid, error message if invalid
     */
    public String validate() {
        if (format == null || format.trim().isEmpty()) {
            return "Model format is required";
        }
        if (content == null || content.trim().isEmpty()) {
            return "Model content is required";
        }
        String[] validFormats = {"jsim", "jsimg", "jsimw", "lqnx", "xml"};
        boolean validFormat = false;
        for (String valid : validFormats) {
            if (valid.equals(format)) {
                validFormat = true;
                break;
            }
        }
        if (!validFormat) {
            return "Invalid format: " + format + ". Valid formats: jsim, jsimg, jsimw, lqnx, xml";
        }
        return null;
    }
}
