/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.rest.model;

/**
 * Standard error response for the REST API.
 */
public class ErrorResponse {

    /**
     * Error code (e.g., "VALIDATION_ERROR", "SOLVER_ERROR").
     */
    private String code;

    /**
     * Human-readable error message.
     */
    private String error;

    /**
     * Optional additional details.
     */
    private String details;

    /**
     * Default constructor for JSON deserialization.
     */
    public ErrorResponse() {
    }

    /**
     * Create an error response with a message.
     * @param error The error message
     */
    public ErrorResponse(String error) {
        this.error = error;
    }

    /**
     * Create an error response with code and message.
     * @param code The error code
     * @param error The error message
     */
    public ErrorResponse(String code, String error) {
        this.code = code;
        this.error = error;
    }

    /**
     * Create an error response with code, message, and details.
     * @param code The error code
     * @param error The error message
     * @param details Additional details
     */
    public ErrorResponse(String code, String error, String details) {
        this.code = code;
        this.error = error;
        this.details = details;
    }

    public String getCode() {
        return code;
    }

    public void setCode(String code) {
        this.code = code;
    }

    public String getError() {
        return error;
    }

    public void setError(String error) {
        this.error = error;
    }

    public String getDetails() {
        return details;
    }

    public void setDetails(String details) {
        this.details = details;
    }
}
