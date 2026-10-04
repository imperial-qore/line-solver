/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.rest.util;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import spark.ResponseTransformer;

/**
 * JSON response transformer using Gson.
 * Converts Java objects to JSON for REST API responses.
 */
public class JsonTransformer implements ResponseTransformer {

    private final Gson gson;

    /**
     * Create a new JSON transformer with default settings.
     */
    public JsonTransformer() {
        this.gson = new GsonBuilder()
                .serializeNulls()
                .serializeSpecialFloatingPointValues()
                .create();
    }

    /**
     * Create a new JSON transformer with pretty printing.
     * @param prettyPrint Whether to format JSON with indentation
     */
    public JsonTransformer(boolean prettyPrint) {
        GsonBuilder builder = new GsonBuilder()
                .serializeNulls()
                .serializeSpecialFloatingPointValues();
        if (prettyPrint) {
            builder.setPrettyPrinting();
        }
        this.gson = builder.create();
    }

    @Override
    public String render(Object model) {
        return gson.toJson(model);
    }

    /**
     * Parse JSON string to object.
     * @param json The JSON string
     * @param clazz The target class
     * @param <T> The target type
     * @return The parsed object
     */
    public <T> T fromJson(String json, Class<T> clazz) {
        return gson.fromJson(json, clazz);
    }

    /**
     * Convert object to JSON string.
     * @param object The object to convert
     * @return The JSON string
     */
    public String toJson(Object object) {
        return gson.toJson(object);
    }

    /**
     * Get the underlying Gson instance.
     * @return The Gson instance
     */
    public Gson getGson() {
        return gson;
    }
}
