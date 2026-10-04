/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang;

import java.io.Serializable;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;

/**
 * Metadata container held by every {@link Model}, twin of the untyped
 * {@code attribute} struct of {@code matlab/src/lang/Model.m}.
 *
 * MATLAB stores model metadata as dynamically-named struct fields, so the Java
 * twin is a string-keyed map. Subclasses needing a typed container extend this
 * one ({@link NetworkAttribute} does), which keeps the covariant
 * {@code getAttribute()} override on {@link Network} legal.
 */
public class ModelAttribute implements Serializable {
    private static final long serialVersionUID = 1L;

    private final Map<String, Object> metadata;

    public ModelAttribute() {
        this.metadata = new LinkedHashMap<String, Object>();
    }

    /**
     * Reads a metadata field.
     *
     * @param key the field name
     * @return the stored value, or null when the field is unset
     */
    public Object get(String key) {
        return this.metadata.get(key);
    }

    /**
     * Writes a metadata field, overwriting any previous value.
     *
     * @param key   the field name
     * @param value the value to store
     */
    public void put(String key, Object value) {
        this.metadata.put(key, value);
    }

    /**
     * Tests whether a metadata field is set.
     *
     * @param key the field name
     * @return true if the field is present
     */
    public boolean has(String key) {
        return this.metadata.containsKey(key);
    }

    /**
     * Removes a metadata field.
     *
     * @param key the field name
     * @return the removed value, or null when the field was unset
     */
    public Object remove(String key) {
        return this.metadata.remove(key);
    }

    /**
     * Lists the metadata fields currently set, in insertion order.
     *
     * @return the field names
     */
    public Set<String> keys() {
        return this.metadata.keySet();
    }

    /**
     * Removes every metadata field.
     */
    public void clear() {
        this.metadata.clear();
    }
}
