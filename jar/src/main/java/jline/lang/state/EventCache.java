/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.state;

import jline.io.Ret;

import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

/**
 * A class storing events
 */
public class EventCache {

    private final Map<EventCacheKey, Ret.EventResult> results;
    private boolean enabled;

    public EventCache(boolean para, boolean enabled) {
        if (para) {
            results = new ConcurrentHashMap<>();
        } else {
            results = new HashMap<>();
        }
        this.enabled = enabled;
    }

    public boolean contains(EventCacheKey key) {
        return results.containsKey(key);
    }

    public Ret.EventResult get(EventCacheKey key) {
        return results.get(key);
    }

    public boolean isEnabled() {
        return enabled;
    }

    public void setEnabled(boolean enabled) {
        this.enabled = enabled;
    }

    public void put(EventCacheKey key, Ret.EventResult result) {
        results.put(key, result);
    }
}
