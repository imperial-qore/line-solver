package jline.lang.state;

import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

public class EventCache {

    private boolean enabled;
    private Map<EventCacheKey, EventResult> results;

    public EventCache(boolean para, boolean enabled) {
        if (para) {
            results = new ConcurrentHashMap<>();
        }
        else {
            results = new HashMap<>();
        }
        this.enabled = enabled;
    }

    public boolean contains(EventCacheKey key) {
        return results.containsKey(key);
    }

    public EventResult get(EventCacheKey key) {
        return results.get(key);
    }

    public void put(EventCacheKey key, EventResult result) {
        results.put(key, result);
    }

    public void setEnabled(boolean enabled) {
        this.enabled = enabled;
    }

    public boolean isEnabled(){
        return enabled;
    }
}
