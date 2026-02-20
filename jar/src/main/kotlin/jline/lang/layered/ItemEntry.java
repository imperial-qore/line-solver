/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.layered;

import jline.lang.processes.Distribution;
import org.apache.commons.lang3.SerializationUtils;

/**
 * A caching service that gives access to items
 */
public class ItemEntry extends Entry {
    protected int cardinality;
    protected Distribution popularity;

    public ItemEntry(LayeredNetwork model, String name, int cardinality, Distribution distribution) {
        super(model, name);
        this.cardinality = cardinality;
        if (distribution.isDiscrete()) {
            this.popularity = SerializationUtils.clone(distribution);
        }
    }

    public int getCardinality() {
        return cardinality;
    }

    public Distribution getPopularity() {
        return popularity;
    }

    @Override
    public ItemEntry on(Task newParent) {
        newParent.addEntry(this);
        this.parent = newParent;
        // Update CacheTask's items to match ItemEntry's cardinality
        if (newParent instanceof CacheTask) {
            CacheTask cacheTask = (CacheTask) newParent;
            if (cacheTask.getItems() < this.cardinality) {
                cacheTask.setItems(this.cardinality);
            }
        }
        return this;
    }
}
