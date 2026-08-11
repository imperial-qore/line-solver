/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang;

import jline.lang.nodes.Cache;

/**
 * A set of cacheable items
 */
public class ItemSet extends NetworkElement {
    private final int nitems;
    private final Cache reference;
    private final boolean replicable;
    private int index;

    public ItemSet(Network model, String name, int nitems, Cache reference) {
        super(name);
        this.nitems = nitems;
        this.index = 0;
        this.replicable = false;
        if (reference == null) {
            throw new RuntimeException("ItemSet must be pinned to a Cache");
        }
        this.reference = reference;
        model.addItemSet(this);
    }

    public int getIndex() {
        return index;
    }

    public void setIndex(int index) {
        this.index = index;
    }

    public int getNumberOfItems() {
        return this.nitems;
    }

    public boolean hasReplicableItems() {
        return this.replicable;
    }
}
