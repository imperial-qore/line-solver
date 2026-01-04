/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.layered;

import jline.lang.constant.ReplacementStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.processes.Distribution;
import jline.lang.processes.Immediate;

/**
 * A task that offers caching services
 */
public class CacheTask extends Task {
    protected int items;
    protected int[] itemLevelCap;  // Changed to array to support multi-level caches
    protected ReplacementStrategy replacestrategy;

    public CacheTask(LayeredNetwork model, String name) {
        super(model, name);
        this.items = 1;
        this.itemLevelCap = new int[]{1};
        this.replacestrategy = ReplacementStrategy.FIFO; // Default replacement strategy
    }

    // Backward compatible constructor with single int capacity
    public CacheTask(LayeredNetwork model, String name, int items, int itemLevelCap, ReplacementStrategy replacestrategy) {
        super(model, name);
        this.items = items;
        this.itemLevelCap = new int[]{itemLevelCap};  // Convert single value to array
        this.replacestrategy = replacestrategy;
    }

    // New constructor supporting multi-level cache with int[] capacity
    public CacheTask(LayeredNetwork model, String name, int items, int[] itemLevelCap, ReplacementStrategy replacestrategy) {
        super(model, name);
        this.items = items;
        this.itemLevelCap = itemLevelCap;
        this.replacestrategy = replacestrategy;
    }

    // Backward compatible constructor with single int capacity
    public CacheTask(LayeredNetwork model, String name, int items, int itemLevelCap, ReplacementStrategy replacestrategy,
                     int multiplicity, SchedStrategy scheduling) {
        super(model, name, multiplicity, scheduling);
        this.items = items;
        this.itemLevelCap = new int[]{itemLevelCap};  // Convert single value to array
        this.replacestrategy = replacestrategy;
    }

    // New constructor supporting multi-level cache with int[] capacity
    public CacheTask(LayeredNetwork model, String name, int items, int[] itemLevelCap, ReplacementStrategy replacestrategy,
                     int multiplicity, SchedStrategy scheduling) {
        super(model, name, multiplicity, scheduling);
        this.items = items;
        this.itemLevelCap = itemLevelCap;
        this.replacestrategy = replacestrategy;
    }

    // Backward compatible constructor with single int capacity
    public CacheTask(LayeredNetwork model, String name, int items, int itemLevelCap, ReplacementStrategy replacestrategy,
                     int multiplicity, SchedStrategy scheduling, Distribution thinkTime) {
        super(model, name, multiplicity, scheduling, thinkTime);
        this.items = items;
        this.itemLevelCap = new int[]{itemLevelCap};  // Convert single value to array
        this.replacestrategy = replacestrategy;
    }

    // New constructor supporting multi-level cache with int[] capacity
    public CacheTask(LayeredNetwork model, String name, int items, int[] itemLevelCap, ReplacementStrategy replacestrategy,
                     int multiplicity, SchedStrategy scheduling, Distribution thinkTime) {
        super(model, name, multiplicity, scheduling, thinkTime);
        this.items = items;
        this.itemLevelCap = itemLevelCap;
        this.replacestrategy = replacestrategy;
    }

    // Backward compatible constructor for Python/user convenience
    public CacheTask(LayeredNetwork model, String name, int items, int itemLevelCap,
                     ReplacementStrategy replacementStrategy, int multiplicity) {
        super(model, name, multiplicity, SchedStrategy.FCFS);
        this.items = items;
        this.itemLevelCap = new int[]{itemLevelCap};  // Convert single value to array
        this.replacestrategy = replacementStrategy;
    }

    // New constructor for Python/user convenience with multi-level support
    public CacheTask(LayeredNetwork model, String name, int items, int[] itemLevelCap,
                     ReplacementStrategy replacementStrategy, int multiplicity) {
        super(model, name, multiplicity, SchedStrategy.FCFS);
        this.items = items;
        this.itemLevelCap = itemLevelCap;
        this.replacestrategy = replacementStrategy;
    }

    public int getItems() {
        return items;
    }

    public void setItems(int items) {
        this.items = items;
    }

    public int[] getItemLevelCap() {
        return itemLevelCap;
    }

    // Backward compatible getter returning total capacity (sum of all levels)
    public int getTotalCapacity() {
        int total = 0;
        for (int cap : itemLevelCap) {
            total += cap;
        }
        return total;
    }

    public void setItemLevelCap(int[] itemLevelCap) {
        this.itemLevelCap = itemLevelCap;
    }

    // Backward compatible setter with single int
    public void setItemLevelCap(int itemLevelCap) {
        this.itemLevelCap = new int[]{itemLevelCap};
    }

    public ReplacementStrategy getReplacestrategy() {
        return replacestrategy;
    }

    public void setReplacestrategy(ReplacementStrategy replacestrategy) {
        this.replacestrategy = replacestrategy;
    }
}
