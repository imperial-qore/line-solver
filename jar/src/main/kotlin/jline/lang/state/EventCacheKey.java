/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.state;

import jline.lang.constant.EventType;
import jline.util.matrix.Matrix;

import java.util.Objects;

/**
 * A data structure acting as a key to the EventCache
 */
public class EventCacheKey {

    // used as a key into the afterEvent cache
    public final int ind;
    public final Matrix inspace;
    public final EventType event;
    public final int jobClass;
    public final boolean isSimulation;


    public EventCacheKey(int ind, Matrix inspace, EventType event, int jobClass, boolean isSimulation) {
        this.ind = ind;
        this.inspace = inspace;
        this.event = event;
        this.jobClass = jobClass;
        this.isSimulation = isSimulation;
    }


    @Override
    public boolean equals(Object obj) {
        if (this == obj) {
            return true;
        }
        if (obj == null || getClass() != obj.getClass()) {
            return false;
        }
        EventCacheKey that = (EventCacheKey) obj;
        return ind == that.ind &&
                jobClass == that.jobClass &&
                isSimulation == that.isSimulation &&
                inspace.equals(that.inspace) &&
                event == that.event;
    }

    @Override
    public int hashCode() {
        return Objects.hash(ind, inspace, event, jobClass, isSimulation);
    }


}
