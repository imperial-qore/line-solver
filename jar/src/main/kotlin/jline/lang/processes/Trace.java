/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

/**
 * Alias for the Replayer class
 */
public class Trace extends Replayer {

    /**
     * Creates a new trace from the specified data.
     * The trace acts as a replayer for the provided data sequence.
     * 
     * @param data the trace data to be replayed
     */
    public Trace(Object data) {
        super(data);
    }
}