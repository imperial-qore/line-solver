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

    /**
     * "Trace", not the inherited "Replayer".
     *
     * getName() stays "Replayer" because it resolves the ProcessType and the
     * JSON wire type, and a Trace replays samples exactly as a Replayer does.
     * The registry name is what makes the Trace entry reachable; a solver that
     * declares only "Replayer" still accepts the model, through the
     * {@link jline.lang.FeatureSet} generalization table.
     *
     * @return the FeatureSet entry naming this distribution
     */
    @Override
    public String getFeatureName() {
        return "Trace";
    }
}