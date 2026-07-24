/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang;

/**
 * An environment model defined by a collection of network sub-models coupled with an environment transition rule
 * that selects the active sub-model.
 * @deprecated Use {@link Environment} instead.
 */
@Deprecated
public class Env extends Environment {
    public Env(String name, int numStages) {
        super(name, numStages);
    }
}
