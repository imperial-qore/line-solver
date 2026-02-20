/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.mam;

import jline.lang.constant.SolverType;
import jline.solvers.SolverOptions;

public class MAMOptions extends SolverOptions {
    public int maxStates = 100;

    public MAMOptions() {
        super(SolverType.MAM);
    }

    public MAMOptions maxStates(int maxStates) {
        this.maxStates = maxStates;
        return this;
    }
}
