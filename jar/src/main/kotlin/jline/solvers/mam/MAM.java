/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.mam;

import jline.lang.Network;
import jline.solvers.SolverOptions;

/**
 * MAM is an alias for SolverMAM (Matrix Analytic Methods solver).
 */
public class MAM extends SolverMAM {

    public MAM(Network model, String method) {
        super(model, method);
    }

    public MAM(Network model, SolverOptions options) {
        super(model, options);
    }

    public MAM(Network model, Object... varargin) {
        super(model, varargin);
    }

    public MAM(Network model) {
        super(model);
    }
}
