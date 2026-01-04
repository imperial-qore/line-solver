/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.qns;

import jline.lang.Network;
import jline.solvers.SolverOptions;

/**
 * QNS is an alias for SolverQNS (Queueing Network Solver).
 */
public class QNS extends SolverQNS {

    public QNS(Network model, String method) {
        super(model, method);
    }

    public QNS(Network model, SolverOptions options) {
        super(model, options);
    }

    public QNS(Network model, Object... varargin) {
        super(model, varargin);
    }

    public QNS(Network model) {
        super(model);
    }
}
