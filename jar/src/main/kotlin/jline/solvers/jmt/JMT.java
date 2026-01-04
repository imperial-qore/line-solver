/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.jmt;

import jline.lang.Network;
import jline.solvers.SolverOptions;

/**
 * JMT is an alias for SolverJMT (Java Modelling Tools solver).
 */
public class JMT extends SolverJMT {

    public JMT(Network model, String method) {
        super(model, method);
    }

    public JMT(Network model, SolverOptions options) {
        super(model, options);
    }

    public JMT(Network model, Object... varargin) {
        super(model, varargin);
    }

    public JMT(Network model) {
        super(model);
    }
}
