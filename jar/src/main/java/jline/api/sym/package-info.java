/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * Computer algebra backend.
 *
 * <p>Java has no computer algebra system of its own, so symbolic work is
 * delegated to SageMath running behind the line-sage-rest service. This package
 * holds the operation interface ({@link jline.api.sym.SymEngine}), its SageMath
 * client ({@link jline.api.sym.SageRestEngine}) and the backend resolution and
 * container lifecycle ({@link jline.api.sym.SymEngines}).
 *
 * <p>Expressions cross the boundary as plain ASCII infix strings and are read
 * over the exact rationals, never as floating point. The same protocol backs
 * the MATLAB (SAGE.m) and Python (line_solver.api.sym) clients.
 *
 * @since LINE 2.0
 */
package jline.api.sym;
