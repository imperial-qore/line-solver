/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.lib.ode;

/**
 * Functional interface for the ODE right-hand side: dy/dt = f(t, y).
 * The function writes derivatives into the yp array.
 */
public interface OdeFunction {
    void evaluate(double t, double[] y, double[] yp);
}
