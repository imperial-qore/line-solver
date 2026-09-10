/*
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
package jline.lib.fjcodes;

import java.util.Objects;

import jline.util.matrix.Matrix;

/**
 * Arrival representation for FJ_codes.
 */
public final class FJArrival {
    public final double lambda;     // Arrival rate
    public final Matrix lambda0;    // D0 matrix (transitions without arrivals)
    public final Matrix lambda1;    // D1 matrix (transitions with arrivals)
    private final int ArrChoice;     // Arrival type: 1=Exponential, other=MAP

    public FJArrival(double lambda, Matrix lambda0, Matrix lambda1, int ArrChoice) {
        this.lambda = lambda;
        this.lambda0 = lambda0;
        this.lambda1 = lambda1;
        this.ArrChoice = ArrChoice;
    }

    public double getLambda() { return lambda; }
    public Matrix getLambda0() { return lambda0; }
    public Matrix getLambda1() { return lambda1; }
    public int getArrChoice() { return ArrChoice; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof FJArrival)) return false;
        FJArrival that = (FJArrival) o;
        return Double.compare(that.lambda, lambda) == 0
                && ArrChoice == that.ArrChoice
                && Objects.equals(lambda0, that.lambda0)
                && Objects.equals(lambda1, that.lambda1);
    }

    @Override
    public int hashCode() {
        return Objects.hash(lambda, lambda0, lambda1, ArrChoice);
    }

    @Override
    public String toString() {
        return "FJArrival(lambda=" + lambda + ", lambda0=" + lambda0
                + ", lambda1=" + lambda1 + ", ArrChoice=" + ArrChoice + ")";
    }
}
