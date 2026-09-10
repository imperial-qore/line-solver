/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.lib.butools.queues;

import java.util.Arrays;
import java.util.Objects;

import jline.lib.butools.mam.GeneralFluidSolution;

/**
 * Result class for FluFluQueue containing ME distribution parameters
 * and computed performance measures.
 */
public final class FluFluResult {
    private final GeneralFluidSolution fluidSolution;
    private final GeneralFluidSolution sojournSolution;
    private final double[] fluidMoments;
    private final double[] sojournMoments;
    private final double lambda;
    private final double mu;

    public FluFluResult(GeneralFluidSolution fluidSolution, GeneralFluidSolution sojournSolution,
                        double[] fluidMoments, double[] sojournMoments, double lambda, double mu) {
        this.fluidSolution = fluidSolution;
        this.sojournSolution = sojournSolution;
        this.fluidMoments = fluidMoments;
        this.sojournMoments = sojournMoments;
        this.lambda = lambda;
        this.mu = mu;
    }

    public GeneralFluidSolution getFluidSolution() { return fluidSolution; }
    public GeneralFluidSolution getSojournSolution() { return sojournSolution; }
    public double[] getFluidMoments() { return fluidMoments; }
    public double[] getSojournMoments() { return sojournMoments; }
    public double getLambda() { return lambda; }
    public double getMu() { return mu; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof FluFluResult)) return false;
        FluFluResult that = (FluFluResult) o;
        return Double.compare(that.lambda, lambda) == 0
                && Double.compare(that.mu, mu) == 0
                && Objects.equals(fluidSolution, that.fluidSolution)
                && Objects.equals(sojournSolution, that.sojournSolution)
                && Arrays.equals(fluidMoments, that.fluidMoments)
                && Arrays.equals(sojournMoments, that.sojournMoments);
    }

    @Override
    public int hashCode() {
        int result = Objects.hash(fluidSolution, sojournSolution, lambda, mu);
        result = 31 * result + Arrays.hashCode(fluidMoments);
        result = 31 * result + Arrays.hashCode(sojournMoments);
        return result;
    }

    @Override
    public String toString() {
        return "FluFluResult(lambda=" + lambda + ", mu=" + mu + ")";
    }
}
