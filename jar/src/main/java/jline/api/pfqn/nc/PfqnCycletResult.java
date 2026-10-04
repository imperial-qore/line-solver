/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.pfqn.nc;

import java.util.List;

/** Density, distribution and moments of a passage time along an overtake-free path. */
public class PfqnCycletResult {
    /** Density on the requested grid. */
    public final double[] f;
    /** Cumulative distribution on the requested grid. */
    public final double[] F;
    /** The first {@code nmom} moments. */
    public final double[] mom;
    /** Which route ran for each path: "exact" (Theorem 2) or "lt" (Theorem 1). */
    public final List<String> method;
    /** Log of the network normalizing constant at population N-1. */
    public final double lG;

    public PfqnCycletResult(double[] f, double[] F, double[] mom, List<String> method, double lG) {
        this.f = f;
        this.F = F;
        this.mom = mom;
        this.method = method;
        this.lG = lG;
    }
}
