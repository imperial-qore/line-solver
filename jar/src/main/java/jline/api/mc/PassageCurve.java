/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.mc;

/** A passage-time law on a time grid. */
public class PassageCurve {
    /** The time grid. */
    public final double[] t;
    /** The cumulative distribution on that grid. */
    public final double[] F;
    /** The density on that grid. */
    public final double[] f;
    /** Mass of the initial law already inside the target: F(0). */
    public final double atom;

    public PassageCurve(double[] t, double[] F, double[] f, double atom) {
        this.t = t;
        this.F = F;
        this.f = f;
        this.atom = atom;
    }
}
