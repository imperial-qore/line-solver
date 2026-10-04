/*
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
package jline.inference.api;

/**
 * Inference problem handed to {@link Infer_variational}.
 *
 * <p>Station-class pairs are flattened column-major, so that pair (m,r) sits at
 * index r*M+m, matching the MATLAB, Python and C++ specifications. Station and
 * class indices inside {@code arcs} are one-based; index 0 marks the external
 * source or sink.</p>
 */
public class VariationalSpec {
    /** (T x 3) transitions [i j c]; i==0 external source, j==0 sink. */
    public int[][] arcs;
    /** (M x R) initial queue lengths. */
    public double[][] x0;
    /** (M) discipline codes: 0 = infinite server, 1 = shared server, 2 = external. */
    public int[] sched;
    /** (M) number of servers. */
    public double[] nservers;
    /** (T) routing probability of each transition. */
    public double[] routeprob;
    /** (T) index in 1..P of the rate governing the transition, 0 when known. */
    public int[] arcparam;
    /** (T) known rate for transitions with arcparam == 0. */
    public double[] arcrate;
    /** (P) Gamma prior shapes. */
    public double[] alpha0;
    /** (P) Gamma prior rates. */
    public double[] beta0;
    /** (K) observation epochs. */
    public double[] obsTimes;
    /** (K x M*R) observed queue lengths, NaN where not observed. */
    public double[][] obsData;
    /** (M*R) support size of the uniform contamination. */
    public double[] obsRange;
    /** probability that a reading is faulty. */
    public double epsilon;
    /**
     * (M*R) upper bound on the queue length, infinite by default. In a closed
     * network this is the chain population, and clamping the load there keeps
     * the expanded state space from crediting a station with more jobs than the
     * network holds.
     */
    public double[] capacity;

    public VariationalSpec() {}

    public int nstations() {
        return x0.length;
    }

    public int nclasses() {
        return x0[0].length;
    }

    public int narcs() {
        return arcs.length;
    }

    public int nparams() {
        return alpha0.length;
    }

    /** Flattened (M*R) copy of a station-class array, column-major. */
    public static double[] flatten(double[][] a) {
        int m = a.length;
        int r = a[0].length;
        double[] out = new double[m * r];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < r; j++) {
                out[j * m + i] = a[i][j];
            }
        }
        return out;
    }
}
