/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.pfqn.ld;

import jline.util.SerializableFunction;
import jline.util.matrix.Matrix;

/**
 * Peak of a class-dependence handle over the population lattice.
 *
 * Ports MATLAB cd_peak_scaling.m. This is the single normalizer used to
 * report utilization at stations with limited class dependence when
 * options.config.cd_peak_norm is enabled, U = T*S/bmax, so that every solver
 * follows the same convention as solver_ncld does for lldscaling
 * (U/max(lldscaling)).
 */
public final class CdPeakScaling {
    private CdPeakScaling() {}

    /**
     * Warning emitted when utilization at a class-dependent station is
     * reported as the unscaled offered load (the default convention); the
     * text matches the MATLAB solvers verbatim.
     */
    public static final String CD_PEAK_NORM_WARNING =
            "Utilization at class-dependent stations is the unscaled offered load T*S; "
            + "set options.config.cd_peak_norm=true to normalize by the peak scaling (bmax).";

    /**
     * Peak of the class-dependence function over the reachable population
     * lattice 0 &lt;= n[r] &lt;= NK[r]. The function returns either a scalar
     * (shared by every class) or a length-K vector, so the peak is taken over
     * both the states and the classes: utilization is a per-station quantity,
     * so the whole station shares one normalizer, as it does for
     * max(lldscaling(ist,:)).
     *
     * @param beta class-dependence handle beta(n)
     * @param NK   per-class (finite) population vector
     * @param K    number of classes
     * @return the lattice peak of the handle (0 if the handle never returns a
     *         finite positive value)
     */
    public static double cd_peak_scaling(SerializableFunction<Matrix, Matrix> beta, int[] NK, int K) {
        double bmax = 0.0;
        int[] n = new int[K];
        while (true) {
            int tot = 0;
            for (int r = 0; r < K; r++) tot += n[r];
            if (tot > 0) {
                Matrix nv = new Matrix(1, K);
                for (int r = 0; r < K; r++) nv.set(0, r, n[r]);
                Matrix bv = beta.apply(nv);
                for (int r = 0; r < bv.getNumElements(); r++) {
                    double v = bv.get(r);
                    if (Double.isFinite(v) && v > bmax) bmax = v;
                }
            }
            if (!pprodNext(n, NK)) break;
        }
        return bmax;
    }

    private static boolean pprodNext(int[] n, int[] N) {
        int R = n.length;
        boolean atMax = true;
        for (int i = 0; i < R; i++) {
            if (n[i] != N[i]) { atMax = false; break; }
        }
        if (atMax) return false;
        int s = R - 1;
        while (s >= 0 && n[s] == N[s]) {
            n[s] = 0;
            s--;
        }
        if (s >= 0) n[s]++;
        return true;
    }
}
