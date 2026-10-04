package jline.api.sn;

import jline.GlobalConstants;
import jline.lang.NetworkStruct;

/**
 * Detects routing that is not identical across job classes.
 *
 * <p>Returns true when the routing probabilities differ between classes, either
 * because a class switches class on a hop or because two classes leave the same
 * station with different probabilities. False for a single-class model and for a
 * multiclass model in which every class traverses the network identically.
 *
 * <p>This is the condition under which per-class visit ratios diverge, so a
 * method that aggregates classes into a per-chain demand vector stops being
 * exact. Used to gate Marie's aggregation-decomposition in SolverMVA.
 *
 * <p>{@code sn.rt} is indexed station-major, {@code (i-1)*K+r}. Ported from
 * matlab/src/api/sn/sn_has_classdep_routing.m.
 *
 * @since LINE 3.0
 */
public final class SnHasClassdepRouting {
    private SnHasClassdepRouting() {}

    /**
     * Checks whether the classes of the network are routed differently.
     *
     * @param sn the NetworkStruct object for the queueing network model
     * @return true if routing is class-dependent, false if all classes route alike
     */
    public static boolean snHasClassdepRouting(NetworkStruct sn) {
        int K = sn.nclasses;
        int M = sn.nstations;
        if (K <= 1) {
            return false;
        }
        double tol = GlobalConstants.FineTol;
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < M; j++) {
                boolean haveShared = false;
                double shared = 0.0;
                for (int r = 0; r < K; r++) {
                    // Class switching makes the routing class-dependent outright.
                    for (int s = 0; s < K; s++) {
                        if (r != s && sn.rt.get(i * K + r, j * K + s) > tol) {
                            return true;
                        }
                    }
                    double p = sn.rt.get(i * K + r, j * K + r);
                    if (!haveShared) {
                        shared = p;
                        haveShared = true;
                    } else if (Math.abs(p - shared) > tol) {
                        return true;
                    }
                }
            }
        }
        return false;
    }
}
