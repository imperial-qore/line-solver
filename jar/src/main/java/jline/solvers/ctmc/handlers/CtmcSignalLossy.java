/**
 * Classes annihilated in place by a G-network signal at a given station.
 *
 * A job removed by a signal leaves the station without a service completion,
 * so the arrival-based utilization estimator counts offered rather than
 * carried load and only the departure-based estimator UN = T * E[S] / c is
 * meaningful for such a class. On the textbook single-queue G-network
 * (lambda+ = 0.5, lambda- = 0.4, mu = 1) the arrival estimator yields 0.5
 * against the exact rho = lambda+/(mu + lambda-) = 0.35714.
 *
 * Mirrors matlab/src/solvers/CTMC/ctmc_signal_lossy.m.
 *
 * @since LINE 3.0
 */
package jline.solvers.ctmc.handlers;

import jline.lang.NetworkStruct;

public final class CtmcSignalLossy {
    private CtmcSignalLossy() {}

    /**
     * @param sn              network structure
     * @param arvRateByClass  stationary arrival rate of each class at the station
     *                        of interest (length nclasses)
     * @return per-class flags, true when the class can be removed by a signal there
     */
    public static boolean[] signalLossyClasses(NetworkStruct sn, double[] arvRateByClass) {
        int K = arvRateByClass.length;
        boolean[] lossy = new boolean[K];
        if (sn.issignal == null || sn.issignal.isEmpty()) {
            return lossy;
        }
        boolean[] issignal = new boolean[K];
        boolean anySignal = false;
        for (int r = 0; r < K; r++) {
            issignal[r] = r < sn.issignal.length() && sn.issignal.get(r) > 0;
            anySignal = anySignal || issignal[r];
        }
        if (!anySignal) {
            return lossy;
        }
        for (int r = 0; r < K; r++) {
            if (!issignal[r] || arvRateByClass[r] <= 0) {
                continue;
            }
            // sn.signaltarget is 0-based in the JAR; -1 means untargeted.
            int tgt = -1;
            if (sn.signaltarget != null && r < sn.signaltarget.length()) {
                tgt = (int) sn.signaltarget.get(r);
            }
            if (tgt >= 0 && tgt < K) {
                lossy[tgt] = true;
            } else {
                // Untargeted signals are class-agnostic: they remove any
                // non-signal job, matching SignalRemoval, MAM and LDES.
                for (int j = 0; j < K; j++) {
                    if (!issignal[j]) {
                        lossy[j] = true;
                    }
                }
            }
        }
        return lossy;
    }
}
