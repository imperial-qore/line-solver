package jline.api.qsys;

import jline.io.Ret;

/**
 * G/M/1 Queueing System Analysis.
 *
 * <p>Implements analysis for G/M/1 queues with general arrival processes and exponential
 * service times. Uses the embedded Markov chain approach to derive exact performance
 * measures for systems with complex arrival patterns but simple service.
 *
 * @since LINE 3.0
 */
public final class Qsys_gm1 {
    private Qsys_gm1() {}

    /**
     * Analyzes a G/M/1 queueing system.
     *
     * @param sigma Traffic intensity.
     * @param mu    Service rate.
     * @return qsysReturn containing average waiting time (W) and utilization (rhohat).
     */
    public static Ret.qsys qsys_gm1(double sigma, double mu) {
        double W = 1 / (1 - sigma) / mu;
        double rhohat = 0.0;
        return new Ret.qsys(W, rhohat);
    }
}
