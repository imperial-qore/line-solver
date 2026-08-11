package jline.api.qsys;

import jline.io.Ret;

/**
 * M/M/1 queueing system analysis.
 *
 * <p>Implements exact analytical solutions for the M/M/1 queue (Poisson arrivals, exponential
 * service times, single server). This fundamental queueing model provides closed-form
 * expressions for key performance metrics used in capacity planning and system design.
 *
 * @since LINE 3.0
 */
public final class Qsys_mm1 {
    private Qsys_mm1() {}

    /**
     * Analyzes an M/M/1 queueing system.
     *
     * @param lambda Arrival rate.
     * @param mu     Service rate.
     * @return qsysReturn containing average waiting time (W) and utilization (rho).
     */
    public static Ret.qsys qsys_mm1(double lambda, double mu) {
        double rho = lambda / mu;
        double W = rho / (1 - rho) / lambda;
        return new Ret.qsys(W, rho);
    }
}
