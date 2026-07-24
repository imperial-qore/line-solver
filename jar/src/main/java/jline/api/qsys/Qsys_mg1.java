/**
 * @file M/G/1 queueing system analysis
 *
 * Implements the Pollaczek-Khinchine formula for M/G/1 queues with Poisson arrivals
 * and general service time distributions. Uses the first two moments of service time
 * to compute exact performance measures via transform methods.
 *
 * @since LINE 3.0
 */
package jline.api.qsys;

import jline.io.Ret;

public final class Qsys_mg1 {
    private Qsys_mg1() {}

    /**
     * Analyzes an M/G/1 queueing system.
     *
     * @param lambda Arrival rate.
     * @param mu     Service rate.
     * @param cs     Coefficient of variation of the service time.
     * @return qsysReturn containing average waiting time (W) and modified utilization (rhohat).
     */
    public static Ret.qsys qsys_mg1(double lambda, double mu, double cs) {
        double rho = lambda / mu;
        double Q = rho + rho * rho / (2 * (1 - rho)) + lambda * lambda * cs * cs / (mu * mu) / (2 * (1 - rho));
        double W = Q / lambda;
        double rhohat = Q / (1 + Q);
        return new Ret.qsys(W, rhohat);
    }
}
