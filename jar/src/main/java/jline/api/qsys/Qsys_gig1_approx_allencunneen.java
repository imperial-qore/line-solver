/**
 * @file Allen-Cunneen approximation for G/G/1 queues
 *
 * Implements the widely-used Allen-Cunneen approximation for general G/G/1 queueing
 * systems. This two-moment approximation provides excellent accuracy for most practical
 * applications and is considered one of the best general-purpose G/G/1 approximations.
 *
 * @since LINE 3.0
 */
package jline.api.qsys;

import jline.io.Ret;

public final class Qsys_gig1_approx_allencunneen {
    private Qsys_gig1_approx_allencunneen() {}

    /**
     * Analyzes a G/G/1 queueing system using the Allen-Cunneen approximation.
     *
     * @param lambda Arrival rate.
     * @param mu     Service rate.
     * @param ca     Coefficient of variation of the arrival process.
     * @param cs     Coefficient of variation of the service time.
     * @return qsysReturn containing average waiting time (W) and modified utilization (rhohat).
     */
    public static Ret.qsys qsys_gig1_approx_allencunneen(double lambda, double mu, double ca, double cs) {
        double rho = lambda / mu;
        double W = (rho / (1 - rho)) / mu * ((cs * cs + ca * ca) / 2) + 1.0 / mu;
        double rhohat = W * lambda / (1 + W * lambda);
        return new Ret.qsys(W, rhohat);
    }
}
