package jline.api.qsys;

import jline.io.Ret;

public final class Qsys_gig1_approx_heyman {
    private Qsys_gig1_approx_heyman() {}

    /**
     * Analyzes a G/G/1 queueing system using Heyman's approximation.
     *
     * @param lambda Arrival rate.
     * @param mu     Service rate.
     * @param ca     Coefficient of variation of the arrival process.
     * @param cs     Coefficient of variation of the service time.
     * @return qsysReturn containing average waiting time (W) and modified utilization (rhohat).
     */
    public static Ret.qsys qsys_gig1_approx_heyman(double lambda, double mu, double ca, double cs) {
        double rho = lambda / mu;
        double W = rho / (1 - rho) / mu * (ca * ca + cs * cs) / 2 + 1.0 / mu;
        double rhohat = W * lambda / (1 + W * lambda);
        return new Ret.qsys(W, rhohat);
    }
}
