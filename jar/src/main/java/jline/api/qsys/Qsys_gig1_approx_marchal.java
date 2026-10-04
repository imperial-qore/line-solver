package jline.api.qsys;

import jline.io.Ret;

public final class Qsys_gig1_approx_marchal {
    private Qsys_gig1_approx_marchal() {}

    /**
     * Analyzes a G/G/1 queueing system using Marchal's approximation.
     *
     * @param lambda Arrival rate.
     * @param mu     Service rate.
     * @param ca     Coefficient of variation of the arrival process.
     * @param cs     Coefficient of variation of the service time.
     * @return qsysReturn containing average waiting time (W) and modified utilization (rhohat).
     */
    public static Ret.qsys qsys_gig1_approx_marchal(double lambda, double mu, double ca, double cs) {
        double rho = lambda / mu;
        double Wmm1 = rho / (1 - rho);
        double W = Wmm1 * (1 + cs * cs) / 2 / mu * (ca + rho * rho * cs * cs) / (1 + rho * rho * cs * cs) + 1 / mu;
        double rhohat = W * lambda / (1 + W * lambda);
        return new Ret.qsys(W, rhohat);
    }
}
