package jline.api.qsys;

import jline.io.Ret;

public final class Qsys_gigk_approx_kingman {
    private Qsys_gigk_approx_kingman() {}

    /**
     * Analyzes a G/G/k queueing system using Kingman's approximation.
     *
     * @param lambda Arrival rate.
     * @param mu     Service rate.
     * @param ca     Coefficient of variation of the arrival process.
     * @param cs     Coefficient of variation of the service time.
     * @param k      Number of servers.
     * @return qsysReturn containing average waiting time (W) and modified utilization (rhohat).
     */
    public static Ret.qsys qsys_gigk_approx_kingman(double lambda, double mu, double ca, double cs, int k) {
        // Note: qsys_mmk populates static fields Ret.qsys.W and Ret.qsys.rho
        Qsys_mmk.qsys_mmk(lambda, mu, k);
        double W = (ca * ca + cs * cs) / 2 * (Ret.qsys.W - 1 / mu) + 1 / mu;
        double rhohat = W * lambda / (1 + W * lambda);
        return new Ret.qsys(W, rhohat);
    }
}
