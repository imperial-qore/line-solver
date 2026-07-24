package jline.api.qsys;

import org.apache.commons.math3.util.FastMath;

import jline.io.Ret;

public final class Qsys_gigk_approx {
    private Qsys_gigk_approx() {}

    /**
     * Analyzes a G/G/k queueing system using an approximation method.
     *
     * @param lambda Arrival rate.
     * @param mu     Service rate.
     * @param ca     Coefficient of variation of the arrival process.
     * @param cs     Coefficient of variation of the service time.
     * @param k      Number of servers.
     * @return qsysReturn containing average waiting time (W) and modified utilization (rhohat).
     */
    public static Ret.qsys qsys_gigk_approx(double lambda, double mu, double ca, double cs, int k) {
        double rho = lambda / (mu * k);
        double alpha;
        if (rho > 0.7) {
            alpha = (Math.pow(rho, (double) k) + rho) / 2;
        } else {
            alpha = FastMath.pow(rho, (k + 1) / 2.0);
        }
        double W = (alpha / mu) * (1 / (1 - rho)) * (ca * ca + cs * cs) / (2 * k) + 1.0 / mu;
        double rhohat = W * lambda / (1 + W * lambda);
        return new Ret.qsys(W, rhohat);
    }
}
