package jline.api.qsys;

import org.apache.commons.math3.util.FastMath;

import jline.io.Ret;

public final class Qsys_gig1_approx_klb {
    private Qsys_gig1_approx_klb() {}

    /**
     * Analyzes a G/G/1 queueing system using the Kramer-Langenbach-Belz (KLB) approximation.
     *
     * @param lambda Arrival rate.
     * @param mu     Service rate.
     * @param ca     Coefficient of variation of the arrival process.
     * @param cs     Coefficient of variation of the service time.
     * @return qsysReturn containing average waiting time (W) and modified utilization (rhohat).
     */
    public static Ret.qsys qsys_gig1_approx_klb(double lambda, double mu, double ca, double cs) {
        // Kramer-Langenbach-Belz formula
        double rho = lambda / mu;
        double g;
        if (ca <= 1) {
            g = FastMath.exp(-2 * (1 - rho) * FastMath.pow(1 - ca * ca, 2) / (3 * rho * (ca * ca + cs * cs)));
        } else {
            g = FastMath.exp(-(1 - rho) * (ca * ca - 1) / (ca * ca + 4 * cs * cs));
        }
        double W = 1.0 / mu * ((rho / (1 - rho)) * ((cs * cs + ca * ca) / 2) * g + 1);
        double rhohat = W * lambda / (1 + W * lambda);
        return new Ret.qsys(W, rhohat);
    }
}
