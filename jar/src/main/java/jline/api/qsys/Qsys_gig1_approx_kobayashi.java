package jline.api.qsys;

import jline.io.Ret;
import org.apache.commons.math3.util.FastMath;

public final class Qsys_gig1_approx_kobayashi {
    private Qsys_gig1_approx_kobayashi() {}

    /**
     * Analyzes a G/G/1 queueing system using Kobayashi's approximation.
     *
     * @param lambda Arrival rate.
     * @param mu     Service rate.
     * @param ca     Coefficient of variation of the arrival process.
     * @param cs     Coefficient of variation of the service time.
     * @return qsysReturn containing average waiting time (W) and modified utilization (rhohat).
     */
    public static Ret.qsys qsys_gig1_approx_kobayashi(double lambda, double mu, double ca, double cs) {
        double rho = lambda / mu;
        double rhohat = FastMath.exp(-2 * (1 - rho) / (rho * (ca * ca + cs * cs / rho)));
        double W = rhohat / (1 - rhohat) / lambda;
        return new Ret.qsys(W, rhohat);
    }
}
