package jline.api.qsys;

import jline.io.Ret;

public final class Qsys_gig1_ubnd_kingman {
    private Qsys_gig1_ubnd_kingman() {}

    /**
     * Kingman's upper bound on the mean waiting time of a G/G/1 queue:
     * Wq &lt;= lambda*(sa^2+ss^2)/(2*(1-rho)), with sa^2=ca^2/lambda^2 and
     * ss^2=cs^2/mu^2. The returned W adds the mean service time, so it
     * upper-bounds the mean response time (time in system).
     *
     * <p>Reference: Kingman, J.F.C. (1962). Some inequalities for the queue
     * GI/G/1. Biometrika 49(3/4), 315-324.
     *
     * @param lambda Arrival rate.
     * @param mu     Service rate.
     * @param ca     Coefficient of variation of the arrival process.
     * @param cs     Coefficient of variation of the service time.
     * @return qsysReturn containing upper bound on average response time (W) and modified utilization (rhohat).
     */
    public static Ret.qsys qsys_gig1_ubnd_kingman(double lambda, double mu, double ca, double cs) {
        double rho = lambda / mu;
        double Wq = lambda * (ca * ca / (lambda * lambda) + cs * cs / (mu * mu)) / (2 * (1 - rho));
        double W = Wq + 1 / mu;
        double rhohat = W * lambda / (1 + W * lambda);
        return new Ret.qsys(W, rhohat);
    }
}
