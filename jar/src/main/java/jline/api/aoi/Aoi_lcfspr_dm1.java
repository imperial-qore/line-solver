/**
 * @file D/M/1 preemptive LCFS Age of Information analysis
 *
 * @since LINE 3.0
 */
package jline.api.aoi;

public final class Aoi_lcfspr_dm1 {
    private Aoi_lcfspr_dm1() {}

    /**
     * Mean, variance, and peak AoI for D/M/1 preemptive LCFS queue.
     *
     * <p>Exact: the age is the backward recurrence time U(0,tau) plus an
     * independent Exp(mu), so E[A] = tau/2 + 1/mu and Var[A] = tau^2/12 + 1/mu^2;
     * E[Apeak] = E[S given success] + tau/q with q = 1 - exp(-mu*tau).</p>
     *
     * @param tau Deterministic interarrival time, must be positive
     * @param mu Service rate (exponential service), must be positive
     * @return AoiResult containing meanAoI, varAoI, peakAoI
     * @throws IllegalArgumentException if parameters are invalid or system is unstable
     */
    public static AoiResult aoi_lcfspr_dm1(double tau, double mu) {
        if (!(tau > 0)) throw new IllegalArgumentException("Interarrival time tau must be positive");
        if (!(mu > 0)) throw new IllegalArgumentException("Service rate mu must be positive");

        double lambda = 1.0 / tau;
        double rho = lambda / mu;
        if (!(rho < 1)) {
            throw new IllegalArgumentException(String.format("System unstable: rho = 1/(tau*mu) = %.4f >= 1", rho));
        }

        // Mean AoI for preemptive LCFS: backward recurrence U(0,tau) plus Exp(mu)
        double meanAoI = tau / 2.0 + 1.0 / mu;
        // Peak AoI (exact, as MATLAB): success prob q = 1 - exp(-mu*tau)
        double q = 1.0 - Math.exp(-mu * tau);
        double esSucc = (1.0 / mu - tau * Math.exp(-mu * tau) - Math.exp(-mu * tau) / mu) / q;
        double peakAoI = esSucc + tau / q;

        // Var[A] = Var[U(0,tau)] + Var[Exp(mu)]
        double varAoI = tau * tau / 12.0 + 1.0 / (mu * mu);

        return new AoiResult(meanAoI, varAoI, peakAoI);
    }
}
