/**
 * @file M/D/1 preemptive LCFS Age of Information analysis
 *
 * @since LINE 3.0
 */
package jline.api.aoi;

public final class Aoi_lcfspr_md1 {
    private Aoi_lcfspr_md1() {}

    /**
     * Mean, variance, and peak AoI for M/D/1 preemptive LCFS queue.
     *
     * <p>Exact: E[A] = exp(lambda*d)/lambda, Var[A] = (exp(2*lambda*d)
     * - 2*lambda*d*exp(lambda*d))/lambda^2 and E[Apeak] = d + exp(lambda*d)/lambda,
     * since deliveries occur at rate lambda*exp(-lambda*d).</p>
     *
     * @param lambda Arrival rate (Poisson arrivals), must be positive
     * @param d Deterministic service time, must be positive
     * @return AoiResult containing meanAoI, varAoI, peakAoI
     * @throws IllegalArgumentException if parameters are invalid or system is unstable
     */
    public static AoiResult aoi_lcfspr_md1(double lambda, double d) {
        if (!(lambda > 0)) throw new IllegalArgumentException("Arrival rate lambda must be positive");
        if (!(d > 0)) throw new IllegalArgumentException("Service time d must be positive");

        double rho = lambda * d;
        if (!(rho < 1)) {
            throw new IllegalArgumentException(String.format("System unstable: rho = lambda*d = %.4f >= 1", rho));
        }

        // Mean AoI for preemptive LCFS: the delivery rate is lambda*exp(-lambda*d)
        double meanAoI = Math.exp(lambda * d) / lambda;
        // Peak AoI (exact, as MATLAB): deliveries at rate lambda*exp(-lambda*d)
        double peakAoI = d + Math.exp(lambda * d) / lambda;

        // Var[A] = E[A^2] - E[A]^2 from A*(s) = f(s)/(s + f(s)), f(s) = lambda*exp(-(s+lambda)*d)
        double varAoI = (Math.exp(2.0 * lambda * d) - 2.0 * lambda * d * Math.exp(lambda * d)) / (lambda * lambda);

        return new AoiResult(meanAoI, varAoI, peakAoI);
    }
}
