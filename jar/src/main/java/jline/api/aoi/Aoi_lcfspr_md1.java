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
     * <p>Exact peak: E[Apeak] = d + exp(lambda*d)/lambda, since deliveries
     * occur at rate lambda*exp(-lambda*d).</p>
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

        double E_Y = 1.0 / lambda;
        double E_S = d;

        double meanAoI = E_Y + E_S;
        // Peak AoI (exact, as MATLAB): deliveries at rate lambda*exp(-lambda*d)
        double peakAoI = d + Math.exp(lambda * d) / lambda;

        // Var[A] = Var[Y] + Var[S] = 1/lambda^2 + 0
        double varAoI = 1.0 / (lambda * lambda);

        return new AoiResult(meanAoI, varAoI, peakAoI);
    }
}
