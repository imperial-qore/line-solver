/**
 * @file M/GI/1 FCFS Age of Information analysis
 *
 * @since LINE 3.0
 */
package jline.api.aoi;

public final class Aoi_fcfs_mgi1 {
    private Aoi_fcfs_mgi1() {}

    /**
     * Mean AoI and LST for M/GI/1 FCFS queue.
     *
     * <p>Exact mean: E[A] = E[H] + E[T] + (1-2*rho)/lambda - dT*(s)/ds at
     * s=lambda, with T*(s) = H*(s)*W*(s) the Pollaczek-Khinchine system-time
     * LST; E[Apeak] = E[T] + E[Y].</p>
     *
     * @param lambda Arrival rate (Poisson arrivals), must be positive
     * @param H_lst  LST of service time distribution
     * @param E_H    Mean service time (first moment), must be positive
     * @param E_H2   Second moment of service time, must be >= E_H^2
     * @return AoiLstResult containing meanAoI, lstAoI, peakAoI
     * @throws IllegalArgumentException if parameters are invalid or system is unstable
     */
    public static AoiLstResult aoi_fcfs_mgi1(final double lambda, final LstFunction H_lst,
                                             double E_H, double E_H2) {
        if (!(lambda > 0)) {
            throw new IllegalArgumentException("Arrival rate lambda must be positive");
        }
        if (!(E_H > 0)) {
            throw new IllegalArgumentException("Mean service time E_H must be positive");
        }
        if (!(E_H2 >= E_H * E_H)) {
            throw new IllegalArgumentException("Second moment E_H2 must be >= E_H^2");
        }

        final double rho = lambda * E_H;
        if (!(rho < 1)) {
            throw new IllegalArgumentException(String.format(
                    "System unstable: rho = lambda*E_H = %.4f >= 1", rho));
        }

        // Mean interarrival time
        double E_Y = 1.0 / lambda;

        // Mean waiting time (Pollaczek-Khinchine formula)
        double E_W = lambda * E_H2 / (2.0 * (1.0 - rho));

        // Mean system time (sojourn time)
        double E_T = E_W + E_H;

        // see _kb/03-api-layer.md for rationale
        double hstep = 1e-6 * Math.max(1.0, lambda);
        double dTstar = (tstarMGI1(lambda + hstep, lambda, rho, H_lst)
                - tstarMGI1(lambda - hstep, lambda, rho, H_lst)) / (2.0 * hstep);
        double meanAoI = E_H + E_T + (1.0 - 2.0 * rho) / lambda - dTstar;

        // Mean Peak AoI
        double peakAoI = E_T + E_Y;

        // LST of AoI (Inoue et al. 2019, Theorem 2), via the general age formula
        //   A*(s) = (lambda/s) * ( T*(s) - Apeak*(s) )
        // the cycle average of exp(-s*age) over a departure interval: the age
        // starts each cycle at the system time T of the packet just delivered
        // and grows to the peak T + (next interarrival) at the next delivery.
        // For M/GI/1 FCFS Lindley gives W' = max(0, T - Y), so W' + Y is
        // max(Y, T), and with Y ~ Exp(lambda) independent of T,
        //   E[exp(-s*max(Y,T))] = T*(s) - (s/(s+lambda)) * T*(s+lambda).
        //
        // THE PREVIOUS FORM WAS NOT AN LST: (lambda*H*(s))/(s+lambda-lambda*H*(s))
        // diverges as s -> 0, so A*(0) was +Inf instead of 1 and the value
        // exceeded 1 for small s. Checked against simulation on M/E2/1: at
        // s = 0.3 the old form gave 1.3062, the form below 0.56935, and the
        // sample path 0.56948.
        LstFunction lstAoI = new LstFunction() {
            @Override
            public double evaluate(double s) {
                if (Math.abs(s) < 1e-12) {
                    return 1.0; // A*(0) = 1 for any proper LST
                }
                double H_s = H_lst.evaluate(s);
                double T_s = tstarMGI1(s, lambda, rho, H_lst);
                double T_sl = tstarMGI1(s + lambda, lambda, rho, H_lst);
                double peak_s = H_s * (T_s - (s / (s + lambda)) * T_sl);
                return (lambda / s) * (T_s - peak_s);
            }
        };

        return new AoiLstResult(meanAoI, lstAoI, peakAoI);
    }

    /** System-time LST T*(s) = H*(s)*W*(s) with W* the P-K waiting LST. */
    private static double tstarMGI1(double s, double lambda, double rho, LstFunction H_lst) {
        if (Math.abs(s) < 1e-12) {
            return 1.0; // removable 0/0 at the origin
        }
        double H_s = H_lst.evaluate(s);
        return H_s * (1.0 - rho) * s / (s - lambda + lambda * H_s);
    }
}
