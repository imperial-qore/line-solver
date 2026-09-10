/**
 * @file Canonical RBM correlation weight function for the RQNA
 *
 * Canonical reflected-Brownian-motion correlation weight w*(t) used by the
 * Robust Queueing Network Analyzer (W. Whitt and W. You, 2018).
 *
 * @since LINE 3.0
 */
package jline.api.npfqn;

import org.apache.commons.math3.special.Erf;
import org.apache.commons.math3.util.FastMath;

public final class Npfqn_rqna_weight {
    private Npfqn_rqna_weight() {}

    /**
     * Canonical RBM correlation weight function w*(t).
     *
     * w*(t) = 1 - (1 - c*(t))/(2 t), where c*(t) is the correlation function of
     * the stationary version of canonical reflected Brownian motion (drift -1,
     * diffusion coefficient 1),
     *
     *   c*(t) = 2(1 - 2t - t^2) Phi^c(sqrt(t)) + 2 sqrt(t) phi(sqrt(t)) (1 + t),
     *
     * with Phi^c the standard-normal complementary cdf and phi its density. The
     * weight is monotonically increasing with w*(0)=0 and w*(Inf)=1.
     *
     * Reference: Whitt and You (2018), eqs. (24)-(25).
     *
     * @param t nonnegative time argument
     * @return the weight w*(t) in [0,1]
     */
    public static double npfqn_rqna_weight(double t) {
        if (t <= 0) {
            return 0.0;
        }
        if (Double.isInfinite(t)) {
            return 1.0;
        }
        double st = FastMath.sqrt(t);
        double phic = 0.5 * Erf.erfc(st / FastMath.sqrt(2.0));      // complementary std-normal cdf
        double phi = FastMath.exp(-t / 2.0) / FastMath.sqrt(2.0 * Math.PI);  // std-normal density at sqrt(t)
        double cstar = 2.0 * (1.0 - 2.0 * t - t * t) * phic + 2.0 * st * phi * (1.0 + t);
        double w;
        if (t < 1e-6) {
            // limit w*(t) -> 0 as t -> 0; series avoids catastrophic cancellation
            w = 0.0;
        } else {
            w = 1.0 - (1.0 - cstar) / (2.0 * t);
        }
        if (w < 0.0) {
            w = 0.0;
        }
        if (w > 1.0) {
            w = 1.0;
        }
        return w;
    }

    /**
     * Vector form of the RBM correlation weight function.
     *
     * @param t array of nonnegative time arguments
     * @return array of weights w*(t), same length as t
     */
    public static double[] npfqn_rqna_weight(double[] t) {
        double[] w = new double[t.length];
        for (int i = 0; i < t.length; i++) {
            w[i] = npfqn_rqna_weight(t[i]);
        }
        return w;
    }
}
