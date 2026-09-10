/**
 * @file Service variability parameter of the Robust Queueing Theory framework
 *
 * Adaptation of the RQT uncertainty sets to the first two moments of a
 * stochastic queue, per C. Bandi, D. Bertsimas, N. Youssef (2015), "Robust
 * Queueing Theory", Operations Research 63(3), 676-700, Section 7.1, Table 1.
 *
 * @since LINE 3.0
 */
package jline.api.qsys;

import org.apache.commons.math3.util.FastMath;

public final class Qsys_gigk_rqt_gamma {
    private Qsys_gigk_rqt_gamma() {}

    /**
     * Service variability parameter of the RQT framework, obtained from the
     * first two moments by the adaptation of Section 7.1,
     * Gamma_s = (2 (theta0 + theta1 sigma_s^2/k + theta2 Gamma_a^2 rho^2 k))^((a-1)/a)
     * - Gamma_a k^((a-1)/a), where (theta0,theta1,theta2) are regressed so that
     * the worst-case system time of Theorem 3 approximates the MEAN system time
     * of the corresponding stochastic queue. The arrival side needs no
     * adaptation: Gamma_a = sigma_a for an external renewal stream. Since the
     * last term cancels Gamma_a at alpha=2, the adaptation acts on the sum
     * Gamma_a + Gamma_s/k^(1/alpha) that Theorem 3 reads.
     *
     * <p>THE FACTOR 2 IS NOT IN THE PRINTED FORMULA and is restored here.
     * Section 7.1 states that the functional form is motivated by Kingman's
     * bound, which the alpha=2 bound of Theorem 3 reproduces when
     * (Gamma_a+Gamma_s)^2 = 2(sigma_a^2+sigma_s^2); the published thetas are all
     * near unity, i.e. corrections to that bound rather than a substitute for
     * its factor 2. Dropping the factor puts M/M/1 about 40% BELOW its exact
     * mean system time at rho=0.9, contradicting the errors of at most 9.5% that
     * Tables 2-3 report; restoring it gives +4.7%.
     *
     * <p>CAUTION: the form is not dimensionally homogeneous, since theta0 is an
     * additive constant on a scale of variances, so it is only valid in the time
     * unit the regression was run in. It is evaluated here in units of the mean
     * service time, 1/mu = 1, and converted back.
     *
     * @param rho     traffic intensity lambda/(k*mu)
     * @param mu      service rate of each server, which sets the time unit
     * @param Gamma_a variability parameter of the arrival uncertainty set
     * @param sigma_s standard deviation of the service time
     * @param k       number of servers
     * @param alpha_a effective arrival tail coefficient in (1,2]
     * @param regime  adaptation regime of Table 1: "independent" (service
     *                distribution unknown), "normal" or "pareto"
     * @return the service variability parameter Gamma_s
     */
    public static double qsys_gigk_rqt_gamma(double rho, double mu, double Gamma_a, double sigma_s,
                                             int k, double alpha_a, String regime) {
        double t0;
        double t1;
        double t2;
        String reg = regime == null ? "independent" : regime.toLowerCase();
        if ("pareto".equals(reg)) {
            t0 = -0.05; t1 = 1.09; t2 = 1.11;
        } else if ("normal".equals(reg)) {
            t0 = -0.02; t1 = 1.03; t2 = 1.04;
        } else if ("independent".equals(reg) || "default".equals(reg)) {
            t0 = -0.06; t1 = 1.07; t2 = 1.07;
        } else {
            throw new RuntimeException("Unknown RQT adaptation regime: " + regime);
        }

        // evaluate in units of the mean service time, then convert back
        double ga = Gamma_a * mu;
        double ss = sigma_s * mu;
        double e = (alpha_a - 1) / alpha_a;
        double b = 2 * (t0 + t1 * ss * ss / k + t2 * ga * ga * rho * rho * k);
        double gs = FastMath.pow(FastMath.max(b, 0.0), e) - ga * FastMath.pow(k, e);
        return gs / mu;
    }

    /**
     * Service variability parameter under the service-distribution-independent
     * regime and the finite-variance tail coefficient alpha = 2.
     *
     * @param rho     traffic intensity lambda/(k*mu)
     * @param mu      service rate of each server
     * @param Gamma_a variability parameter of the arrival uncertainty set
     * @param sigma_s standard deviation of the service time
     * @param k       number of servers
     * @return the service variability parameter Gamma_s
     */
    public static double qsys_gigk_rqt_gamma(double rho, double mu, double Gamma_a, double sigma_s, int k) {
        return qsys_gigk_rqt_gamma(rho, mu, Gamma_a, sigma_s, k, 2.0, "independent");
    }
}
