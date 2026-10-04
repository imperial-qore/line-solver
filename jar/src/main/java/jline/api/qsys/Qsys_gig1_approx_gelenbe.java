package jline.api.qsys;

import java.util.HashMap;

public final class Qsys_gig1_approx_gelenbe {
    private Qsys_gig1_approx_gelenbe() {}

    /**
     * G/G/1 queue approximation using Gelenbe's diffusion method with
     * instantaneous-return boundary:
     * p(0)=1-rho, p(n)=rho*(1-rhat)*rhat^(n-1) with
     * rhat=exp(-2*(1-rho)/(rho*ca2+cs2)), hence E[N]=rho/(1-rhat) and the
     * mean response time (time in system) is W=E[N]/lambda=1/(mu*(1-rhat)).
     *
     * <p>Reference: Gelenbe, E. (1975). On approximate computer system models.
     * Journal of the ACM 22(2), 261-269.
     *
     * @param lambda arrival rate
     * @param mu service rate
     * @param ca coefficient of variation of the inter-arrival time
     * @param cs coefficient of variation of the service time
     * @return HashMap containing W (mean response time) and rhohat
     */
    public static HashMap<String, Object> qsys_gig1_approx_gelenbe(double lambda, double mu, double ca, double cs) {
        HashMap<String, Object> result = new HashMap<String, Object>();

        double rho = lambda / mu;
        // ca and cs are COEFFICIENTS of variation, as everywhere else in this
        // package (see Qsys_gig1_approx_heyman); the formula squares them.
        double rhat = Math.exp(-2.0 * (1.0 - rho) / (rho * ca * ca + cs * cs));
        double W = 1.0 / (mu * (1.0 - rhat));
        double rhohat = W * lambda / (1.0 + W * lambda);

        result.put("W", W);
        result.put("rhohat", rhohat);

        return result;
    }
}
