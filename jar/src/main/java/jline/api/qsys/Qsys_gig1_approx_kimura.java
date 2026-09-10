package jline.api.qsys;

import java.util.HashMap;

public final class Qsys_gig1_approx_kimura {
    private Qsys_gig1_approx_kimura() {}

    /**
     * G/G/1 queue approximation using Kimura's diffusion-interpolation method:
     * Wq = rho*(ca2+cs2)/(mu*(1-rho)*(1+ca2)), exact for M/M/1 and M/G/1.
     * The returned W adds the mean service time (response time, time in system).
     *
     * <p>Reference: Kimura, T. (1986). A two-moment approximation for the mean
     * waiting time in the GI/G/s queue. Management Science 32(6), 751-763.
     *
     * @param lambda arrival rate
     * @param mu service rate
     * @param ca coefficient of variation of the inter-arrival time
     * @param cs coefficient of variation of the service time
     * @return HashMap containing W (mean response time) and rhohat
     */
    public static HashMap<String, Object> qsys_gig1_approx_kimura(double lambda, double mu, double ca, double cs) {
        HashMap<String, Object> result = new HashMap<String, Object>();

        double rho = lambda / mu;
        // ca and cs are COEFFICIENTS of variation, as everywhere else in this
        // package (see Qsys_gig1_approx_heyman); the formula squares them.
        double ca2 = ca * ca;
        double Wq = rho * (ca2 + cs * cs) / mu / (1.0 - rho) / (1.0 + ca2);
        double W = Wq + 1.0 / mu;
        double rhohat = W * lambda / (1.0 + W * lambda);

        result.put("W", W);
        result.put("rhohat", rhohat);

        return result;
    }
}
