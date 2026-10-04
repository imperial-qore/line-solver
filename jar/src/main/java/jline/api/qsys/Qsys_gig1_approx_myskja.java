package jline.api.qsys;

import java.util.HashMap;

public final class Qsys_gig1_approx_myskja {
    private Qsys_gig1_approx_myskja() {}

    /**
     * G/G/1 queue approximation using Myskja's third-moment method:
     * Wq = rho/(2*mu*(1-rho))*((1+cs2)+(q0/qa)^(1/rho-rho)*(1/rho)*(ca2-1)),
     * exact for M/G/1 (ca2=1). The returned W adds the mean service time
     * (response time, time in system).
     *
     * @param lambda arrival rate
     * @param mu service rate
     * @param ca squared coefficient of variation of inter-arrival time
     * @param cs squared coefficient of variation of service time
     * @param q0 lowest value of the relative third moment for a given mean and SCV
     * @param qa third relative moment E[X^3]/6/E[X]^3, X=inter-arrival time r.v.
     * @return HashMap containing W (mean response time) and rhohat
     */
    public static HashMap<String, Object> qsys_gig1_approx_myskja(double lambda, double mu, double ca, double cs, double q0, double qa) {
        HashMap<String, Object> result = new HashMap<String, Object>();

        double rho = lambda / mu;

        double Wq = rho / (2.0 * mu * (1.0 - rho))
                * ((1.0 + cs) + Math.pow(q0 / qa, 1.0 / rho - rho) * (1.0 / rho) * (ca - 1.0));
        double W = Wq + 1.0 / mu;
        double rhohat = W * lambda / (1.0 + W * lambda);

        result.put("W", W);
        result.put("rhohat", rhohat);

        return result;
    }
}
