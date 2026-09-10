package jline.api.qsys;

import java.util.HashMap;

import jline.io.Ret;

public final class Qsys_gig1_approx_myskja2 {
    private Qsys_gig1_approx_myskja2() {}

    /**
     * G/G/1 queue approximation using Myskja's enhanced third-moment method,
     * returning the mean response time (time in system). For ca2=1 the
     * interpolation parameter theta is a 0/0 form, so the exact M/G/1 result
     * is returned instead (also the interpolation anchor of the method).
     *
     * @param lambda arrival rate
     * @param mu service rate
     * @param ca squared coefficient of variation of inter-arrival time
     * @param cs squared coefficient of variation of service time
     * @param q0 lowest value of the relative third moment for a given mean and SCV
     * @param qa third relative moment E[X^3]/6/E[X]^3, X=inter-arrival time r.v.
     * @return HashMap containing W (mean response time) and rhohat
     */
    public static HashMap<String, Object> qsys_gig1_approx_myskja2(double lambda, double mu,
                                                                   double ca, double cs,
                                                                   double q0, double qa) {
        HashMap<String, Object> result = new HashMap<String, Object>();

        if (Math.abs(ca - 1.0) < 1e-8) {
            // M/G/1 case: exact (qsys_mg1 expects the coefficient of variation)
            Qsys_mg1.qsys_mg1(lambda, mu, Math.sqrt(cs));
            double W = Ret.qsys.W;
            result.put("W", W);
            result.put("rhohat", Ret.qsys.rho);
            return result;
        }

        double ra = (1.0 + ca) / 2.0;
        double rs = (1.0 + cs) / 2.0;
        double rho = lambda / mu;

        double theta = (rho * (qa - ra) - (qa - ra * ra)) / (2.0 * rho * (ra - 1.0));
        double d = (1.0 + 1.0 / ra) * (1.0 - rs) * (1.0 - Math.pow(q0 / qa, 3.0)) * (1.0 - Math.pow(rho, 3.0));
        double D = Math.pow(rs - theta, 2.0) + (2.0 * rs - 1.0 + d) * (ra - 1.0);
        D = Math.max(D, 0.0); // guard small negative values due to round-off

        double W = (rho / (1.0 - rho)) / lambda * (rs + (1.0 / rho) * (Math.sqrt(D) - (rs - theta)));
        double rhohat = W * lambda / (1.0 + W * lambda);

        result.put("W", W);
        result.put("rhohat", rhohat);

        return result;
    }
}
