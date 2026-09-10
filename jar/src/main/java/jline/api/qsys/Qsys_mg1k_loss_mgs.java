package jline.api.qsys;

import java.util.HashMap;

public final class Qsys_mg1k_loss_mgs {
    private Qsys_mg1k_loss_mgs() {}

    /**
     * M/G/1/K loss probability using MacGregor Smith approximation.
     *
     * Calculates the loss probability for an M/G/1/K queue using the
     * MacGregor Smith approximation method.
     * Reference: J. MacGregor Smith - Optimal Design and Performance Modelling
     * of M/G/1/K Queueing Systems
     *
     * @param lambda arrival rate
     * @param mu mean service rate (1/mean service time)
     * @param mu_scv squared coefficient of variation of service time
     * @param K system capacity (maximum number of customers)
     * @return HashMap containing lossprob (loss probability) and rho (utilization)
     */
    public static HashMap<String, Object> qsys_mg1k_loss_mgs(double lambda, double mu, double mu_scv, int K) {
        HashMap<String, Object> result = new HashMap<String, Object>();

        double rho = lambda / mu;
        double s = Math.sqrt(mu_scv);
        double sqrt_rho = Math.sqrt(rho);

        double exponent1 = (sqrt_rho * s * s - sqrt_rho + 2.0 * K) / (2.0 + sqrt_rho * s * s - sqrt_rho);
        double exponent2 = 2.0 * (1.0 + sqrt_rho * s * s - sqrt_rho + K) / (2.0 + sqrt_rho * s * s - sqrt_rho);

        double lossprob_num = Math.pow(rho, exponent1) * (rho - 1.0);
        double lossprob_den = Math.pow(rho, exponent2) - 1.0;

        double lossprob = lossprob_num / lossprob_den;

        result.put("lossprob", lossprob);
        result.put("rho", rho);

        return result;
    }
}
