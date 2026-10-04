package jline.api.qsys;

import java.util.HashMap;

public final class Qsys_mm1k_loss {
    private Qsys_mm1k_loss() {}

    /**
     * M/M/1/K loss probability calculation.
     * <p>
     * Calculates the loss probability for an M/M/1/K queue (finite capacity K).
     * Based on: Niu-Cooper, Transform-Free Analysis of M/G/1/K and Related Queues,
     * Mathematics of Operations Research Vol. 18, No. 2 (May, 1993), pp. 486-510.
     *
     * @param lambda arrival rate
     * @param mu service rate
     * @param K system capacity (maximum number of customers)
     * @return HashMap containing lossprob (loss probability) and rho (utilization)
     */
    public static HashMap<String, Object> qsys_mm1k_loss(double lambda, double mu, int K) {
        HashMap<String, Object> result = new HashMap<String, Object>();

        double rho = lambda / mu;

        double lossprob;
        if (Math.abs(rho - 1.0) < 1e-10) {
            // Special case when rho = 1
            lossprob = 1.0 / (K + 1.0);
        } else {
            // General case
            lossprob = (1.0 - rho) / (1.0 - Math.pow(rho, K + 1.0)) * Math.pow(rho, (double) K);
        }

        result.put("lossprob", lossprob);
        result.put("rho", rho);

        return result;
    }
}
