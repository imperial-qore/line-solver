package jline.api.qsys;

import java.util.HashMap;

public final class Qsys_gig1_lbnd {
    private Qsys_gig1_lbnd() {}

    /**
     * G/G/1 queue lower bounds.
     *
     * Computes fundamental theoretical lower bounds for G/G/1 queues.
     * These are the minimum possible values that performance measures
     * cannot fall below for any realization of the arrival and service processes.
     *
     * @param lambda arrival rate
     * @param mu service rate
     * @param ca2 squared coefficient of variation of inter-arrival time
     * @param cs2 squared coefficient of variation of service time
     * @return HashMap containing:
     *         - L: lower bound on average number in system
     *         - Lq: lower bound on average number in queue
     *         - W: lower bound on average time in system
     *         - Wq: lower bound on average waiting time in queue
     *         - p0: upper bound on probability of empty system
     */
    public static HashMap<String, Object> qsys_gig1_lbnd(double lambda, double mu, double ca2, double cs2) {
        HashMap<String, Object> result = new HashMap<String, Object>();

        double rho = lambda / mu;

        // Fundamental lower bounds
        double L = rho;        // At least the average number being served
        double W = 1.0 / mu;   // At least the service time
        double Lq = 0.0;       // Queue length is non-negative
        double Wq = 0.0;       // Waiting time is non-negative
        double p0 = 1.0 - rho; // Upper bound on empty system probability

        result.put("L", L);
        result.put("Lq", Lq);
        result.put("W", W);
        result.put("Wq", Wq);
        result.put("p0", p0);

        return result;
    }
}
