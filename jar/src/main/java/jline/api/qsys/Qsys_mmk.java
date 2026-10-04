/**
 * @file M/M/k queueing system analysis
 *
 * Implements exact analytical solutions for M/M/k queues using the Erlang-C formula.
 *
 * @since LINE 3.0
 */
package jline.api.qsys;

import jline.io.Ret;
import org.apache.commons.math3.util.FastMath;

public final class Qsys_mmk {
    private Qsys_mmk() {}

    /**
     * Analyzes an M/M/k queueing system.
     *
     * @param lambda Arrival rate.
     * @param mu     Service rate.
     * @param k      Number of servers.
     * @return qsysReturn containing average waiting time (W) and utilization (rho).
     */
    public static Ret.qsys qsys_mmk(double lambda, double mu, int k) {
        double rho = lambda / mu / k;
        double Q = rho / (1 - rho) * ErlangC(rho, k) + k * rho;
        double W = Q / lambda;
        return new Ret.qsys(W, rho);
    }

    /**
     * Calculates the probability that an arriving customer is forced to join the queue
     * (i.e., all servers are occupied) in an M/M/k system.
     *
     * @param nu Utilization.
     * @param C  The number of servers.
     * @return Probability that an arriving customer is forced to join the queue.
     */
    public static double ErlangC(double nu, int C) {
        // j! is carried as a double, as MATLAB's factorial(j) is: an int holds only
        // up to 12!, so 13! wrapped NEGATIVE and every k >= 13 read a nonsense S and
        // a nonsense k!. An M/M/18 at rho = 8/9 reported Q = 23.8323 for 20.2494.
        double S = 0.0;
        double factj;
        double factj_1 = 1.0;
        for (int j = 0; j <= C - 1; j++) {
            if (j == 0) {
                S += FastMath.pow(C * nu, j);
            } else {
                factj = j * factj_1;
                S += FastMath.pow(C * nu, j) / factj;
                factj_1 = factj;
            }
        }
        return 1.0 / (1 + (1 - nu) * (C * factj_1) / FastMath.pow(C * nu, C) * S);
    }
}
