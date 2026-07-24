/**
 * @file M/M/c/K queueing system analysis
 *
 * Exact closed-form analysis of the M/M/c/K queue (c servers, total system
 * capacity K, blocked arrivals lost). Port of MATLAB qsys_mmck.m and of the
 * python-native qsys_mmck.
 *
 * @since LINE 3.0
 */
package jline.api.qsys;

import org.apache.commons.math3.util.FastMath;

public final class Qsys_mmck {
    private Qsys_mmck() {}

    /** Result container for the M/M/c/K metrics. */
    public static class Result {
        /** Mean number of jobs in the system, L. */
        public final double meanQueueLength;
        /** Mean number waiting in queue, Lq. */
        public final double meanQueueLengthQ;
        /** Mean waiting time in queue, Wq (Little's law on the effective rate). */
        public final double meanWaitingTime;
        /** Mean sojourn time, W. */
        public final double meanSojournTime;
        /** Per-server utilization, lambdaEff/(c*mu). */
        public final double utilization;
        /** Effective throughput, lambda*(1 - pK). */
        public final double throughput;
        /** Blocking probability, pK. */
        public final double lossProbability;

        Result(double L, double Lq, double Wq, double W, double util, double tput, double pK) {
            this.meanQueueLength = L;
            this.meanQueueLengthQ = Lq;
            this.meanWaitingTime = Wq;
            this.meanSojournTime = W;
            this.utilization = util;
            this.throughput = tput;
            this.lossProbability = pK;
        }
    }

    /**
     * Analyzes an M/M/c/K queueing system exactly.
     *
     * Stationary distribution (truncated Erlang form): with a = lambda/mu and
     * rho = a/c, p_n prop a^n/n! for n &lt;= c and a^c/c! * rho^(n-c) for
     * c &lt;= n &lt;= K, normalized over the K+1 levels.
     *
     * @param lambda Poisson arrival rate (&gt; 0)
     * @param mu     Per-server exponential service rate (&gt; 0)
     * @param c      Number of servers (&gt;= 1)
     * @param K      System capacity, total jobs allowed (K &gt;= c)
     * @return Result with mean counts, times, utilization, throughput and loss
     */
    public static Result qsys_mmck(double lambda, double mu, int c, int K) {
        if (lambda <= 0 || mu <= 0 || c < 1 || K < c) {
            throw new IllegalArgumentException("qsys_mmck: require lambda>0, mu>0, c>=1, K>=c");
        }
        double a = lambda / mu;
        double rho = a / c;

        double[] p = new double[K + 1];
        double fact = 1.0; // n! running product
        for (int n = 0; n < c; n++) {
            if (n > 0) fact *= n;   // fact == n! after this line
            p[n] = FastMath.pow(a, n) / fact;
        }
        // fact == (c-1)! on exit; c! = (c-1)! * c
        double cfact = fact * c;
        double acOverCfact = FastMath.pow(a, c) / cfact;
        for (int n = c; n <= K; n++) {
            p[n] = acOverCfact * FastMath.pow(rho, n - c);
        }
        double S = 0.0;
        for (int n = 0; n <= K; n++) S += p[n];
        for (int n = 0; n <= K; n++) p[n] /= S;

        double L = 0.0, Lq = 0.0;
        for (int n = 0; n <= K; n++) {
            L += n * p[n];
            if (n > c) Lq += (n - c) * p[n];
        }
        double pK = p[K];
        double lambdaEff = lambda * (1.0 - pK);
        double util = lambdaEff / (c * mu);
        double Wq = lambdaEff > 0 ? Lq / lambdaEff : 0.0;
        double W = lambdaEff > 0 ? L / lambdaEff : 0.0;
        return new Result(L, Lq, Wq, W, util, lambdaEff, pK);
    }
}
