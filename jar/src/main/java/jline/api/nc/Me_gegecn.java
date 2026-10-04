/**
 * Censored GE/GE/c/K;N queue by entropy maximisation.
 *
 * Maximum Entropy solution of the single-class censored FCFS GE/GE/c/K;N
 * queue of Kouvatsos (1994) "Entropy Maximisation and Queueing Network
 * Models", Section 4.1, equations (4.1)-(4.3). The queue holds at most N
 * jobs and never fewer than K; arrivals finding N jobs present are turned
 * away and departures are not allowed from state K. For a queue embedded
 * in an open network K is always 0; a positive K arises in closed
 * networks, where it records the minimum occupancy forced by the
 * remaining stations being full.
 *
 * The ME state probabilities subject to normalisation, the marginal
 * utilisations, the mean queue length excluding J jobs and the full-buffer
 * probability coincide with the global balance solution
 * p(n) = p(K) G_n x^h(n) y^f(n). The Lagrangian coefficients are invariant
 * to N and K, so letting K -&gt; 0 and N -&gt; infinity recovers the stable
 * GE/GE/c solution used by {@link Me_oqn}.
 *
 * @since LINE 3.0
 */
package jline.api.nc;

import jline.io.InputOutput;

public final class Me_gegecn {
    private Me_gegecn() {}

    /**
     * Solves a censored GE/GE/c/K;N queue by entropy maximisation.
     *
     * @param lambda arrival rate offered to the queue, the arrivals that
     *               are turned away included
     * @param Ca     squared coefficient of variation of the interarrival
     *               times (Ca &gt;= 1, the GE distribution being undefined
     *               below 1)
     * @param mu     service rate of one server
     * @param Cs     squared coefficient of variation of the service times
     *               (Cs &gt;= 1)
     * @param c      number of servers, finite and at least one
     * @param K      minimum number of jobs in the queue
     * @param N      buffer capacity in jobs, in service included (N &gt; K)
     * @return the queue length distribution and its first moments
     */
    public static MeGegecnResult me_gegecn(double lambda, double Ca, double mu, double Cs,
                                           int c, int K, int N) {
        if (c < 1) {
            InputOutput.line_error("me_gegecn", "me_gegecn requires a finite number of servers c >= 1.");
        }
        if (N <= K) {
            InputOutput.line_error("me_gegecn", "me_gegecn requires N > K.");
        }
        if (Ca < 1 - 1e-12 || Cs < 1 - 1e-12) {
            InputOutput.line_error("me_gegecn",
                    "me_gegecn requires Ca >= 1 and Cs >= 1: the GE distribution is not defined for scv < 1.");
        }
        if (mu <= 0) {
            InputOutput.line_error("me_gegecn", "me_gegecn requires a positive service rate.");
        }

        double tau = 2.0 / (Ca + 1.0);
        double sigma = 2.0 / (Cs + 1.0);
        double rho = lambda / (c * mu);

        int J = Math.max(c, K + 1);
        double den1 = sigma * (1 - tau) + tau;
        double den2 = tau * rho * (1 - sigma) + sigma;

        // Lagrangian coefficients g(l), l = K+1,...,J, stored at index l-1
        double[] g = new double[J];
        for (int l = 0; l < J; l++) {
            g[l] = 1.0;
        }
        if (K < c - 1) {
            g[K] = tau * c * rho / ((K + 1) * den1);
        } else if (K == c - 1) {
            g[K] = tau * sigma * rho / den2;
        } else {
            g[K] = (den1 / den2) * tau * rho;
        }
        for (int l = K + 2; l <= J; l++) {
            if (l < J) {
                g[l - 1] = (tau * c * rho + (l - 1) * sigma * (1 - tau)) / (l * den1);
            } else {
                g[l - 1] = sigma * (tau * c * rho + (J - 1) * sigma * (1 - tau)) / (J * den2);
            }
        }

        double x = (tau * rho + sigma * (1 - tau)) / den2;
        double y = 1.0 / (1.0 - (1 - sigma) * x);

        // see _kb/03-api-layer.md for rationale
        int len = N - K + 1;
        double[] cumlogg = new double[J - K];
        double acc = 0.0;
        for (int l = K + 1; l <= J; l++) {
            acc += Math.log(g[l - 1]);
            cumlogg[l - K - 1] = acc;
        }
        double[] logp = new double[len];
        double maxlog = Double.NEGATIVE_INFINITY;
        for (int idx = 0; idx < len; idx++) {
            int n = K + idx;
            double v = 0.0;
            if (n > K) {
                int m = Math.max(K + 1, Math.min(c, n));
                v = cumlogg[m - K - 1];
            }
            v += Math.max(0, n - J) * Math.log(x);
            v += Math.max(0, n - N + 1) * Math.log(y);
            logp[idx] = v;
            if (v > maxlog) {
                maxlog = v;
            }
        }
        double[] p = new double[len];
        double sum = 0.0;
        for (int idx = 0; idx < len; idx++) {
            p[idx] = Math.exp(logp[idx] - maxlog);
            sum += p[idx];
        }
        double L = 0.0;
        double Ebusy = 0.0;
        for (int idx = 0; idx < len; idx++) {
            p[idx] /= sum;
            int n = K + idx;
            L += n * p[idx];
            Ebusy += Math.min(n, c) * p[idx];
        }
        double PB = me_gegecn_pb(p, K, N, c, Cs, Ca);
        return new MeGegecnResult(p, K, N, L, Ebusy / c, PB, L - Ebusy);
    }

    /**
     * Blocking probability seen by one arrival stream of a censored
     * GE/GE/c/K;N queue, equation (4.3) of the source.
     *
     * Because a GE arrival process is a batch process, an arrival can be
     * blocked while the queue holds fewer than N jobs: the factor
     * (1-tau)^(N-n) is the probability that the batch overflows the
     * residual room, and the extra factor of the first sum accounts for
     * the servers still idle. With a Poisson stream (Ca = 1) every term
     * but n = N vanishes and PB reduces to the PASTA value p(N).
     *
     * The stream scv is a per-stream quantity, so the same node solution
     * yields a different blocking probability for the external arrivals,
     * for the flow from each upstream station and for the flow released by
     * each holding node, which is how PBe, PB^i_j and PB^h_j are obtained
     * in the transfer-blocking algorithm of Tahilramani, Manjunath and
     * Bose (1999).
     *
     * @param p  ME queue length distribution, p[idx] = Pr{n = K+idx}
     * @param K  minimum number of jobs in the queue
     * @param N  buffer capacity in jobs
     * @param c  number of servers
     * @param Cs squared coefficient of variation of the service times
     * @param Ca squared coefficient of variation of the interarrival times
     *           of the stream whose blocking probability is requested
     * @return probability that an arrival of this stream finds the queue full
     */
    public static double me_gegecn_pb(double[] p, int K, int N, int c, double Cs, double Ca) {
        double tau = 2.0 / (Ca + 1.0);
        double sigma = 2.0 / (Cs + 1.0);
        int len = N - K + 1;
        double PB = 0.0;
        if (K < c) {
            int last = Math.min(c - K, len);
            for (int idx = 0; idx < last; idx++) {
                int n = K + idx;
                double fac = Math.pow(sigma / (sigma * (1 - tau) + tau), c - n);
                PB += Math.pow(1 - tau, N - n) * fac * p[idx];
            }
        }
        int lo = Math.max(c, K);
        if (lo <= N) {
            for (int idx = lo - K; idx < len; idx++) {
                int n = K + idx;
                PB += Math.pow(1 - tau, N - n) * p[idx];
            }
        }
        return PB;
    }
}
