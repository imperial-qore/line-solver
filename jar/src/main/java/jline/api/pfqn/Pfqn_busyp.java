/**
 * Mean Busy Period of Order n for a Subnetwork of a Product-Form Network
 *
 * @since LINE 3.0
 */
package jline.api.pfqn;

import jline.util.matrix.Matrix;

import java.util.ArrayList;
import java.util.List;

/**
 * Mean busy period of order n for a subnetwork, after H. Daduna, "Busy Periods
 * for Subnetworks in Stochastic Networks: Mean Value Analysis", J. ACM 35(3),
 * 1988: Theorem 1 for a closed Gordon-Newell network and Theorem 3 for an open
 * Jackson network. Both are evaluated in the log domain, which serves the same
 * purpose as the ratio recursions of Corollaries 2 and 4, namely avoiding the
 * overflow of the individual normalizing constants.
 *
 * The paper is single-chain: alpha is the stochastic solution of x*P = x for a
 * closed network and the solution of x = gamma + x*P for an open one, and every
 * node is a state-dependent single-server FCFS station. By the insensitivity
 * noted in Section 5 of the paper the result depends on the service processes
 * only through the rates mu.
 */
public final class Pfqn_busyp {
    private Pfqn_busyp() {}

    /** Default relative tolerance of the open-network tail truncation. */
    public static final double DEFAULT_TOL = 1e-12;

    /**
     * The load-dependent rate of one node, as a FUNCTION of the jobs it holds.
     *
     * A DENSE TABLE CANNOT EXPRESS AN INFINITE SERVER. The open-network form
     * below truncates its tail adaptively and may need the rate at an order far
     * beyond any table a caller would build, and the Matrix overloads answer
     * such a query with their last column -- correct for the multiserver
     * staircase, which saturates at c/S, and wrong for a delay station, whose
     * rate is k/S and grows without bound. The MATLAB and C++ ports of this
     * routine have always taken a function for that reason; this is its Java
     * spelling, and the Matrix overloads are now a wrapper over it.
     */
    public interface RateFunction {
        /**
         * @param node zero-based node index in the ORIGINAL numbering handed to
         *             pfqn_busyp, not the subnetwork's
         * @param jobs number of jobs held, 1 or more
         * @return the total service rate of the node at that occupancy
         */
        double rate(int node, int jobs);
    }

    /**
     * Mean busy period of order n for the subnetwork, that is the time from the
     * instant a job entering the subnetwork finds n-1 jobs in it up to the next
     * instant when fewer than n jobs remain in it.
     *
     * @param alpha  relative arrival rates (1xJ)
     * @param mu     load-dependent service rates (JxK), mu(j,k-1) with k jobs at node j;
     *               a table shorter than the population keeps its last rate
     * @param P      routing matrix (JxJ)
     * @param N      population, Double.POSITIVE_INFINITY for an open network
     * @param subnet zero-based indexes of the nodes forming the subnetwork
     * @param n      busy period orders, 1 &lt;= n &lt;= N
     * @param gamma  external arrival rates (1xJ), null for a closed network
     * @param tol    relative tolerance of the open-network tail truncation
     * @return mean busy period duration for each requested order
     */
    public static double[] pfqn_busyp(Matrix alpha, Matrix mu, Matrix P, double N,
                                      int[] subnet, int[] n, Matrix gamma, double tol) {
        return pfqn_busyp(alpha, asRateFunction(mu), P, N, subnet, n, gamma, tol);
    }

    /**
     * Rate-function form, the one the MATLAB and C++ ports take. Identical to
     * the Matrix form above except that mu is queried at the occupancy the
     * algorithm actually reaches, so a delay station in the subnetwork is exact
     * rather than clamped to the last column of a table.
     *
     * @param alpha  relative arrival rates (1xJ)
     * @param mu     load-dependent service rate of each node at each occupancy
     * @param P      routing matrix (JxJ)
     * @param N      population, Double.POSITIVE_INFINITY for an open network
     * @param subnet zero-based indexes of the nodes forming the subnetwork
     * @param n      busy period orders, 1 &lt;= n &lt;= N
     * @param gamma  external arrival rates (1xJ), null for a closed network
     * @param tol    relative tolerance of the open-network tail truncation
     * @return mean busy period duration for each requested order
     */
    public static double[] pfqn_busyp(Matrix alpha, RateFunction mu, Matrix P, double N,
                                      int[] subnet, int[] n, Matrix gamma, double tol) {
        int J = alpha.length();
        boolean isClosed = !Double.isInfinite(N);

        int[] target = sortedUnique(subnet);
        if (target.length == 0) {
            throw new IllegalArgumentException("The subnetwork must be non-empty.");
        }
        if (isClosed && target.length >= J) {
            // a closed network needs jobs outside the subnetwork to start a busy period
            throw new IllegalArgumentException(
                    "In a closed network the subnetwork must be a proper subset of the nodes.");
        }
        if (target[0] < 0 || target[target.length - 1] >= J) {
            throw new IllegalArgumentException("The subnetwork indexes are out of range.");
        }
        int[] compl = complement(J, target);

        int nmax = 0;
        for (int t = 0; t < n.length; t++) {
            if (n[t] < 1) {
                throw new IllegalArgumentException("The busy period order must be a positive integer.");
            }
            if (isClosed && n[t] > N) {
                throw new IllegalArgumentException("The busy period order must be an integer in 1..N.");
            }
            nmax = Math.max(nmax, n[t]);
        }
        if (!isClosed && gamma == null) {
            throw new IllegalArgumentException(
                    "An open network requires the external arrival rates gamma.");
        }

        // A(I) for a closed network, C(I) for an open one: both are the total rate at
        // which jobs enter the subnetwork from outside it, which is what starts a busy
        // period. The closed network has no external stream.
        double inflow = 0.0;
        for (int i = 0; i < compl.length; i++) {
            for (int j = 0; j < target.length; j++) {
                inflow += alpha.get(compl[i]) * P.get(compl[i], target[j]);
            }
        }
        if (gamma != null) {
            for (int j = 0; j < target.length; j++) {
                inflow += gamma.get(target[j]);
            }
        }
        if (inflow <= 0) {
            throw new IllegalArgumentException(
                    "No job ever enters the subnetwork, its busy period is undefined.");
        }

        double[] b = new double[n.length];
        if (isClosed) {
            int pop = (int) Math.round(N);
            double[] lG = lgvec(select(alpha, target), rates(mu, target, pop), pop);
            double[] lH = lgvec(select(alpha, compl), rates(mu, compl, pop), pop);
            for (int t = 0; t < n.length; t++) {
                // Theorem 1: sum_{m=n}^{N} G(m,I) H(N-m,I) over G(n-1,I) H(N-n,I) A(I)
                double[] terms = new double[pop - n[t] + 1];
                for (int m = n[t]; m <= pop; m++) {
                    terms[m - n[t]] = lG[m] + lH[pop - m];
                }
                b[t] = Math.exp(lse(terms) - lG[n[t] - 1] - lH[pop - n[t]] - Math.log(inflow));
            }
        } else {
            int K = truncate(select(alpha, target), mu, target, nmax, tol);
            double[] lG = lgvec(select(alpha, target), rates(mu, target, K), K);
            for (int t = 0; t < n.length; t++) {
                // Theorem 3: sum_{m=n}^{Inf} G(m,I) over G(n-1,I) C(I), the tail summed
                // up to the truncation order
                double[] terms = new double[K - n[t] + 1];
                for (int m = n[t]; m <= K; m++) {
                    terms[m - n[t]] = lG[m];
                }
                b[t] = Math.exp(lse(terms) - lG[n[t] - 1] - Math.log(inflow));
            }
        }
        return b;
    }

    /** Closed-network form with the default tolerance. */
    public static double[] pfqn_busyp(Matrix alpha, Matrix mu, Matrix P, double N,
                                      int[] subnet, int[] n) {
        return pfqn_busyp(alpha, mu, P, N, subnet, n, null, DEFAULT_TOL);
    }

    /** Single-order form. */
    public static double pfqn_busyp(Matrix alpha, Matrix mu, Matrix P, double N,
                                    int[] subnet, int n, Matrix gamma) {
        return pfqn_busyp(alpha, mu, P, N, subnet, new int[]{n}, gamma, DEFAULT_TOL)[0];
    }

    /** Rate-function form with the default tolerance. */
    public static double[] pfqn_busyp(Matrix alpha, RateFunction mu, Matrix P, double N,
                                      int[] subnet, int[] n, Matrix gamma) {
        return pfqn_busyp(alpha, mu, P, N, subnet, n, gamma, DEFAULT_TOL);
    }

    /**
     * Log normalizing constants of orders 0..K of a set of nodes: lg[m] is the log
     * of the sum over the compositions n_1+...+n_L = m of the product over the nodes
     * of prod_{k=1}^{n_i} alpha_i/mu_i(k), the G(m,I) and H(m,I) of the paper. The
     * nodes are convolved one at a time in the log domain.
     */
    private static double[] lgvec(double[] alpha, double[][] mu, int K) {
        double[] lg = new double[K + 1];
        for (int m = 1; m <= K; m++) {
            lg[m] = Double.NEGATIVE_INFINITY;
        }
        for (int i = 0; i < alpha.length; i++) {
            double[] li = new double[K + 1];
            double acc = 0.0;
            for (int m = 1; m <= K; m++) {
                acc += Math.log(alpha[i]) - Math.log(mu[i][m - 1]);
                li[m] = acc;
            }
            double[] lgnew = new double[K + 1];
            for (int m = 0; m <= K; m++) {
                double[] terms = new double[m + 1];
                for (int k = 0; k <= m; k++) {
                    terms[k] = lg[m - k] + li[k];
                }
                lgnew[m] = lse(terms);
            }
            lg = lgnew;
        }
        return lg;
    }

    /**
     * Truncation order of the open-network sum sum_{m&gt;=n} G(m,I). The subnetwork
     * terms decay geometrically once the rates saturate, so the truncation grows
     * until the geometric tail estimate is negligible against the partial sum.
     */
    private static int truncate(double[] alpha, RateFunction mu, int[] subnet, int nmax, double tol) {
        int K = Math.max(nmax + 8, 16);
        double[][] rows = rates(mu, subnet, K);
        double rhoMax = 0.0;
        for (int i = 0; i < alpha.length; i++) {
            rhoMax = Math.max(rhoMax, alpha[i] / rows[i][K - 1]);
        }
        if (rhoMax >= 1) {
            throw new IllegalArgumentException(
                    "The subnetwork is not stable, its busy period is infinite.");
        }
        while (true) {
            double[] lg = lgvec(alpha, rates(mu, subnet, K), K);
            // decay rate read off the last two orders, the exact ratio for a saturated
            // single-server subnetwork and an upper estimate otherwise
            double r = Math.exp(lg[K] - lg[K - 1]);
            if (!(r < 1)) {
                r = rhoMax;
            }
            double ltail = lg[K] + Math.log(r) - Math.log1p(-r);
            double[] partial = new double[K - nmax + 1];
            for (int m = nmax; m <= K; m++) {
                partial[m - nmax] = lg[m];
            }
            if (ltail - lse(partial) < Math.log(tol)) {
                return K;
            }
            K = 2 * K;
            if (K > 1e6) {
                throw new IllegalArgumentException(
                        "The open busy period sum did not converge, the subnetwork is nearly saturated.");
            }
        }
    }

    /**
     * Rates of the selected nodes for populations 1..K. A rate table shorter than K
     * keeps its last rate, the saturated-server convention.
     */
    private static double[][] rates(RateFunction mu, int[] idx, int K) {
        double[][] out = new double[idx.length][K];
        for (int i = 0; i < idx.length; i++) {
            for (int k = 0; k < K; k++) {
                out[i][k] = mu.rate(idx[i], k + 1);
            }
        }
        return out;
    }

    /**
     * The dense-table reading of mu: column k-1 holds the rate at k jobs and a
     * table shorter than the query keeps its last column, which is what every
     * Matrix-taking overload above has always meant.
     */
    private static RateFunction asRateFunction(final Matrix mu) {
        final int cols = mu.getNumCols();
        return new RateFunction() {
            @Override
            public double rate(int node, int jobs) {
                return mu.get(node, Math.min(jobs - 1, cols - 1));
            }
        };
    }

    /** log-sum-exp, stable when every entry is -infinity. */
    private static double lse(double[] v) {
        double m = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < v.length; i++) {
            if (v[i] > m) {
                m = v[i];
            }
        }
        if (Double.isInfinite(m) || Double.isNaN(m)) {
            return m;
        }
        double s = 0.0;
        for (int i = 0; i < v.length; i++) {
            s += Math.exp(v[i] - m);
        }
        return m + Math.log(s);
    }

    private static double[] select(Matrix alpha, int[] idx) {
        double[] out = new double[idx.length];
        for (int i = 0; i < idx.length; i++) {
            out[i] = alpha.get(idx[i]);
        }
        return out;
    }

    private static int[] sortedUnique(int[] v) {
        int[] copy = v.clone();
        java.util.Arrays.sort(copy);
        List<Integer> out = new ArrayList<Integer>();
        for (int i = 0; i < copy.length; i++) {
            if (i == 0 || copy[i] != copy[i - 1]) {
                out.add(copy[i]);
            }
        }
        int[] res = new int[out.size()];
        for (int i = 0; i < res.length; i++) {
            res[i] = out.get(i);
        }
        return res;
    }

    private static int[] complement(int J, int[] subnet) {
        boolean[] inSubnet = new boolean[J];
        for (int i = 0; i < subnet.length; i++) {
            inSubnet[subnet[i]] = true;
        }
        int count = 0;
        for (int j = 0; j < J; j++) {
            if (!inSubnet[j]) {
                count++;
            }
        }
        int[] out = new int[count];
        int t = 0;
        for (int j = 0; j < J; j++) {
            if (!inSubnet[j]) {
                out[t++] = j;
            }
        }
        return out;
    }
}
