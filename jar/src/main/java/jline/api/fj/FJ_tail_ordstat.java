package jline.api.fj;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.solvers.BrentSolver;

/**
 * Tail latency of a k-of-n (quorum) fork-join request.
 *
 * <p>The request forks into N parallel tasks and joins on the KREQ-th of them. KREQ = N is the
 * ordinary AND-join and reproduces {@link FJ_tail_forktail} exactly; KREQ = 1 is the first
 * completion.</p>
 *
 * <p>Each branch is the same black box ForkTail uses: its task response time is fitted by a
 * generalized exponential law F_i(x) = (1-exp(-x/beta_i))^alpha_i matched on the branch mean and
 * variance by {@link FJ_tail_forktail#geFit}. The request completes once KREQ of the N branches
 * have, so its law is the KREQ-th ORDER STATISTIC of independent, not identically distributed
 * branch times,</p>
 *
 * <pre>    F_X(x) = P(at least KREQ of the N branches are done by x),</pre>
 *
 * <p>the upper tail of a Poisson-binomial with success probabilities F_i(x). It is evaluated by
 * the convolution recurrence, which adds no cancellation, and inverted by bisection. At KREQ = N
 * it reduces term by term to prod_i F_i(x), the product of the branch CDFs that ForkTail
 * inverts, so a full join evaluates exactly as it did before this class existed.</p>
 *
 * <p>BRANCH INDEPENDENCE is assumed, as in ForkTail: the branches of one request are positively
 * correlated through their shared arrival instant, so the true quorum percentile is somewhat
 * larger than this one. The same heavy-traffic caveat applies.</p>
 *
 * <p>Port of matlab/src/api/fj/fj_tail_ordstat.m.</p>
 *
 * <p>References: M. Nguyen, S. Alesawi, N. Li, H. Che, H. Jiang, "ForkTail: A Black-Box Fork-Join
 * Tail Latency Prediction Model for User-Facing Datacenter Workloads", ACM HPDC 2018, for the
 * branch law; A. Thomasian, "Analysis of Fork/Join and Related Queueing Systems", ACM Computing
 * Surveys 47(2), Article 17, 2014, Sec. 3, for the quorum.</p>
 */
public class FJ_tail_ordstat {

    private FJ_tail_ordstat() {}

    /**
     * Tail latency of a quorum request over K identical branches.
     *
     * @param ET mean task response time of a branch
     * @param VT variance of the task response time of a branch
     * @param K number of branches
     * @param p percentile, a fraction in (0,1) or a percentage in (0,100)
     * @param kreq the join fires on the kreq-th branch
     * @return the predicted p-th percentile of the request response time
     */
    public static double fj_tail_ordstat(double ET, double VT, int K, double p, int kreq) {
        double[] et = new double[K];
        double[] vt = new double[K];
        for (int i = 0; i < K; i++) {
            et[i] = ET;
            vt[i] = VT;
        }
        // The AND-join bound comes from the HOMOGENEOUS closed form, not from
        // root-finding the identical-branch product: the two agree only to about
        // 1e-12, and at kreq = K this function must return the ForkTail value bit
        // for bit, as the other three codebases do.
        return core(et, vt, p, kreq,
                FJ_tail_forktail.fj_tail_forktail(ET, VT, K, normalizePercentile(p)));
    }

    /**
     * Tail latency of a quorum request over heterogeneous branches, one entry per branch.
     *
     * @param ET mean task response time per branch
     * @param VT variance of the task response time per branch
     * @param p percentile, a fraction in (0,1) or a percentage in (0,100)
     * @param kreq the join fires on the kreq-th branch
     * @return the predicted p-th percentile of the request response time
     */
    public static double fj_tail_ordstat(double[] ET, double[] VT, double p, final int kreq) {
        if (ET.length != VT.length) {
            throw new IllegalArgumentException("ET and VT must have the same number of entries.");
        }
        final double pp0 = normalizePercentile(p);
        final double andjoin = (ET.length == 1)
                ? FJ_tail_forktail.fj_tail_forktail(ET[0], VT[0], 1, pp0)
                : FJ_tail_forktail.fj_tail_forktail(ET, VT, pp0);
        return core(ET, VT, p, kreq, andjoin);
    }

    /**
     * The order-statistic inversion, given the AND-join percentile that bounds it above.
     *
     * @param ET mean task response time per branch
     * @param VT variance of the task response time per branch
     * @param p percentile, a fraction in (0,1) or a percentage in (0,100)
     * @param kreq the join fires on the kreq-th branch
     * @param andjoin the AND-join percentile of the same branches, the upper bracket
     * @return the predicted p-th percentile of the request response time
     */
    private static double core(double[] ET, double[] VT, double p, final int kreq,
                               double andjoin) {
        final int n = ET.length;
        if (kreq < 1 || kreq > n) {
            throw new IllegalArgumentException(
                    "The quorum must satisfy 1 <= kreq <= " + n + ". Got " + kreq + ".");
        }
        for (int i = 0; i < n; i++) {
            if (ET[i] <= 0 || VT[i] <= 0) {
                throw new IllegalArgumentException(
                        "The task response time mean and variance must be positive.");
            }
        }
        final double pp = normalizePercentile(p);

        final double[] alpha = new double[n];
        final double[] beta = new double[n];
        double xloBranch = Double.POSITIVE_INFINITY;
        for (int i = 0; i < n; i++) {
            FJ_tail_forktail.GEFit fit = FJ_tail_forktail.geFit(ET[i], VT[i]);
            alpha[i] = fit.alpha;
            beta[i] = fit.beta;
            // A SINGLE branch percentile is a lower bound for every quorum.
            xloBranch = Math.min(xloBranch,
                    -beta[i] * Math.log(1.0 - Math.pow(pp, 1.0 / alpha[i])));
        }

        // The AND-join percentile is an upper bound for every quorum, so the two
        // bracket the root without a search.
        final double xhi = andjoin;
        if (kreq == n) {
            return xhi;
        }

        UnivariateFunction residual = new UnivariateFunction() {
            @Override
            public double value(double x) {
                double[] u = new double[n];
                for (int i = 0; i < n; i++) {
                    u[i] = geCdf(x, alpha[i], beta[i]);
                }
                return poissbinUpper(u, kreq) - pp;
            }
        };

        double lo = xloBranch;
        // A single branch may already meet the percentile; the quorum is met earlier.
        while (residual.value(lo) > 0.0 && lo > Double.MIN_NORMAL) {
            lo /= 2.0;
        }
        return new BrentSolver(1e-12, 1e-12).solve(200, residual, lo, xhi);
    }

    /**
     * The generalized-exponential CDF (1-exp(-x/beta))^alpha, clamped into [0,1] so that a
     * rounding excursion cannot leave the binomial recurrence out of domain.
     */
    private static double geCdf(double x, double alpha, double beta) {
        double u = Math.exp(alpha * Math.log1p(-Math.exp(-x / beta)));
        if (Double.isNaN(u) || Double.isInfinite(u)) {
            return 0.0;
        }
        return Math.min(1.0, Math.max(0.0, u));
    }

    /**
     * P(at least kreq successes) for independent Bernoulli trials of success probabilities u, by
     * the convolution recurrence over the trials. Every term is non-negative, so no cancellation
     * is introduced.
     */
    private static double poissbinUpper(double[] u, int kreq) {
        final int n = u.length;
        double[] pmf = new double[n + 1];
        pmf[0] = 1.0;
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j >= 1; j--) {
                pmf[j] = pmf[j] * (1.0 - u[i]) + pmf[j - 1] * u[i];
            }
            pmf[0] *= (1.0 - u[i]);
        }
        double acc = 0.0;
        for (int j = kreq; j <= n; j++) {
            acc += pmf[j];
        }
        return acc;
    }

    private static double normalizePercentile(double p) {
        double pp = p;
        if (pp > 1.0) {
            pp = pp / 100.0;
        }
        if (pp <= 0.0 || pp >= 1.0) {
            throw new IllegalArgumentException(
                    "The percentile must lie strictly between 0 and 1 (or 0 and 100).");
        }
        return pp;
    }
}
