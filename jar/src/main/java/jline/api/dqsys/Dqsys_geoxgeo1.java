/**
 * @file Geo^X/Geo/1 discrete-time batch-arrival queueing system analysis
 *
 * Exact stationary mean analysis of the slotted single-server queue with
 * Bernoulli batch arrivals and geometric service, derived from the queue-length
 * probability generating function.
 */
package jline.api.dqsys;

/**
 * Geo^X/Geo/1 discrete-time queueing system analysis.
 *
 * <p>Generalizes {@link Dqsys_geogeo1} from single to batch arrivals. In every
 * slot a batch arrives with probability {@code a}; the batch size {@code X} is
 * supported on {@code {1,2,...}} and drawn independently. Service is geometric
 * with per-slot completion probability {@code s} and a single server, so the
 * arrival rate is {@code lambda = a E[X]} and stability requires
 * {@code lambda < s}.
 *
 * <p>With {@code A(z) = 1 - a + a X(z)} the pgf of the number of jobs arriving
 * in one slot, the slot-boundary queue length obeys
 * {@code X(t+1) = X(t) - D(t) + A(t)} with the departure resolved first
 * (Daduna's LA- and D/A-rules), giving
 *
 * <pre>
 *   P(z) = p0 s (z-1) A(z) / ( z - A(z)(s + (1-s)z) ),   p0 = 1 - lambda/s
 * </pre>
 *
 * <p>Differentiating at {@code z = 1} yields the mean, which needs only the
 * first two factorial moments of the batch:
 *
 * <pre>
 *   E[N] = lambda + ( a E[X(X-1)]/2 + lambda(1-s) ) / (s - lambda)
 * </pre>
 *
 * <p>At {@code X == 1} this collapses to {@code a(1-a)/(s-a)}, the Geo/Geo/1
 * result, so {@link Dqsys_geogeo1} is the degenerate-batch special case. The
 * batch enters only through its first two factorial moments, so the primary
 * entry point takes those directly and any batch law is supported; the
 * three-argument overload is the geometric-batch convenience form that gives
 * the model its name.
 *
 * <p>Reference: H. Daduna, Queueing Networks with Discrete Time Scale, LNCS
 * 2046, Springer 2001, chapter 6 (networks with batch movements) for the
 * discrete-time batch framework; the single-node pgf above is the standard
 * discrete-time M/G/1-type derivation.
 *
 * @since LINE 3.1.0
 */
public final class Dqsys_geoxgeo1 {
    private Dqsys_geoxgeo1() {}

    /**
     * Analyzes a Geo^X/Geo/1 queue with geometrically distributed batch sizes
     * under the LAS-DA convention.
     *
     * @param a    per-slot probability that a batch arrives, 0 &lt; a &lt;= 1
     * @param beta batch-size geometric parameter; the batch is supported on
     *             {1,2,...} with mean 1/beta, 0 &lt; beta &lt;= 1
     * @param s    per-slot service completion probability, 0 &lt; s &lt;= 1
     * @return the stationary metrics
     */
    public static GeoXGeo1Result dqsys_geoxgeo1(double a, double beta, double s) {
        return dqsys_geoxgeo1(a, beta, s, GeoGeo1Convention.LAS_DA);
    }

    /**
     * Analyzes a Geo^X/Geo/1 queue with geometrically distributed batch sizes.
     *
     * <p>For a geometric batch on {1,2,...}, {@code E[X] = 1/beta} and
     * {@code E[X(X-1)] = 2(1-beta)/beta^2}. At {@code beta = 1} the batch is
     * always a single job and the result equals {@link Dqsys_geogeo1}.
     *
     * @param a          per-slot probability that a batch arrives
     * @param beta       batch-size geometric parameter
     * @param s          per-slot service completion probability
     * @param convention observation epoch
     * @return the stationary metrics
     */
    public static GeoXGeo1Result dqsys_geoxgeo1(double a, double beta, double s,
                                               GeoGeo1Convention convention) {
        if (!(beta > 0.0) || beta > 1.0) {
            throw new IllegalArgumentException("Batch geometric parameter beta=" + beta
                    + " must lie in (0,1]");
        }
        double batchMean = 1.0 / beta;
        double batchSecondFactorial = 2.0 * (1.0 - beta) / (beta * beta);
        return dqsys_geoxgeo1_moments(a, batchMean, batchSecondFactorial, s, convention);
    }

    /**
     * Analyzes a Geo^X/Geo/1 queue for an arbitrary batch-size law, specified by
     * its first two factorial moments.
     *
     * @param a                    per-slot probability that a batch arrives,
     *                             0 &lt; a &lt;= 1
     * @param batchMean            E[X], must be at least 1 since a batch that
     *                             arrives carries at least one job
     * @param batchSecondFactorial E[X(X-1)], must be non-negative and at least
     *                             {@code batchMean^2 - batchMean} by
     *                             Cauchy-Schwarz
     * @param s                    per-slot service completion probability,
     *                             0 &lt; s &lt;= 1
     * @param convention           observation epoch
     * @return the stationary metrics
     */
    public static GeoXGeo1Result dqsys_geoxgeo1_moments(double a, double batchMean,
                                                       double batchSecondFactorial, double s,
                                                       GeoGeo1Convention convention) {
        if (convention == null) {
            throw new IllegalArgumentException("Convention must not be null");
        }
        if (!(a > 0.0) || a > 1.0) {
            throw new IllegalArgumentException("Batch arrival probability a=" + a
                    + " must lie in (0,1]");
        }
        if (!(s > 0.0) || s > 1.0) {
            throw new IllegalArgumentException("Service probability s=" + s
                    + " must lie in (0,1]");
        }
        if (!(batchMean >= 1.0)) {
            throw new IllegalArgumentException("Mean batch size E[X]=" + batchMean
                    + " must be at least 1: a batch that arrives carries at least one job");
        }
        if (batchSecondFactorial < 0.0) {
            throw new IllegalArgumentException("E[X(X-1)]=" + batchSecondFactorial
                    + " must be non-negative");
        }
        // E[X^2] = E[X(X-1)] + E[X] >= E[X]^2, so a smaller second factorial
        // moment describes no random variable at all.
        double minSecondFactorial = batchMean * batchMean - batchMean;
        if (batchSecondFactorial < minSecondFactorial - 1e-9 * Math.max(1.0, minSecondFactorial)) {
            throw new IllegalArgumentException("E[X(X-1)]=" + batchSecondFactorial
                    + " is below E[X]^2-E[X]=" + minSecondFactorial
                    + ", so the batch moments describe no random variable");
        }

        double lambda = a * batchMean;
        if (!(lambda < s)) {
            throw new IllegalArgumentException("Load lambda/s=" + (lambda / s)
                    + " must be strictly less than 1");
        }

        double rho = lambda / s;
        double boundaryEmptyProb = 1.0 - rho;

        // P'(1) from the generating function; the batch enters only through its
        // first two factorial moments.
        double meanAtBoundary = lambda
                + (a * batchSecondFactorial / 2.0 + lambda * (1.0 - s)) / (s - lambda);

        // The queueing delay does not depend on the observation epoch; only the
        // accounting of the slot in which service takes place does.
        double meanSojournAtBoundary = meanAtBoundary / lambda;
        double meanWaitingTime = meanSojournAtBoundary - 1.0 / s;
        double meanWaitingQueue = lambda * meanWaitingTime;

        double meanQueueLength;
        double meanSojournTime;
        double meanServiceTime;
        if (convention == GeoGeo1Convention.LAS_DA) {
            meanQueueLength = meanAtBoundary;
            meanSojournTime = meanSojournAtBoundary;
            meanServiceTime = 1.0 / s;
        } else {
            // One departure earlier: the epoch drops exactly the departures of
            // the slot, whose rate is lambda, hence one slot of sojourn.
            meanQueueLength = meanAtBoundary - lambda;
            meanSojournTime = meanSojournAtBoundary - 1.0;
            meanServiceTime = (1.0 - s) / s;
        }

        return new GeoXGeo1Result(convention, a, batchMean, batchSecondFactorial, s, lambda,
                rho, boundaryEmptyProb, meanQueueLength, meanWaitingQueue, meanSojournTime,
                meanWaitingTime, meanServiceTime);
    }
}
