/**
 * @file Geo/Geo/1 discrete-time queueing system analysis
 *
 * Exact stationary analysis of the single-server slotted queue with Bernoulli
 * arrivals and geometric service, in both the late-arrival-delayed-access and
 * the early-arrival conventions.
 */
package jline.api.qsys;

/**
 * Geo/Geo/1 discrete-time queueing system analysis.
 *
 * <p>Time advances in slots of unit length. In every slot an arrival occurs with
 * probability {@code a} and, if the server is engaged, a service completion
 * occurs with probability {@code s}, both independently of everything else. The
 * system content at slot boundaries is then a discrete birth-death chain whose
 * stationary distribution is geometric, so all mean measures are available in
 * closed form.
 *
 * <p>Under {@link GeoGeo1Convention#LAS_DA} the chain is
 *
 * <pre>
 *   0 -&gt; 1        with probability a
 *   n -&gt; n-1      with probability s(1-a)          (n &gt;= 1)
 *   n -&gt; n+1      with probability a(1-s)          (n &gt;= 1)
 * </pre>
 *
 * giving {@code pi_0 = 1-rho} and {@code pi_n = (1-rho) rho/(1-a) r^(n-1)} with
 * {@code rho = a/s} and {@code r = a(1-s)/(s(1-a))}. Under
 * {@link GeoGeo1Convention#EAS} an arrival may also depart within its own slot,
 * so the upward transition has probability {@code a(1-s)} from every state
 * including the empty one, giving the plain geometric {@code pi_n = (1-r) r^n}.
 *
 * <p>The two conventions describe the same queueing delay and differ only in
 * whether the arrival slot is counted as part of the sojourn: the mean sojourn
 * times differ by exactly one slot and the mean waiting times coincide. See
 * {@link GeoGeo1Convention} for the exact observation epochs.
 *
 * <p>The LAS-DA branch reproduces H. Daduna, Queueing Networks with Discrete
 * Time Scale, LNCS 2046, Springer 2001, corollary 2.7 term for term: his
 * {@code pi(n) = (1 - b/p) (bq/(cp))^n (1/q)^[n>0]}, with {@code b = a},
 * {@code p = s}, {@code c = 1-a} and {@code q = 1-s}, is the form below, and
 * his tail ratio {@code bq/(cp)} is {@code r}. He likewise states the service
 * time is geometric on {1,2,...}, which is the support jline's
 * {@code Geometric} uses.
 *
 * @since LINE 3.1.0
 */
public final class Qsys_geogeo1 {
    private Qsys_geogeo1() {}

    /**
     * Analyzes a Geo/Geo/1 queue under the late-arrival-delayed-access
     * convention, which is the one realized by a slotted LDES simulation of a
     * {@code Geometric(a)} source feeding a {@code Geometric(s)} FCFS queue.
     *
     * @param a per-slot arrival probability, 0 &lt; a &lt; s
     * @param s per-slot service completion probability, 0 &lt; s &lt;= 1
     * @return the stationary metrics
     */
    public static GeoGeo1Result qsys_geogeo1(double a, double s) {
        return qsys_geogeo1(a, s, GeoGeo1Convention.LAS_DA);
    }

    /**
     * Analyzes a Geo/Geo/1 queue under the requested slot-boundary convention.
     *
     * @param a          per-slot arrival probability, 0 &lt; a &lt; s
     * @param s          per-slot service completion probability, 0 &lt; s &lt;= 1
     * @param convention slot-boundary convention
     * @return the stationary metrics
     */
    public static GeoGeo1Result qsys_geogeo1(double a, double s, GeoGeo1Convention convention) {
        if (convention == null) {
            throw new IllegalArgumentException("Convention must not be null");
        }
        if (!(a > 0.0) || a > 1.0) {
            throw new IllegalArgumentException("Arrival probability a=" + a + " must lie in (0,1]");
        }
        if (!(s > 0.0) || s > 1.0) {
            throw new IllegalArgumentException("Service probability s=" + s + " must lie in (0,1]");
        }
        if (!(a < s)) {
            throw new IllegalArgumentException("Load a/s=" + (a / s) + " must be strictly less than 1");
        }

        double rho = a / s;
        double ratio = a * (1.0 - s) / (s * (1.0 - a));

        // The queueing delay does not depend on the convention: only the
        // accounting of the slot in which service takes place does.
        double meanWaitingTime = a * (1.0 - s) / (s * (s - a));
        double meanWaitingQueue = a * meanWaitingTime;

        double emptyProb;
        double meanQueueLength;
        double meanSojournTime;
        double meanServiceTime;
        if (convention == GeoGeo1Convention.LAS_DA) {
            emptyProb = 1.0 - rho;
            meanQueueLength = a * (1.0 - a) / (s - a);
            meanSojournTime = (1.0 - a) / (s - a);
            meanServiceTime = 1.0 / s;
        } else {
            emptyProb = 1.0 - ratio;
            meanQueueLength = a * (1.0 - s) / (s - a);
            meanSojournTime = (1.0 - s) / (s - a);
            meanServiceTime = (1.0 - s) / s;
        }

        return new GeoGeo1Result(convention, a, s, rho, a, emptyProb, ratio,
                meanQueueLength, meanWaitingQueue, meanSojournTime, meanWaitingTime,
                meanServiceTime);
    }
}
