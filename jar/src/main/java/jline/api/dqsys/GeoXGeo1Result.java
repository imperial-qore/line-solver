/**
 * @file Result of discrete-time Geo^X/Geo/1 batch-arrival queue analysis
 *
 * @since LINE 3.1.0
 */
package jline.api.dqsys;

/**
 * Result of a discrete-time Geo^X/Geo/1 queue analysis.
 *
 * <p>All time quantities are in slots. Unlike the single-arrival Geo/Geo/1, the
 * stationary queue length has no elementary closed form for a general batch
 * law, so the distribution is exposed through its probability generating
 * function {@link #pgf(double)} rather than a pmf. The mean measures are exact.
 */
public final class GeoXGeo1Result {

    private final GeoGeo1Convention convention;
    private final double batchArrivalProb;
    private final double batchMean;
    private final double batchSecondFactorialMoment;
    private final double serviceProb;
    private final double arrivalRate;
    private final double utilization;
    private final double emptyProb;
    private final double meanQueueLength;
    private final double meanWaitingQueue;
    private final double meanSojournTime;
    private final double meanWaitingTime;
    private final double meanServiceTime;

    public GeoXGeo1Result(GeoGeo1Convention convention, double batchArrivalProb, double batchMean,
                          double batchSecondFactorialMoment, double serviceProb, double arrivalRate,
                          double utilization, double emptyProb, double meanQueueLength,
                          double meanWaitingQueue, double meanSojournTime, double meanWaitingTime,
                          double meanServiceTime) {
        this.convention = convention;
        this.batchArrivalProb = batchArrivalProb;
        this.batchMean = batchMean;
        this.batchSecondFactorialMoment = batchSecondFactorialMoment;
        this.serviceProb = serviceProb;
        this.arrivalRate = arrivalRate;
        this.utilization = utilization;
        this.emptyProb = emptyProb;
        this.meanQueueLength = meanQueueLength;
        this.meanWaitingQueue = meanWaitingQueue;
        this.meanSojournTime = meanSojournTime;
        this.meanWaitingTime = meanWaitingTime;
        this.meanServiceTime = meanServiceTime;
    }

    /** Observation epoch the metrics refer to. */
    public GeoGeo1Convention getConvention() { return convention; }

    /** Per-slot probability that a batch arrives. */
    public double getBatchArrivalProb() { return batchArrivalProb; }

    /** Mean batch size E[X], conditional on a batch arriving. */
    public double getBatchMean() { return batchMean; }

    /** Second factorial moment E[X(X-1)] of the batch size. */
    public double getBatchSecondFactorialMoment() { return batchSecondFactorialMoment; }

    /** Per-slot service completion probability s. */
    public double getServiceProb() { return serviceProb; }

    /** Jobs arriving per slot, a*E[X]. */
    public double getArrivalRate() { return arrivalRate; }

    /** Departures per slot, equal to the arrival rate in steady state. */
    public double getThroughput() { return arrivalRate; }

    /** Fraction of slots in which the server is serving, lambda/s. */
    public double getUtilization() { return utilization; }

    /**
     * Probability that the system is empty AT THE SLOT BOUNDARY, {@code 1 - lambda/s}.
     *
     * <p>This is deliberately not epoch-dependent. It is the constant that
     * appears in the generating function under both conventions. The empty
     * probability at the EAS epoch is {@code p0 + s*pi_1}, which has no
     * elementary form for a general batch law; obtain it from {@link #pgf} in
     * the limit if it is needed.
     */
    public double getBoundaryEmptyProb() { return emptyProb; }

    /** Mean number of jobs in the system. */
    public double getMeanQueueLength() { return meanQueueLength; }

    /** Mean number of jobs waiting, i.e. not in service. */
    public double getMeanWaitingQueue() { return meanWaitingQueue; }

    /** Mean sojourn time in slots, per job (not per batch). */
    public double getMeanSojournTime() { return meanSojournTime; }

    /** Mean waiting time in slots, per job. */
    public double getMeanWaitingTime() { return meanWaitingTime; }

    /** Mean service time in slots, 1/s under LAS_DA and (1-s)/s under EAS. */
    public double getMeanServiceTime() { return meanServiceTime; }

    /**
     * Probability generating function of the stationary queue length at the
     * observation epoch, evaluated by the caller-supplied batch pgf.
     *
     * <p>Under LAS_DA, with {@code A(z)} the pgf of the number of jobs arriving
     * in one slot and {@code p0 = 1 - lambda/s},
     *
     * <pre>
     *   P(z) = p0 s (z-1) A(z) / ( z - A(z)(s + (1-s)z) )
     * </pre>
     *
     * Under EAS the epoch is one departure earlier, which divides out the
     * departure step: {@code P_EAS(z) = P(z) (s/z + 1 - s) + p0 s (1 - 1/z)}.
     *
     * @param z            evaluation point, {@code 0 < z <= 1}
     * @param batchPgfAtZ  the slot-arrival pgf {@code A(z)} at the same z
     * @return {@code P(z)}
     */
    public double pgf(double z, double batchPgfAtZ) {
        if (!(z > 0.0) || z > 1.0) {
            throw new IllegalArgumentException("PGF argument z=" + z + " must lie in (0,1]");
        }
        if (z == 1.0) {
            return 1.0;
        }
        double s = serviceProb;
        double denom = z - batchPgfAtZ * (s + (1.0 - s) * z);
        double boundary = emptyProbAtBoundary() * s * (z - 1.0) * batchPgfAtZ / denom;
        if (convention == GeoGeo1Convention.LAS_DA) {
            return boundary;
        }
        return boundary * (s / z + 1.0 - s) + emptyProbAtBoundary() * s * (1.0 - 1.0 / z);
    }

    /**
     * The empty probability at the slot boundary, {@code 1 - lambda/s}, which is
     * the constant appearing in the pgf regardless of the reported epoch.
     */
    private double emptyProbAtBoundary() {
        return 1.0 - utilization;
    }

    @Override
    public String toString() {
        return "GeoXGeo1Result(convention=" + convention
                + ", a=" + batchArrivalProb + ", E[X]=" + batchMean
                + ", s=" + serviceProb + ", lambda=" + arrivalRate
                + ", utilization=" + utilization
                + ", meanQueueLength=" + meanQueueLength
                + ", meanSojournTime=" + meanSojournTime + ")";
    }
}
