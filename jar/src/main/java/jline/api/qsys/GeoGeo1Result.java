/**
 * @file Result of discrete-time Geo/Geo/1 queue analysis
 *
 * @since LINE 3.1.0
 */
package jline.api.qsys;

/**
 * Result of a discrete-time Geo/Geo/1 queue analysis.
 *
 * <p>All time quantities are expressed in slots. The queue-length distribution
 * is geometric in both conventions but with different forms, so it is exposed
 * through {@link #pmf(int)} rather than as a fixed pair of parameters.
 */
public final class GeoGeo1Result {

    private final GeoGeo1Convention convention;
    private final double arrivalProb;
    private final double serviceProb;
    private final double utilization;
    private final double throughput;
    private final double emptyProb;
    private final double ratio;
    private final double meanQueueLength;
    private final double meanWaitingQueue;
    private final double meanSojournTime;
    private final double meanWaitingTime;
    private final double meanServiceTime;

    public GeoGeo1Result(GeoGeo1Convention convention, double arrivalProb, double serviceProb,
                         double utilization, double throughput, double emptyProb, double ratio,
                         double meanQueueLength, double meanWaitingQueue, double meanSojournTime,
                         double meanWaitingTime, double meanServiceTime) {
        this.convention = convention;
        this.arrivalProb = arrivalProb;
        this.serviceProb = serviceProb;
        this.utilization = utilization;
        this.throughput = throughput;
        this.emptyProb = emptyProb;
        this.ratio = ratio;
        this.meanQueueLength = meanQueueLength;
        this.meanWaitingQueue = meanWaitingQueue;
        this.meanSojournTime = meanSojournTime;
        this.meanWaitingTime = meanWaitingTime;
        this.meanServiceTime = meanServiceTime;
    }

    /** Convention under which the metrics were computed. */
    public GeoGeo1Convention getConvention() { return convention; }

    /** Per-slot arrival probability a. */
    public double getArrivalProb() { return arrivalProb; }

    /** Per-slot service completion probability s. */
    public double getServiceProb() { return serviceProb; }

    /** Fraction of slots in which the server is serving, a/s in both conventions. */
    public double getUtilization() { return utilization; }

    /** Departures per slot, equal to a in steady state. */
    public double getThroughput() { return throughput; }

    /** Probability that the system is empty at a slot boundary. */
    public double getEmptyProb() { return emptyProb; }

    /** Geometric decay ratio of the queue-length tail, a(1-s)/(s(1-a)). */
    public double getRatio() { return ratio; }

    /** Mean number of jobs in the system. */
    public double getMeanQueueLength() { return meanQueueLength; }

    /** Mean number of jobs waiting, i.e. not in service. */
    public double getMeanWaitingQueue() { return meanWaitingQueue; }

    /** Mean sojourn time in slots. */
    public double getMeanSojournTime() { return meanSojournTime; }

    /** Mean waiting time in slots, identical in both conventions. */
    public double getMeanWaitingTime() { return meanWaitingTime; }

    /**
     * Mean service time in slots: 1/s under {@link GeoGeo1Convention#LAS_DA}
     * (support {1,2,...}) and (1-s)/s under {@link GeoGeo1Convention#EAS}
     * (support {0,1,...}).
     */
    public double getMeanServiceTime() { return meanServiceTime; }

    /**
     * Stationary probability of finding n jobs in the system at a slot boundary.
     *
     * @param n number of jobs, n &gt;= 0
     * @return the stationary probability
     */
    public double pmf(int n) {
        if (n < 0) throw new IllegalArgumentException("Queue length must be non-negative");
        if (convention == GeoGeo1Convention.EAS) {
            return (1.0 - ratio) * Math.pow(ratio, n);
        }
        if (n == 0) return emptyProb;
        return emptyProb * (utilization / (1.0 - arrivalProb)) * Math.pow(ratio, n - 1);
    }

    @Override
    public String toString() {
        return "GeoGeo1Result(convention=" + convention
                + ", a=" + arrivalProb + ", s=" + serviceProb
                + ", utilization=" + utilization
                + ", meanQueueLength=" + meanQueueLength
                + ", meanWaitingQueue=" + meanWaitingQueue
                + ", meanSojournTime=" + meanSojournTime
                + ", meanWaitingTime=" + meanWaitingTime + ")";
    }
}
