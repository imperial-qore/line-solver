/**
 * @file Result of PH/M/1 queue analysis
 *
 * @since LINE 3.1.0
 */
package jline.api.qsys;

/**
 * Result of PH/M/1 queue analysis.
 */
public final class PhM1Result {
    private final double meanQueueLength;
    private final double meanWaitingQueue;
    private final double meanWaitingTime;
    private final double meanSojournTime;
    private final double utilization;
    private final double sigma;

    public PhM1Result(double meanQueueLength, double meanWaitingQueue, double meanWaitingTime,
                      double meanSojournTime, double utilization, double sigma) {
        this.meanQueueLength = meanQueueLength;
        this.meanWaitingQueue = meanWaitingQueue;
        this.meanWaitingTime = meanWaitingTime;
        this.meanSojournTime = meanSojournTime;
        this.utilization = utilization;
        this.sigma = sigma;
    }

    public double getMeanQueueLength() { return meanQueueLength; }
    public double getMeanWaitingQueue() { return meanWaitingQueue; }
    public double getMeanWaitingTime() { return meanWaitingTime; }
    public double getMeanSojournTime() { return meanSojournTime; }
    public double getUtilization() { return utilization; }
    public double getSigma() { return sigma; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof PhM1Result)) return false;
        PhM1Result that = (PhM1Result) o;
        return Double.compare(that.meanQueueLength, meanQueueLength) == 0
                && Double.compare(that.meanWaitingQueue, meanWaitingQueue) == 0
                && Double.compare(that.meanWaitingTime, meanWaitingTime) == 0
                && Double.compare(that.meanSojournTime, meanSojournTime) == 0
                && Double.compare(that.utilization, utilization) == 0
                && Double.compare(that.sigma, sigma) == 0;
    }

    @Override
    public int hashCode() {
        int r = Double.hashCode(meanQueueLength);
        r = 31 * r + Double.hashCode(meanWaitingQueue);
        r = 31 * r + Double.hashCode(meanWaitingTime);
        r = 31 * r + Double.hashCode(meanSojournTime);
        r = 31 * r + Double.hashCode(utilization);
        r = 31 * r + Double.hashCode(sigma);
        return r;
    }

    @Override
    public String toString() {
        return "PhM1Result(meanQueueLength=" + meanQueueLength
                + ", meanWaitingQueue=" + meanWaitingQueue
                + ", meanWaitingTime=" + meanWaitingTime
                + ", meanSojournTime=" + meanSojournTime
                + ", utilization=" + utilization
                + ", sigma=" + sigma + ")";
    }
}
