/**
 * @file Result of M/D/c queue analysis via Crommelin's embedded DTMC
 *
 * @since LINE 3.1.0
 */
package jline.api.qsys;

/**
 * Result of M/D/c queue analysis via Crommelin's embedded DTMC.
 */
public final class MDcCrommelinResult {
    private final double meanQueueLength;   // E[N] = number in system
    private final double meanWaitingQueue;  // Lq = E[(N-c)+]
    private final double meanWaitingTime;   // Wq = Lq / lambda
    private final double meanSojournTime;   // W = Wq + s
    private final double utilization;       // rho = lambda * s / c

    public MDcCrommelinResult(double meanQueueLength, double meanWaitingQueue, double meanWaitingTime,
                              double meanSojournTime, double utilization) {
        this.meanQueueLength = meanQueueLength;
        this.meanWaitingQueue = meanWaitingQueue;
        this.meanWaitingTime = meanWaitingTime;
        this.meanSojournTime = meanSojournTime;
        this.utilization = utilization;
    }

    public double getMeanQueueLength() { return meanQueueLength; }
    public double getMeanWaitingQueue() { return meanWaitingQueue; }
    public double getMeanWaitingTime() { return meanWaitingTime; }
    public double getMeanSojournTime() { return meanSojournTime; }
    public double getUtilization() { return utilization; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof MDcCrommelinResult)) return false;
        MDcCrommelinResult that = (MDcCrommelinResult) o;
        return Double.compare(that.meanQueueLength, meanQueueLength) == 0
                && Double.compare(that.meanWaitingQueue, meanWaitingQueue) == 0
                && Double.compare(that.meanWaitingTime, meanWaitingTime) == 0
                && Double.compare(that.meanSojournTime, meanSojournTime) == 0
                && Double.compare(that.utilization, utilization) == 0;
    }

    @Override
    public int hashCode() {
        int r = Double.hashCode(meanQueueLength);
        r = 31 * r + Double.hashCode(meanWaitingQueue);
        r = 31 * r + Double.hashCode(meanWaitingTime);
        r = 31 * r + Double.hashCode(meanSojournTime);
        r = 31 * r + Double.hashCode(utilization);
        return r;
    }

    @Override
    public String toString() {
        return "MDcCrommelinResult(meanQueueLength=" + meanQueueLength
                + ", meanWaitingQueue=" + meanWaitingQueue
                + ", meanWaitingTime=" + meanWaitingTime
                + ", meanSojournTime=" + meanSojournTime
                + ", utilization=" + utilization + ")";
    }
}
