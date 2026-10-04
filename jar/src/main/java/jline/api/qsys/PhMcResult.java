/**
 * @file Result of PH/M/c queue analysis (matrix-geometric)
 *
 * @since LINE 3.1.0
 */
package jline.api.qsys;

/**
 * Result of PH/M/c queue analysis.
 */
public final class PhMcResult {
    private final double meanQueueLength;
    private final double meanWaitingQueue;
    private final double meanWaitingTime;
    private final double meanSojournTime;
    private final double utilization;

    public PhMcResult(double meanQueueLength, double meanWaitingQueue, double meanWaitingTime,
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
        if (!(o instanceof PhMcResult)) return false;
        PhMcResult that = (PhMcResult) o;
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
        return "PhMcResult(meanQueueLength=" + meanQueueLength
                + ", meanWaitingQueue=" + meanWaitingQueue
                + ", meanWaitingTime=" + meanWaitingTime
                + ", meanSojournTime=" + meanSojournTime
                + ", utilization=" + utilization + ")";
    }
}
