/**
 * @file Result of MAP/D/c queue analysis
 *
 * @since LINE 3.1.0
 */
package jline.api.qsys;

import java.util.Objects;

import jline.util.matrix.Matrix;

/**
 * Result of MAP/D/c queue analysis.
 */
public final class QsysMapDcResult {
    private final double meanQueueLength;
    private final double meanWaitingTime;
    private final double meanSojournTime;
    private final double utilization;
    private final Matrix queueLengthDist;
    private final Matrix waitingTimeDist;
    private final String analyzer;

    public QsysMapDcResult(double meanQueueLength, double meanWaitingTime, double meanSojournTime,
                           double utilization, Matrix queueLengthDist, Matrix waitingTimeDist, String analyzer) {
        this.meanQueueLength = meanQueueLength;
        this.meanWaitingTime = meanWaitingTime;
        this.meanSojournTime = meanSojournTime;
        this.utilization = utilization;
        this.queueLengthDist = queueLengthDist;
        this.waitingTimeDist = waitingTimeDist;
        this.analyzer = analyzer;
    }

    public double getMeanQueueLength() { return meanQueueLength; }
    public double getMeanWaitingTime() { return meanWaitingTime; }
    public double getMeanSojournTime() { return meanSojournTime; }
    public double getUtilization() { return utilization; }
    public Matrix getQueueLengthDist() { return queueLengthDist; }
    public Matrix getWaitingTimeDist() { return waitingTimeDist; }
    public String getAnalyzer() { return analyzer; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof QsysMapDcResult)) return false;
        QsysMapDcResult that = (QsysMapDcResult) o;
        return Double.compare(that.meanQueueLength, meanQueueLength) == 0
                && Double.compare(that.meanWaitingTime, meanWaitingTime) == 0
                && Double.compare(that.meanSojournTime, meanSojournTime) == 0
                && Double.compare(that.utilization, utilization) == 0
                && Objects.equals(queueLengthDist, that.queueLengthDist)
                && Objects.equals(waitingTimeDist, that.waitingTimeDist)
                && Objects.equals(analyzer, that.analyzer);
    }

    @Override
    public int hashCode() {
        return Objects.hash(meanQueueLength, meanWaitingTime, meanSojournTime, utilization,
                queueLengthDist, waitingTimeDist, analyzer);
    }

    @Override
    public String toString() {
        return "QsysMapDcResult(meanQueueLength=" + meanQueueLength
                + ", meanWaitingTime=" + meanWaitingTime
                + ", meanSojournTime=" + meanSojournTime
                + ", utilization=" + utilization
                + ", queueLengthDist=" + queueLengthDist
                + ", waitingTimeDist=" + waitingTimeDist
                + ", analyzer=" + analyzer + ")";
    }
}
