/**
 * @file Result of MAP/PH type queue analysis
 *
 * @since LINE 3.1.0
 */
package jline.api.qsys;

import java.util.Objects;

import jline.util.matrix.Matrix;

/**
 * Result of MAP/PH type queue analysis.
 */
public final class QsysMapPhResult {
    private final double meanQueueLength;
    private final double meanWaitingTime;
    private final double meanSojournTime;
    private final double utilization;
    private final Matrix queueLengthDist;
    private final Matrix queueLengthMoments;
    private final Matrix sojournTimeMoments;
    private final String analyzer;

    public QsysMapPhResult(double meanQueueLength, double meanWaitingTime, double meanSojournTime,
                           double utilization, Matrix queueLengthDist, Matrix queueLengthMoments,
                           Matrix sojournTimeMoments, String analyzer) {
        this.meanQueueLength = meanQueueLength;
        this.meanWaitingTime = meanWaitingTime;
        this.meanSojournTime = meanSojournTime;
        this.utilization = utilization;
        this.queueLengthDist = queueLengthDist;
        this.queueLengthMoments = queueLengthMoments;
        this.sojournTimeMoments = sojournTimeMoments;
        this.analyzer = analyzer;
    }

    public double getMeanQueueLength() { return meanQueueLength; }
    public double getMeanWaitingTime() { return meanWaitingTime; }
    public double getMeanSojournTime() { return meanSojournTime; }
    public double getUtilization() { return utilization; }
    public Matrix getQueueLengthDist() { return queueLengthDist; }
    public Matrix getQueueLengthMoments() { return queueLengthMoments; }
    public Matrix getSojournTimeMoments() { return sojournTimeMoments; }
    public String getAnalyzer() { return analyzer; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof QsysMapPhResult)) return false;
        QsysMapPhResult that = (QsysMapPhResult) o;
        return Double.compare(that.meanQueueLength, meanQueueLength) == 0
                && Double.compare(that.meanWaitingTime, meanWaitingTime) == 0
                && Double.compare(that.meanSojournTime, meanSojournTime) == 0
                && Double.compare(that.utilization, utilization) == 0
                && Objects.equals(queueLengthDist, that.queueLengthDist)
                && Objects.equals(queueLengthMoments, that.queueLengthMoments)
                && Objects.equals(sojournTimeMoments, that.sojournTimeMoments)
                && Objects.equals(analyzer, that.analyzer);
    }

    @Override
    public int hashCode() {
        return Objects.hash(meanQueueLength, meanWaitingTime, meanSojournTime, utilization,
                queueLengthDist, queueLengthMoments, sojournTimeMoments, analyzer);
    }

    @Override
    public String toString() {
        return "QsysMapPhResult(meanQueueLength=" + meanQueueLength
                + ", meanWaitingTime=" + meanWaitingTime
                + ", meanSojournTime=" + meanSojournTime
                + ", utilization=" + utilization
                + ", queueLengthDist=" + queueLengthDist
                + ", queueLengthMoments=" + queueLengthMoments
                + ", sojournTimeMoments=" + sojournTimeMoments
                + ", analyzer=" + analyzer + ")";
    }
}
