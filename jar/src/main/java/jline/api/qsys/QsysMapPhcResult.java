/**
 * @file Result of the exact MAP/PH/c analysis
 *
 * @since LINE 3.1.0
 */
package jline.api.qsys;

import jline.util.matrix.Matrix;

/** Return value of qsys_mapphc, mirroring the MATLAB struct. */
public final class QsysMapPhcResult {
    private final double meanQueueLength;
    private final double meanWaitingTime;
    private final double meanSojournTime;
    private final double utilization;
    private final Matrix queueLengthDist;
    private final Matrix waitingTimeMoments;
    private final Matrix waitingTimeCCDF;
    private final Matrix waitingTimePoints;
    private final double probWait;
    private final int phaseCount;
    private final String analyzer;

    public QsysMapPhcResult(double meanQueueLength, double meanWaitingTime, double meanSojournTime,
                            double utilization, Matrix queueLengthDist, Matrix waitingTimeMoments,
                            Matrix waitingTimeCCDF, Matrix waitingTimePoints, double probWait,
                            int phaseCount, String analyzer) {
        this.meanQueueLength = meanQueueLength;
        this.meanWaitingTime = meanWaitingTime;
        this.meanSojournTime = meanSojournTime;
        this.utilization = utilization;
        this.queueLengthDist = queueLengthDist;
        this.waitingTimeMoments = waitingTimeMoments;
        this.waitingTimeCCDF = waitingTimeCCDF;
        this.waitingTimePoints = waitingTimePoints;
        this.probWait = probWait;
        this.phaseCount = phaseCount;
        this.analyzer = analyzer;
    }

    public double getMeanQueueLength() { return meanQueueLength; }
    public double getMeanWaitingTime() { return meanWaitingTime; }
    public double getMeanSojournTime() { return meanSojournTime; }
    public double getUtilization() { return utilization; }
    public Matrix getQueueLengthDist() { return queueLengthDist; }
    public Matrix getWaitingTimeMoments() { return waitingTimeMoments; }
    public Matrix getWaitingTimeCCDF() { return waitingTimeCCDF; }
    public Matrix getWaitingTimePoints() { return waitingTimePoints; }
    public double getProbWait() { return probWait; }
    public int getPhaseCount() { return phaseCount; }
    public String getAnalyzer() { return analyzer; }
}
