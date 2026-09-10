/**
 * @file Result of the MMAP[K]/G[K]/1 per-type waiting time analysis
 *
 * @since LINE 3.1.0
 */
package jline.api.qsys;

import jline.util.matrix.Matrix;

/** Return value of qsys_mmapgk1, mirroring the MATLAB struct. */
public final class QsysMmapGk1Result {
    private final Matrix lambdas;
    private final double arrivalRate;
    private final double utilization;
    private final Matrix idleVector;
    private final Matrix waitMoments;
    private final Matrix meanWaitingTime;
    private final Matrix meanSojournTime;
    private final double meanQueueLength;
    private final Matrix waitCDF;
    private final Matrix waitPoints;
    private final String analyzer;

    public QsysMmapGk1Result(Matrix lambdas, double arrivalRate, double utilization,
                             Matrix idleVector, Matrix waitMoments, Matrix meanWaitingTime,
                             Matrix meanSojournTime, double meanQueueLength, Matrix waitCDF,
                             Matrix waitPoints, String analyzer) {
        this.lambdas = lambdas;
        this.arrivalRate = arrivalRate;
        this.utilization = utilization;
        this.idleVector = idleVector;
        this.waitMoments = waitMoments;
        this.meanWaitingTime = meanWaitingTime;
        this.meanSojournTime = meanSojournTime;
        this.meanQueueLength = meanQueueLength;
        this.waitCDF = waitCDF;
        this.waitPoints = waitPoints;
        this.analyzer = analyzer;
    }

    public Matrix getLambdas() { return lambdas; }
    public double getArrivalRate() { return arrivalRate; }
    public double getUtilization() { return utilization; }
    public Matrix getIdleVector() { return idleVector; }
    public Matrix getWaitMoments() { return waitMoments; }
    public Matrix getMeanWaitingTime() { return meanWaitingTime; }
    public Matrix getMeanSojournTime() { return meanSojournTime; }
    public double getMeanQueueLength() { return meanQueueLength; }
    public Matrix getWaitCDF() { return waitCDF; }
    public Matrix getWaitPoints() { return waitPoints; }
    public String getAnalyzer() { return analyzer; }
}
