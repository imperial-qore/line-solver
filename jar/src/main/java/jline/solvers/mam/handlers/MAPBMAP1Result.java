package jline.solvers.mam.handlers;

import jline.util.matrix.Matrix;

/** Result of MAP/BMAP/1 queue analysis. */
public final class MAPBMAP1Result {
    public final double meanQueueLength;
    public final double utilization;
    public final double meanResponseTime;
    public final double throughput;
    public final Matrix pi;
    public final Matrix R;
    public final double meanBatchSize;

    public MAPBMAP1Result(double meanQueueLength, double utilization, double meanResponseTime,
                          double throughput, Matrix pi, Matrix R, double meanBatchSize) {
        this.meanQueueLength = meanQueueLength;
        this.utilization = utilization;
        this.meanResponseTime = meanResponseTime;
        this.throughput = throughput;
        this.pi = pi;
        this.R = R;
        this.meanBatchSize = meanBatchSize;
    }
}
