/**
 * @file Result of BMAP/MAP/1 queue analysis.
 *
 * @since LINE 3.1.0
 */
package jline.solvers.mam.handlers;

import java.util.Objects;

import jline.util.matrix.Matrix;

/**
 * Result of BMAP/MAP/1 queue analysis.
 */
public final class BMAPMAP1Result {
    /** Mean queue length E[N] */
    public final double meanQueueLength;
    /** Server utilization rho */
    public final double utilization;
    /** Mean response time E[R] */
    public final double meanResponseTime;
    /** Throughput (total customer arrival rate) */
    public final double throughput;
    /** Aggregated stationary probabilities [pi0, pi1, piStar] */
    public final Matrix pi;
    /** G matrix */
    public final Matrix G;
    /** Mean batch size */
    public final double meanBatchSize;

    public BMAPMAP1Result(double meanQueueLength, double utilization, double meanResponseTime,
                          double throughput, Matrix pi, Matrix G, double meanBatchSize) {
        this.meanQueueLength = meanQueueLength;
        this.utilization = utilization;
        this.meanResponseTime = meanResponseTime;
        this.throughput = throughput;
        this.pi = pi;
        this.G = G;
        this.meanBatchSize = meanBatchSize;
    }

    public double getMeanQueueLength() { return meanQueueLength; }
    public double getUtilization() { return utilization; }
    public double getMeanResponseTime() { return meanResponseTime; }
    public double getThroughput() { return throughput; }
    public Matrix getPi() { return pi; }
    public Matrix getG() { return G; }
    public double getMeanBatchSize() { return meanBatchSize; }

    public double component1() { return meanQueueLength; }
    public double component2() { return utilization; }
    public double component3() { return meanResponseTime; }
    public double component4() { return throughput; }
    public Matrix component5() { return pi; }
    public Matrix component6() { return G; }
    public double component7() { return meanBatchSize; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof BMAPMAP1Result)) return false;
        BMAPMAP1Result that = (BMAPMAP1Result) o;
        return Double.compare(that.meanQueueLength, meanQueueLength) == 0
                && Double.compare(that.utilization, utilization) == 0
                && Double.compare(that.meanResponseTime, meanResponseTime) == 0
                && Double.compare(that.throughput, throughput) == 0
                && Double.compare(that.meanBatchSize, meanBatchSize) == 0
                && Objects.equals(pi, that.pi)
                && Objects.equals(G, that.G);
    }

    @Override
    public int hashCode() {
        return Objects.hash(meanQueueLength, utilization, meanResponseTime, throughput, pi, G, meanBatchSize);
    }

    @Override
    public String toString() {
        return "BMAPMAP1Result(meanQueueLength=" + meanQueueLength + ", utilization=" + utilization
                + ", meanResponseTime=" + meanResponseTime + ", throughput=" + throughput
                + ", pi=" + pi + ", G=" + G + ", meanBatchSize=" + meanBatchSize + ")";
    }
}
