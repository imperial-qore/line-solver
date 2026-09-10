/**
 * @file SingleClassCountStats data class
 *
 * @since LINE 3.0
 */
package jline.api.trace;

import java.util.Arrays;
import java.util.Objects;

public final class SingleClassCountStats {
    private final double mean;
    private final double variance;
    private final double idc;
    private final double skewness;
    private final int[] counts;

    public SingleClassCountStats(double mean, double variance, double idc, double skewness, int[] counts) {
        this.mean = mean;
        this.variance = variance;
        this.idc = idc;
        this.skewness = skewness;
        this.counts = counts;
    }

    public double getMean() { return mean; }
    public double getVariance() { return variance; }
    public double getIdc() { return idc; }
    public double getSkewness() { return skewness; }
    public int[] getCounts() { return counts; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof SingleClassCountStats)) return false;
        SingleClassCountStats that = (SingleClassCountStats) o;
        return Double.compare(that.mean, mean) == 0
                && Double.compare(that.variance, variance) == 0
                && Double.compare(that.idc, idc) == 0
                && Double.compare(that.skewness, skewness) == 0
                && Arrays.equals(counts, that.counts);
    }

    @Override
    public int hashCode() {
        int result = Objects.hash(mean, variance, idc, skewness);
        result = 31 * result + Arrays.hashCode(counts);
        return result;
    }

    @Override
    public String toString() {
        return "SingleClassCountStats(mean=" + mean + ", variance=" + variance + ", idc=" + idc
                + ", skewness=" + skewness + ", counts=" + Arrays.toString(counts) + ")";
    }
}
