/**
 * @file AoI result type
 *
 * Result of Age of Information analysis for simple queue types.
 * Contains the three primary AoI metrics: mean, variance, and peak.
 *
 * @since LINE 3.0
 */
package jline.api.aoi;

import java.util.Objects;

/**
 * Result of Age of Information analysis for simple queue types.
 */
public final class AoiResult {
    private final double meanAoI;
    private final double varAoI;
    private final double peakAoI;

    public AoiResult(double meanAoI, double varAoI, double peakAoI) {
        this.meanAoI = meanAoI;
        this.varAoI = varAoI;
        this.peakAoI = peakAoI;
    }

    public double getMeanAoI() {
        return meanAoI;
    }

    public double getVarAoI() {
        return varAoI;
    }

    public double getPeakAoI() {
        return peakAoI;
    }

    public double component1() {
        return meanAoI;
    }

    public double component2() {
        return varAoI;
    }

    public double component3() {
        return peakAoI;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof AoiResult)) return false;
        AoiResult that = (AoiResult) o;
        return Double.compare(that.meanAoI, meanAoI) == 0
                && Double.compare(that.varAoI, varAoI) == 0
                && Double.compare(that.peakAoI, peakAoI) == 0;
    }

    @Override
    public int hashCode() {
        return Objects.hash(meanAoI, varAoI, peakAoI);
    }

    @Override
    public String toString() {
        return "AoiResult(meanAoI=" + meanAoI + ", varAoI=" + varAoI + ", peakAoI=" + peakAoI + ")";
    }
}
