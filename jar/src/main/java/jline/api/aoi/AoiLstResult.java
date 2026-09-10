/**
 * @file AoI LST result type
 *
 * Result of Age of Information analysis for general queue types.
 * Extends AoiResult with an optional LST (Laplace-Stieltjes Transform)
 * function for the AoI distribution.
 *
 * @since LINE 3.0
 */
package jline.api.aoi;

import java.util.Objects;

/**
 * Result of Age of Information analysis for general queue types.
 */
public final class AoiLstResult {
    private final double meanAoI;
    private final LstFunction lstAoI;
    private final double peakAoI;

    public AoiLstResult(double meanAoI, LstFunction lstAoI, double peakAoI) {
        this.meanAoI = meanAoI;
        this.lstAoI = lstAoI;
        this.peakAoI = peakAoI;
    }

    public double getMeanAoI() {
        return meanAoI;
    }

    public LstFunction getLstAoI() {
        return lstAoI;
    }

    public double getPeakAoI() {
        return peakAoI;
    }

    public double component1() {
        return meanAoI;
    }

    public LstFunction component2() {
        return lstAoI;
    }

    public double component3() {
        return peakAoI;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof AoiLstResult)) return false;
        AoiLstResult that = (AoiLstResult) o;
        return Double.compare(that.meanAoI, meanAoI) == 0
                && Double.compare(that.peakAoI, peakAoI) == 0
                && Objects.equals(lstAoI, that.lstAoI);
    }

    @Override
    public int hashCode() {
        return Objects.hash(meanAoI, lstAoI, peakAoI);
    }

    @Override
    public String toString() {
        return "AoiLstResult(meanAoI=" + meanAoI + ", lstAoI=" + lstAoI + ", peakAoI=" + peakAoI + ")";
    }
}
