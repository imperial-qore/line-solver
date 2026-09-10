/**
 * @file Result of MFQ-based AoI solver
 *
 * Contains ME parameters and moments for AoI/PAoI distributions.
 *
 * The AoI/PAoI distributions have CDF: F(t) = 1 - g * expm(A*t) * h
 *
 * @since LINE 3.0
 */
package jline.api.aoi;

import java.util.Objects;

import jline.util.matrix.Matrix;

public final class AoiMfqResult {
    private final Matrix aoiG;
    private final Matrix aoiA;
    private final Matrix aoiH;
    private final double aoiMean;
    private final double aoiVar;
    private final Matrix paoiG;
    private final Matrix paoiA;
    private final Matrix paoiH;
    private final double paoiMean;
    private final double paoiVar;
    private final String systemType;
    private final double preemption;

    public AoiMfqResult(Matrix aoiG, Matrix aoiA, Matrix aoiH, double aoiMean, double aoiVar,
                        Matrix paoiG, Matrix paoiA, Matrix paoiH, double paoiMean, double paoiVar,
                        String systemType, double preemption) {
        this.aoiG = aoiG;
        this.aoiA = aoiA;
        this.aoiH = aoiH;
        this.aoiMean = aoiMean;
        this.aoiVar = aoiVar;
        this.paoiG = paoiG;
        this.paoiA = paoiA;
        this.paoiH = paoiH;
        this.paoiMean = paoiMean;
        this.paoiVar = paoiVar;
        this.systemType = systemType;
        this.preemption = preemption;
    }

    public Matrix getAoiG() { return aoiG; }
    public Matrix getAoiA() { return aoiA; }
    public Matrix getAoiH() { return aoiH; }
    public double getAoiMean() { return aoiMean; }
    public double getAoiVar() { return aoiVar; }
    public Matrix getPaoiG() { return paoiG; }
    public Matrix getPaoiA() { return paoiA; }
    public Matrix getPaoiH() { return paoiH; }
    public double getPaoiMean() { return paoiMean; }
    public double getPaoiVar() { return paoiVar; }
    public String getSystemType() { return systemType; }
    public double getPreemption() { return preemption; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof AoiMfqResult)) return false;
        AoiMfqResult that = (AoiMfqResult) o;
        return Double.compare(that.aoiMean, aoiMean) == 0
                && Double.compare(that.aoiVar, aoiVar) == 0
                && Double.compare(that.paoiMean, paoiMean) == 0
                && Double.compare(that.paoiVar, paoiVar) == 0
                && Double.compare(that.preemption, preemption) == 0
                && Objects.equals(aoiG, that.aoiG)
                && Objects.equals(aoiA, that.aoiA)
                && Objects.equals(aoiH, that.aoiH)
                && Objects.equals(paoiG, that.paoiG)
                && Objects.equals(paoiA, that.paoiA)
                && Objects.equals(paoiH, that.paoiH)
                && Objects.equals(systemType, that.systemType);
    }

    @Override
    public int hashCode() {
        return Objects.hash(aoiG, aoiA, aoiH, aoiMean, aoiVar, paoiG, paoiA, paoiH, paoiMean, paoiVar, systemType, preemption);
    }

    @Override
    public String toString() {
        return "AoiMfqResult(aoiMean=" + aoiMean + ", aoiVar=" + aoiVar
                + ", paoiMean=" + paoiMean + ", paoiVar=" + paoiVar
                + ", systemType=" + systemType + ", preemption=" + preemption + ")";
    }
}
