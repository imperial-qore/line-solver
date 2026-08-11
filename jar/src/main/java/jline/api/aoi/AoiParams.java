/**
 * @file Extracted parameters for AoI analysis
 *
 * @since LINE 3.0
 */
package jline.api.aoi;

import java.util.Objects;

import jline.util.matrix.Matrix;

/**
 * Extracted parameters for AoI analysis.
 */
public final class AoiParams {
    private final Matrix tau;
    private final Matrix T;
    private final Matrix sigma;
    private final Matrix S;
    private final double p;
    private final double lambda;
    private final double r;
    private final String systemType;
    private final String arrivalType;

    public AoiParams(Matrix tau, Matrix T, Matrix sigma, Matrix S, double p,
                     double lambda, double r, String systemType, String arrivalType) {
        this.tau = tau;
        this.T = T;
        this.sigma = sigma;
        this.S = S;
        this.p = p;
        this.lambda = lambda;
        this.r = r;
        this.systemType = systemType;
        this.arrivalType = arrivalType;
    }

    public Matrix getTau() { return tau; }
    public Matrix getT() { return T; }
    public Matrix getSigma() { return sigma; }
    public Matrix getS() { return S; }
    public double getP() { return p; }
    public double getLambda() { return lambda; }
    public double getR() { return r; }
    public String getSystemType() { return systemType; }
    public String getArrivalType() { return arrivalType; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof AoiParams)) return false;
        AoiParams that = (AoiParams) o;
        return Double.compare(that.p, p) == 0
                && Double.compare(that.lambda, lambda) == 0
                && Double.compare(that.r, r) == 0
                && Objects.equals(tau, that.tau)
                && Objects.equals(T, that.T)
                && Objects.equals(sigma, that.sigma)
                && Objects.equals(S, that.S)
                && Objects.equals(systemType, that.systemType)
                && Objects.equals(arrivalType, that.arrivalType);
    }

    @Override
    public int hashCode() {
        return Objects.hash(tau, T, sigma, S, p, lambda, r, systemType, arrivalType);
    }

    @Override
    public String toString() {
        return "AoiParams(tau=" + tau + ", T=" + T + ", sigma=" + sigma + ", S=" + S
                + ", p=" + p + ", lambda=" + lambda + ", r=" + r
                + ", systemType=" + systemType + ", arrivalType=" + arrivalType + ")";
    }
}
