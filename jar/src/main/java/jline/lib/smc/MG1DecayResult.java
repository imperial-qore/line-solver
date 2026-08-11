/**
 * @file Result of MG1 Decay computation
 *
 * @since LINE 3.1.0
 */
package jline.lib.smc;

import java.util.Objects;

import jline.util.matrix.Matrix;

/**
 * Result of MG1 Decay computation.
 */
public final class MG1DecayResult {
    private final double eta;
    private final Matrix uT;

    public MG1DecayResult(double eta, Matrix uT) {
        this.eta = eta;
        this.uT = uT;
    }

    public double getEta() { return eta; }
    public Matrix getUT() { return uT; }

    public double component1() { return eta; }
    public Matrix component2() { return uT; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof MG1DecayResult)) return false;
        MG1DecayResult that = (MG1DecayResult) o;
        return Double.compare(that.eta, eta) == 0 && Objects.equals(uT, that.uT);
    }

    @Override
    public int hashCode() {
        return Objects.hash(eta, uT);
    }

    @Override
    public String toString() {
        return "MG1DecayResult(eta=" + eta + ", uT=" + uT + ")";
    }
}
