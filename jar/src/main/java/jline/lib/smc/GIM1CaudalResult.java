/**
 * @file Result of GIM1 Caudal computation
 *
 * @since LINE 3.1.0
 */
package jline.lib.smc;

import java.util.Objects;

import jline.util.matrix.Matrix;

/**
 * Result of GIM1 Caudal computation.
 */
public final class GIM1CaudalResult {
    private final double eta;
    private final Matrix v;

    public GIM1CaudalResult(double eta, Matrix v) {
        this.eta = eta;
        this.v = v;
    }

    public double getEta() { return eta; }
    public Matrix getV() { return v; }

    public double component1() { return eta; }
    public Matrix component2() { return v; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof GIM1CaudalResult)) return false;
        GIM1CaudalResult that = (GIM1CaudalResult) o;
        return Double.compare(that.eta, eta) == 0 && Objects.equals(v, that.v);
    }

    @Override
    public int hashCode() {
        return Objects.hash(eta, v);
    }

    @Override
    public String toString() {
        return "GIM1CaudalResult(eta=" + eta + ", v=" + v + ")";
    }
}
