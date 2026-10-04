/**
 * @file GK bound factor result for fork-join analysis
 *
 * @since LINE 3.0
 */
package jline.api.fj;

import java.util.Objects;

/**
 * Result container for G(K) bound factors.
 */
public final class GKBoundResult {
    public final int K;
    public final double exponential;
    public final double uniform;
    public final double evd;
    public final double upperBound;

    public GKBoundResult(int K, double exponential, double uniform, double evd, double upperBound) {
        this.K = K;
        this.exponential = exponential;
        this.uniform = uniform;
        this.evd = evd;
        this.upperBound = upperBound;
    }

    public int getK() { return K; }
    public double getExponential() { return exponential; }
    public double getUniform() { return uniform; }
    public double getEvd() { return evd; }
    public double getUpperBound() { return upperBound; }

    public int component1() { return K; }
    public double component2() { return exponential; }
    public double component3() { return uniform; }
    public double component4() { return evd; }
    public double component5() { return upperBound; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof GKBoundResult)) return false;
        GKBoundResult that = (GKBoundResult) o;
        return K == that.K
                && Double.compare(that.exponential, exponential) == 0
                && Double.compare(that.uniform, uniform) == 0
                && Double.compare(that.evd, evd) == 0
                && Double.compare(that.upperBound, upperBound) == 0;
    }

    @Override
    public int hashCode() {
        return Objects.hash(K, exponential, uniform, evd, upperBound);
    }

    @Override
    public String toString() {
        return "GKBoundResult(K=" + K + ", exponential=" + exponential
                + ", uniform=" + uniform + ", evd=" + evd
                + ", upperBound=" + upperBound + ")";
    }
}
