/**
 * @file Result container for fork-join characteristic maximum
 *
 * @since LINE 3.0
 */
package jline.api.fj;

import java.util.Objects;

/**
 * Result container for characteristic maximum.
 */
public final class CharMaxResult {
    public final double MK;
    public final double mK;

    public CharMaxResult(double MK, double mK) {
        this.MK = MK;
        this.mK = mK;
    }

    public double getMK() { return MK; }
    public double getMk() { return mK; }

    public double component1() { return MK; }
    public double component2() { return mK; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof CharMaxResult)) return false;
        CharMaxResult that = (CharMaxResult) o;
        return Double.compare(that.MK, MK) == 0
                && Double.compare(that.mK, mK) == 0;
    }

    @Override
    public int hashCode() {
        return Objects.hash(MK, mK);
    }

    @Override
    public String toString() {
        return "CharMaxResult(MK=" + MK + ", mK=" + mK + ")";
    }
}
