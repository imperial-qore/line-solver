/**
 * @file Result of aoi_dist2ph conversion
 *
 * @since LINE 3.0
 */
package jline.api.aoi;

import java.util.Objects;

import jline.util.matrix.Matrix;

/**
 * Result of aoi_dist2ph conversion.
 */
public final class Aoi_dist2phResult {
    private final Matrix alpha;
    private final Matrix T;

    public Aoi_dist2phResult(Matrix alpha, Matrix T) {
        this.alpha = alpha;
        this.T = T;
    }

    public Matrix getAlpha() { return alpha; }
    public Matrix getT() { return T; }

    public Matrix component1() { return alpha; }
    public Matrix component2() { return T; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof Aoi_dist2phResult)) return false;
        Aoi_dist2phResult that = (Aoi_dist2phResult) o;
        return Objects.equals(alpha, that.alpha) && Objects.equals(T, that.T);
    }

    @Override
    public int hashCode() {
        return Objects.hash(alpha, T);
    }

    @Override
    public String toString() {
        return "Aoi_dist2phResult(alpha=" + alpha + ", T=" + T + ")";
    }
}
