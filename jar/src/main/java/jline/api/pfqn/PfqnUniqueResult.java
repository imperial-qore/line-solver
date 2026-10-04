/**
 * @file PfqnUniqueResult data class
 *
 * @since LINE 3.0
 */
package jline.api.pfqn;

import java.util.Arrays;
import java.util.Objects;

import jline.util.matrix.Matrix;

/**
 * Result class for pfqn_unique containing all output matrices and mapping information.
 */
public final class PfqnUniqueResult {
    private final Matrix L_unique;
    private final Matrix mu_unique;
    private final Matrix gamma_unique;
    private final Matrix mi;
    private final int[] mapping;

    public PfqnUniqueResult(Matrix L_unique, Matrix mu_unique, Matrix gamma_unique, Matrix mi, int[] mapping) {
        this.L_unique = L_unique;
        this.mu_unique = mu_unique;
        this.gamma_unique = gamma_unique;
        this.mi = mi;
        this.mapping = mapping;
    }

    public Matrix getL_unique() { return L_unique; }
    public Matrix getMu_unique() { return mu_unique; }
    public Matrix getGamma_unique() { return gamma_unique; }
    public Matrix getMi() { return mi; }
    public int[] getMapping() { return mapping; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof PfqnUniqueResult)) return false;
        PfqnUniqueResult that = (PfqnUniqueResult) o;
        return Objects.equals(L_unique, that.L_unique)
                && Objects.equals(mu_unique, that.mu_unique)
                && Objects.equals(gamma_unique, that.gamma_unique)
                && Objects.equals(mi, that.mi)
                && Arrays.equals(mapping, that.mapping);
    }

    @Override
    public int hashCode() {
        int result = Objects.hash(L_unique, mu_unique, gamma_unique, mi);
        result = 31 * result + Arrays.hashCode(mapping);
        return result;
    }

    @Override
    public String toString() {
        return "PfqnUniqueResult(L_unique=" + L_unique + ", mu_unique=" + mu_unique
                + ", gamma_unique=" + gamma_unique + ", mi=" + mi
                + ", mapping=" + Arrays.toString(mapping) + ")";
    }
}
