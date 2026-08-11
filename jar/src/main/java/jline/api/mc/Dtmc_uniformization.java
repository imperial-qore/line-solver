package jline.api.mc;

import java.util.Objects;

import jline.util.Pair;
import jline.util.matrix.Matrix;

public final class Dtmc_uniformization {
    private Dtmc_uniformization() {}

    /**
     * Result class for DTMC uniformization analysis containing the probability vector and maximum iterations used.
     */
    public static final class DtmcUniformizationResult {
        public final Matrix pi;
        public final int kmax;

        public DtmcUniformizationResult(Matrix pi, int kmax) {
            this.pi = pi;
            this.kmax = kmax;
        }

        public Matrix component1() { return pi; }
        public int component2() { return kmax; }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (!(o instanceof DtmcUniformizationResult)) return false;
            DtmcUniformizationResult that = (DtmcUniformizationResult) o;
            return kmax == that.kmax && Objects.equals(pi, that.pi);
        }

        @Override
        public int hashCode() {
            return Objects.hash(pi, kmax);
        }
    }

    public static DtmcUniformizationResult dtmc_uniformization(Matrix pi0, Matrix P) {
        return dtmc_uniformization(pi0, P, 1e4, 1e-12, -1);
    }

    public static DtmcUniformizationResult dtmc_uniformization(Matrix pi0, Matrix P, double t) {
        return dtmc_uniformization(pi0, P, t, 1e-12, -1);
    }

    public static DtmcUniformizationResult dtmc_uniformization(Matrix pi0, Matrix P, double t, double tol) {
        return dtmc_uniformization(pi0, P, t, tol, -1);
    }

    /**
     * Compute the transient probability distribution of a DTMC using uniformization.
     *
     * @param maxiter Cap on the truncation depth; pass a nonpositive value, as the
     *                overloads above do, to let the Fox-Glynn right truncation point
     *                size it so the Poisson tail is below tol rather than being cut
     *                at a fixed depth
     */
    public static DtmcUniformizationResult dtmc_uniformization(Matrix pi0, Matrix P, double t, double tol, int maxiter) {
        Matrix Q = Ctmc_makeinfgen.ctmc_makeinfgen(P);

        Pair<Matrix, Integer> result = ctmc_uniformization_with_kmax(pi0, Q, t, tol, maxiter);
        return new DtmcUniformizationResult(result.getLeft(), result.getRight());
    }

    /**
     * Extended CTMC uniformization that returns both the probability vector and kmax.
     */
    private static Pair<Matrix, Integer> ctmc_uniformization_with_kmax(Matrix pi0, Matrix Q, double t, double tol, int maxiter) {
        int n = Q.getNumCols();
        double q = 0.0;

        for (int i = 0; i < n; i++) {
            q = Math.max(q, 1.1 * Math.abs(Q.get(i, i)));
        }

        Matrix pi = Ctmc_foxglynn.ctmc_foxglynn(pi0, Q, t, tol, maxiter);
        int kmax = Ctmc_foxglynn.ctmc_foxglynn_weights(q * t, tol, maxiter).right;

        return new Pair<Matrix, Integer>(pi, kmax);
    }
}
