/**
 * @file Maximum Entropy OQN Algorithm Result
 *
 * Result of the ME OQN algorithm.
 *
 * @since LINE 3.0
 */
package jline.api.nc;

import java.util.Objects;

import jline.util.matrix.Matrix;

/**
 * Result of the ME OQN algorithm.
 */
public final class MeOqnResult {
    private final Matrix L;
    private final Matrix W;
    private final Matrix Ca;
    private final Matrix Cd;
    private final Matrix lambda;
    private final Matrix rho;
    private final int iter;

    public MeOqnResult(Matrix L, Matrix W, Matrix Ca, Matrix Cd, Matrix lambda, Matrix rho) {
        this(L, W, Ca, Cd, lambda, rho, 0);
    }

    public MeOqnResult(Matrix L, Matrix W, Matrix Ca, Matrix Cd, Matrix lambda, Matrix rho, int iter) {
        this.L = L;
        this.W = W;
        this.Ca = Ca;
        this.Cd = Cd;
        this.lambda = lambda;
        this.rho = rho;
        this.iter = iter;
    }

    public Matrix getL() { return L; }
    public Matrix getW() { return W; }
    public Matrix getCa() { return Ca; }
    public Matrix getCd() { return Cd; }
    public Matrix getLambda() { return lambda; }
    public Matrix getRho() { return rho; }
    public int getIter() { return iter; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof MeOqnResult)) return false;
        MeOqnResult that = (MeOqnResult) o;
        return iter == that.iter
                && Objects.equals(L, that.L)
                && Objects.equals(W, that.W)
                && Objects.equals(Ca, that.Ca)
                && Objects.equals(Cd, that.Cd)
                && Objects.equals(lambda, that.lambda)
                && Objects.equals(rho, that.rho);
    }

    @Override
    public int hashCode() {
        return Objects.hash(L, W, Ca, Cd, lambda, rho, iter);
    }

    @Override
    public String toString() {
        return "MeOqnResult(L=" + L + ", W=" + W + ", Ca=" + Ca
                + ", Cd=" + Cd + ", lambda=" + lambda + ", rho=" + rho
                + ", iter=" + iter + ")";
    }
}
