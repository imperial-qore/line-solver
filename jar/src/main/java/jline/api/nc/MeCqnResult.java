/**
 * @file Maximum Entropy CQN Algorithm Result
 *
 * Result of the ME closed queueing network algorithm.
 *
 * @since LINE 3.0
 */
package jline.api.nc;

import jline.util.matrix.Matrix;

/**
 * Result of the ME closed queueing network algorithm (Kouvatsos 1994,
 * Section 3.3).
 */
public final class MeCqnResult {
    private final Matrix L;
    private final Matrix W;
    private final Matrix Ca;
    private final Matrix Cd;
    private final Matrix lambda;
    private final Matrix rho;
    private final Matrix X;
    private final int iter;

    public MeCqnResult(Matrix L, Matrix W, Matrix Ca, Matrix Cd, Matrix lambda,
                       Matrix rho, Matrix X, int iter) {
        this.L = L;
        this.W = W;
        this.Ca = Ca;
        this.Cd = Cd;
        this.lambda = lambda;
        this.rho = rho;
        this.X = X;
        this.iter = iter;
    }

    public Matrix getL() { return L; }
    public Matrix getW() { return W; }
    public Matrix getCa() { return Ca; }
    public Matrix getCd() { return Cd; }
    public Matrix getLambda() { return lambda; }
    public Matrix getRho() { return rho; }
    public Matrix getX() { return X; }
    public int getIter() { return iter; }
}
