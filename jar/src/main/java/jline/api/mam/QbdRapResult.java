/**
 * @file Result of the equilibrium analysis of a QBD with RAP components
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import java.util.Objects;

import jline.util.matrix.Matrix;

/**
 * Return set of {@link Qbd_rap#qbd_rap}. The fields mirror the MATLAB return
 * list of qbd_rap.m and the Python QbdRapResult dataclass.
 */
public final class QbdRapResult {
    private final Matrix levelProb;
    private final double QN;
    private final Matrix R;
    private final Matrix G;
    private final Matrix U;
    private final double spr;
    private final Matrix pqueue;
    private final Matrix pi0;

    public QbdRapResult(Matrix levelProb, double QN, Matrix R, Matrix G, Matrix U, double spr, Matrix pqueue,
                        Matrix pi0) {
        this.levelProb = levelProb;
        this.QN = QN;
        this.R = R;
        this.G = G;
        this.U = U;
        this.spr = spr;
        this.pqueue = pqueue;
        this.pi0 = pi0;
    }

    /** Row vector of marginal level probabilities, levels 0 to numLevels. */
    public Matrix getLevelProb() { return levelProb; }

    /** Mean queue length, computed exactly as pi0*R*inv(I-R)^2*e. */
    public double getQN() { return QN; }

    /** Rate matrix R = A0*inv(-U). */
    public Matrix getR() { return R; }

    /** Matrix G solving A0*G^2 + A1*G + A2 = 0. */
    public Matrix getG() { return G; }

    /** Matrix U = A1 + A0*G. */
    public Matrix getU() { return U; }

    /** Spectral radius Sp(R); the process is positive recurrent iff Sp(R) &lt; 1. */
    public double getSpr() { return spr; }

    /** (numLevels+1) x m matrix whose n-th row is the level vector pi_n. */
    public Matrix getPqueue() { return pqueue; }

    /** Level-0 vector pi_0, the boundary vector of Theorem 7. */
    public Matrix getPi0() { return pi0; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof QbdRapResult)) return false;
        QbdRapResult that = (QbdRapResult) o;
        return Double.compare(that.QN, QN) == 0
                && Double.compare(that.spr, spr) == 0
                && Objects.equals(levelProb, that.levelProb) && Objects.equals(R, that.R)
                && Objects.equals(G, that.G) && Objects.equals(U, that.U)
                && Objects.equals(pqueue, that.pqueue) && Objects.equals(pi0, that.pi0);
    }

    @Override
    public int hashCode() {
        return Objects.hash(levelProb, QN, R, G, U, spr, pqueue, pi0);
    }

    @Override
    public String toString() {
        return "QbdRapResult(levelProb=" + levelProb + ", QN=" + QN + ", R=" + R
                + ", G=" + G + ", U=" + U + ", spr=" + spr + ", pqueue=" + pqueue + ", pi0=" + pi0 + ")";
    }
}
