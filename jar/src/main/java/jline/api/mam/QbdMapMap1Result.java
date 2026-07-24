/**
 * @file QBD MAP/MAP/1 result type
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import java.util.Objects;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * Result of qbd_mapmap1 analysis.
 */
public final class QbdMapMap1Result {
    private final double XN;
    private final double QN;
    private final double UN;
    private final Matrix pqueue;
    private final Matrix R;
    private final Matrix eta;
    private final Matrix G;
    private final Matrix A_1;
    private final Matrix A0;
    private final Matrix A1;
    private final Matrix U;
    private final MatrixCell MAPs;

    public QbdMapMap1Result(double XN, double QN, double UN, Matrix pqueue, Matrix R, Matrix eta,
                            Matrix G, Matrix A_1, Matrix A0, Matrix A1, Matrix U, MatrixCell MAPs) {
        this.XN = XN;
        this.QN = QN;
        this.UN = UN;
        this.pqueue = pqueue;
        this.R = R;
        this.eta = eta;
        this.G = G;
        this.A_1 = A_1;
        this.A0 = A0;
        this.A1 = A1;
        this.U = U;
        this.MAPs = MAPs;
    }

    public double getXN() { return XN; }
    public double getQN() { return QN; }
    public double getUN() { return UN; }
    public Matrix getPqueue() { return pqueue; }
    public Matrix getR() { return R; }
    public Matrix getEta() { return eta; }
    public Matrix getG() { return G; }
    public Matrix getA_1() { return A_1; }
    public Matrix getA0() { return A0; }
    public Matrix getA1() { return A1; }
    public Matrix getU() { return U; }
    public MatrixCell getMAPs() { return MAPs; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof QbdMapMap1Result)) return false;
        QbdMapMap1Result that = (QbdMapMap1Result) o;
        return Double.compare(that.XN, XN) == 0
                && Double.compare(that.QN, QN) == 0
                && Double.compare(that.UN, UN) == 0
                && Objects.equals(pqueue, that.pqueue)
                && Objects.equals(R, that.R)
                && Objects.equals(eta, that.eta)
                && Objects.equals(G, that.G)
                && Objects.equals(A_1, that.A_1)
                && Objects.equals(A0, that.A0)
                && Objects.equals(A1, that.A1)
                && Objects.equals(U, that.U)
                && Objects.equals(MAPs, that.MAPs);
    }

    @Override
    public int hashCode() {
        return Objects.hash(XN, QN, UN, pqueue, R, eta, G, A_1, A0, A1, U, MAPs);
    }

    @Override
    public String toString() {
        return "QbdMapMap1Result(XN=" + XN + ", QN=" + QN + ", UN=" + UN + ")";
    }
}
