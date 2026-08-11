/**
 * @file Result of RAP/RAP/1 QBD analysis
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import java.util.Objects;

import jline.util.matrix.Matrix;

public final class QbdRapRap1Result {
    private final double XN;
    private final double QN;
    private final double UN;
    private final Matrix pqueue;
    private final Matrix R;
    private final Matrix eta;
    private final Matrix G;
    private final Matrix B;
    private final Matrix L;
    private final Matrix F;

    public QbdRapRap1Result(double XN, double QN, double UN, Matrix pqueue, Matrix R, Matrix eta,
                            Matrix G, Matrix B, Matrix L, Matrix F) {
        this.XN = XN;
        this.QN = QN;
        this.UN = UN;
        this.pqueue = pqueue;
        this.R = R;
        this.eta = eta;
        this.G = G;
        this.B = B;
        this.L = L;
        this.F = F;
    }

    public double getXN() { return XN; }
    public double getQN() { return QN; }
    public double getUN() { return UN; }
    public Matrix getPqueue() { return pqueue; }
    public Matrix getR() { return R; }
    public Matrix getEta() { return eta; }
    public Matrix getG() { return G; }
    public Matrix getB() { return B; }
    public Matrix getL() { return L; }
    public Matrix getF() { return F; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof QbdRapRap1Result)) return false;
        QbdRapRap1Result that = (QbdRapRap1Result) o;
        return Double.compare(that.XN, XN) == 0
                && Double.compare(that.QN, QN) == 0
                && Double.compare(that.UN, UN) == 0
                && Objects.equals(pqueue, that.pqueue) && Objects.equals(R, that.R)
                && Objects.equals(eta, that.eta) && Objects.equals(G, that.G)
                && Objects.equals(B, that.B) && Objects.equals(L, that.L) && Objects.equals(F, that.F);
    }

    @Override
    public int hashCode() {
        return Objects.hash(XN, QN, UN, pqueue, R, eta, G, B, L, F);
    }

    @Override
    public String toString() {
        return "QbdRapRap1Result(XN=" + XN + ", QN=" + QN + ", UN=" + UN
                + ", pqueue=" + pqueue + ", R=" + R + ", eta=" + eta + ", G=" + G
                + ", B=" + B + ", L=" + L + ", F=" + F + ")";
    }
}
