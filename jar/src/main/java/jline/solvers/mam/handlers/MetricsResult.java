package jline.solvers.mam.handlers;

import jline.util.matrix.Matrix;

public final class MetricsResult {
    public final Matrix QN;
    public final Matrix UN;
    public final Matrix RN;
    public final Matrix TN;
    public final Matrix CN;
    public final Matrix XN;

    public MetricsResult(Matrix QN, Matrix UN, Matrix RN, Matrix TN, Matrix CN, Matrix XN) {
        this.QN = QN;
        this.UN = UN;
        this.RN = RN;
        this.TN = TN;
        this.CN = CN;
        this.XN = XN;
    }

    public Matrix getQN() { return QN; }
    public Matrix getUN() { return UN; }
    public Matrix getRN() { return RN; }
    public Matrix getTN() { return TN; }
    public Matrix getCN() { return CN; }
    public Matrix getXN() { return XN; }
}
