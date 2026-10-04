package jline.lib.kpctoolbox.erchmm;

import jline.util.matrix.MatrixCell;

/** Result class for ER-CHMM EM fitting. */
public final class ERCHMMFitResult {
    public final MatrixCell MAP;
    public final double logLikelihood;
    public final int[] orders;

    public ERCHMMFitResult(MatrixCell MAP, double logLikelihood, int[] orders) {
        this.MAP = MAP;
        this.logLikelihood = logLikelihood;
        this.orders = orders;
    }

    public MatrixCell getMAP() { return MAP; }
    public double getLogLikelihood() { return logLikelihood; }
    public int[] getOrders() { return orders; }
}
