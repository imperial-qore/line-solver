/**
 * @file Result of Fork-Join MAM analysis
 *
 * @since LINE 3.0
 */
package jline.solvers.mam;

import java.util.List;
import java.util.Objects;

import jline.lib.fjcodes.MainFJ.FJPercentileResult;
import jline.util.matrix.Matrix;

/**
 * Result from Fork-Join analysis.
 */
public final class MAMFJResult {
    private final Matrix QN;
    private final Matrix UN;
    private final Matrix RN;
    private final Matrix TN;
    private final Matrix XN;
    private final List<FJPercentileResult> percentileResults;

    public MAMFJResult(Matrix QN, Matrix UN, Matrix RN, Matrix TN, Matrix XN,
                       List<FJPercentileResult> percentileResults) {
        this.QN = QN;
        this.UN = UN;
        this.RN = RN;
        this.TN = TN;
        this.XN = XN;
        this.percentileResults = percentileResults;
    }

    public Matrix getQN() { return QN; }
    public Matrix getUN() { return UN; }
    public Matrix getRN() { return RN; }
    public Matrix getTN() { return TN; }
    public Matrix getXN() { return XN; }
    public List<FJPercentileResult> getPercentileResults() { return percentileResults; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof MAMFJResult)) return false;
        MAMFJResult that = (MAMFJResult) o;
        return Objects.equals(QN, that.QN) && Objects.equals(UN, that.UN)
                && Objects.equals(RN, that.RN) && Objects.equals(TN, that.TN)
                && Objects.equals(XN, that.XN)
                && Objects.equals(percentileResults, that.percentileResults);
    }

    @Override
    public int hashCode() {
        return Objects.hash(QN, UN, RN, TN, XN, percentileResults);
    }

    @Override
    public String toString() {
        return "MAMFJResult(QN=" + QN + ", UN=" + UN + ", RN=" + RN
                + ", TN=" + TN + ", XN=" + XN + ", percentileResults=" + percentileResults + ")";
    }
}
