/**
 * Kpctoolbox facade class providing static-style access to KPC-Toolbox functions.
 */
package jline.lib.kpctoolbox;

import java.util.ArrayList;
import java.util.List;

import jline.lib.kpctoolbox.basic.BasicUtils;
import jline.lib.kpctoolbox.kpcfit.KPCFit.KPCFitOptions;
import jline.lib.kpctoolbox.kpcfit.KPCFit.KPCFitResult;
import jline.lib.kpctoolbox.mc.CTMC;
import jline.lib.kpctoolbox.mc.DTMC;
import jline.lib.kpctoolbox.mmpp.MMPP;
import jline.lib.kpctoolbox.mvph.MVPH;
import jline.lib.kpctoolbox.trace.TraceAnalysis;
import jline.lib.kpctoolbox.trace.TraceAnalysis.TraceSummary;
import jline.util.Pair;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Kpctoolbox {

    private Kpctoolbox() {}

    // ========== TRACE ANALYSIS ==========

    public static double trace_mean(Matrix S) {
        return TraceAnalysis.trace_mean(S.toArray1D());
    }

    public static double trace_var(Matrix S) {
        return TraceAnalysis.trace_var(S.toArray1D());
    }

    public static double trace_scv(Matrix S) {
        return TraceAnalysis.trace_scv(S.toArray1D());
    }

    public static Matrix trace_acf(Matrix S, int maxLag) {
        double[] data = S.toArray1D();
        int[] lags = new int[maxLag];
        for (int i = 0; i < maxLag; i++) lags[i] = i + 1;
        double[] acf = TraceAnalysis.trace_acf(data, lags);
        Matrix result = new Matrix(acf.length, 1);
        for (int i = 0; i < acf.length; i++) {
            result.set(i, 0, acf[i]);
        }
        return result;
    }

    public static double trace_skew(Matrix S) {
        return TraceAnalysis.trace_skew(S.toArray1D());
    }

    public static double trace_joint(Matrix S, Matrix lags, Matrix orders) {
        double[] sArr = S.toArray1D();
        int[] lagsArr = new int[lags.getNumRows() * lags.getNumCols()];
        int[] ordersArr = new int[orders.getNumRows() * orders.getNumCols()];
        double[] lags1D = lags.toArray1D();
        double[] orders1D = orders.toArray1D();
        for (int i = 0; i < lagsArr.length; i++) {
            lagsArr[i] = (int) lags1D[i];
        }
        for (int i = 0; i < ordersArr.length; i++) {
            ordersArr[i] = (int) orders1D[i];
        }
        return TraceAnalysis.trace_joint(sArr, lagsArr, ordersArr);
    }

    public static Matrix trace_bicov(Matrix S1, Matrix S2, int maxLag) {
        double[] data = S1.toArray1D();
        int[] grid = new int[maxLag];
        for (int i = 0; i < maxLag; i++) grid[i] = i + 1;
        Pair<double[], int[][]> result = TraceAnalysis.trace_bicov(data, grid);
        double[] bicov = result.getFirst();
        Matrix m = new Matrix(bicov.length, 1);
        for (int i = 0; i < bicov.length; i++) {
            m.set(i, 0, bicov[i]);
        }
        return m;
    }

    public static Matrix trace_bicov(Matrix S, int[] grid) {
        Pair<double[], int[][]> result = TraceAnalysis.trace_bicov(S.toArray1D(), grid);
        double[] bicov = result.getFirst();
        Matrix m = new Matrix(bicov.length, 1);
        for (int i = 0; i < bicov.length; i++) {
            m.set(i, 0, bicov[i]);
        }
        return m;
    }

    public static double trace_idi(Matrix S, int numIntervals) {
        double[] data = S.toArray1D();
        double[] result = TraceAnalysis.trace_idi(data, new int[] { numIntervals });
        return (result.length > 0) ? result[0] : Double.NaN;
    }

    public static double trace_idc(Matrix S, double timeWindow) {
        double[] data = S.toArray1D();
        double mean = TraceAnalysis.trace_mean(data);
        int k;
        if (mean > 0) {
            k = Math.max(1, (int) Math.ceil(timeWindow / mean));
        } else {
            k = 1;
        }
        double[] result = TraceAnalysis.trace_idi(data, new int[] { k });
        return (result.length > 0) ? result[0] : Double.NaN;
    }

    public static Matrix trace_gamma(Matrix S, int limit) {
        double[] data = S.toArray1D();
        double gamma = TraceAnalysis.trace_gamma(data, limit);
        Matrix result = new Matrix(1, 1);
        result.set(0, 0, gamma);
        return result;
    }

    public static Matrix trace_shuffle(Matrix S) {
        double[] shuffled = TraceAnalysis.trace_shuffle(S.toArray1D());
        Matrix result = new Matrix(shuffled.length, 1);
        for (int i = 0; i < shuffled.length; i++) {
            result.set(i, 0, shuffled[i]);
        }
        return result;
    }

    public static Matrix trace_iat2counts(Matrix S, int numBins) {
        double[] data = S.toArray1D();
        double scale = TraceAnalysis.trace_mean(data) * numBins;
        int[] counts = TraceAnalysis.trace_iat2counts(data, scale);
        Matrix result = new Matrix(counts.length, 1);
        for (int i = 0; i < counts.length; i++) {
            result.set(i, 0, (double) counts[i]);
        }
        return result;
    }

    public static Matrix trace_iat2bins(Matrix S, int numBins) {
        double[] data = S.toArray1D();
        double totalTime = 0.0;
        for (double v : data) totalTime += v;
        double scale = totalTime / numBins;
        Pair<int[], int[]> result = TraceAnalysis.trace_iat2bins(data, scale);
        int[] counts = result.getFirst();
        Matrix m = new Matrix(counts.length, 1);
        for (int i = 0; i < counts.length; i++) {
            m.set(i, 0, (double) counts[i]);
        }
        return m;
    }

    public static Matrix trace_pmf(Matrix S) {
        Pair<double[], double[]> result = TraceAnalysis.trace_pmf(S.toArray1D());
        double[] pmf = result.getFirst();
        double[] values = result.getSecond();
        Matrix m = new Matrix(pmf.length, 2);
        for (int i = 0; i < pmf.length; i++) {
            m.set(i, 0, pmf[i]);
            m.set(i, 1, values[i]);
        }
        return m;
    }

    public static TraceSummary trace_summary(Matrix S) {
        return TraceAnalysis.trace_summary(S.toArray1D());
    }

    public static Matrix autocov(Matrix S) {
        double[] acv = TraceAnalysis.autocov(S.toArray1D());
        Matrix result = new Matrix(acv.length, 1);
        for (int i = 0; i < acv.length; i++) {
            result.set(i, 0, acv[i]);
        }
        return result;
    }

    // ========== MULTI-TRACE ANALYSIS ==========

    public static double mtrace_mean(java.util.List<Matrix> traces) {
        ArrayList<double[]> traceArrays = new ArrayList<double[]>();
        for (int i = 0; i < traces.size(); i++) {
            traceArrays.add(traces.get(i).toArray1D());
        }
        double[] means = TraceAnalysis.mtrace_mean(traceArrays);
        double sum = 0.0;
        for (double m : means) sum += m;
        return sum / means.length;
    }

    public static double mtrace_var(java.util.List<Matrix> traces) {
        ArrayList<double[]> traceArrays = new ArrayList<double[]>();
        for (int i = 0; i < traces.size(); i++) {
            traceArrays.add(traces.get(i).toArray1D());
        }
        double[] means = TraceAnalysis.mtrace_mean(traceArrays);
        double total = 0.0;
        for (double m : means) total += m;
        double overallMean = total / means.length;
        double sumSqDiff = 0.0;
        for (double m : means) {
            sumSqDiff += (m - overallMean) * (m - overallMean);
        }
        return (means.length > 1) ? sumSqDiff / (means.length - 1) : 0.0;
    }

    // ========== BASIC UTILITIES ==========

    public static Matrix eye(int n) {
        return Matrix.eye(n);
    }

    public static Matrix ones(int m, int n) {
        Matrix result = new Matrix(m, n);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                result.set(i, j, 1.0);
            }
        }
        return result;
    }

    public static Matrix zeros(int m, int n) {
        return new Matrix(m, n);
    }

    public static int maxpos(Matrix data) {
        return BasicUtils.maxpos(data.toArray1D());
    }

    public static int minpos(Matrix data) {
        return BasicUtils.minpos(data.toArray1D());
    }

    // ========== MVPH (MULTIVARIATE PHASE-TYPE) ==========

    public static double mvph_mean_x(Matrix alpha, Matrix A, Matrix B) {
        return MVPH.mvph_mean_x(alpha.toArray1D(), A, B, Matrix.eye(B.getNumRows()));
    }

    public static double mvph_mean_y(Matrix alpha, Matrix A, Matrix C) {
        return MVPH.mvph_mean_y(alpha.toArray1D(), A, C, Matrix.eye(A.getNumRows()));
    }

    public static double mvph_cov(Matrix alpha, Matrix A, Matrix B, Matrix C) {
        return MVPH.mvph_cov(alpha.toArray1D(), A, B, C);
    }

    public static double mvph_corr(Matrix alpha, Matrix A, Matrix B, Matrix C) {
        return MVPH.mvph_corr(alpha.toArray1D(), A, B, C);
    }

    public static double mvph_joint(Matrix alpha, Matrix A, Matrix B, Matrix C, Matrix x, Matrix y) {
        int n1 = (int) x.get(0, 0);
        int n2 = (int) y.get(0, 0);
        return MVPH.mvph_joint(alpha.toArray1D(), A, B, C, n1, n2);
    }
}
