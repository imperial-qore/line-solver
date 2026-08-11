/**
 * @file M3PP superposition trace-based counting process fitting
 *
 * Fits a marked MMAP to a multi-class trace by matching per-class counting-process
 * statistics, following the MATLAB reference {@code m3pp_superpos_fitc_trace.m}.
 *
 * @since LINE 3.0
 */
package jline.api.mam.m3pp;

import java.util.ArrayList;
import java.util.List;
import java.util.TreeSet;

import jline.api.trace.Mtrace_iat2counts;
import jline.util.Pair;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class M3pp_superpos_fitc_trace {
    private M3pp_superpos_fitc_trace() {}

    /**
     * Superposes k MMPP[2] processes to fit a multi-class trace with m classes.
     *
     * <p>Port of MATLAB {@code m3pp_superpos_fitc_trace.m}: per-class rates,
     * IDC(t)/IDC(inf) and third central moment of counts are estimated from the
     * per-class counting processes of the trace (via
     * {@link Mtrace_iat2counts#mtrace_iat2counts}) and fed to
     * {@link M3pp_superpos_fitc#m3pp_superpos_fitc}.</p>
     *
     * @param T    inter-arrival times
     * @param A    per-event class labels
     * @param t    finite time scale (null -> 10*mean(T))
     * @param tinf near-infinite time scale (null -> max(10*t,(sum(T)-T(0))/100))
     * @return pair of the fitted marked MMAP (flat block array) and the per-class components
     */
    public static Pair<Matrix[], List<Matrix[]>> m3pp_superpos_fitc_trace(double[] T, int[] A,
                                                                          Double t, Double tinf) {
        double meanT = mean(T);
        double sumT = sum(T);
        double tUsed = (t != null) ? t : 10.0 * meanT;
        double tinfUsed = (tinf != null) ? tinf : Math.max(10.0 * tUsed, (sumT - T[0]) / 100.0);

        double a = 1.0 / meanT;

        // per-class rates (classes in sorted label order, matching mtrace_iat2counts)
        TreeSet<Integer> distinct = new TreeSet<Integer>();
        for (int v : A) distinct.add(v);
        List<Integer> L = new ArrayList<Integer>(distinct);
        int m = L.size();
        double[] pv = new double[m];
        for (int i = 0; i < m; i++) {
            int count = 0;
            for (int v : A) if (v == L.get(i)) count++;
            pv[i] = (double) count / A.length;
        }
        double[] av = new double[m];
        for (int i = 0; i < m; i++) av[i] = pv[i] * a;

        // per-class counting processes at the two time scales
        Matrix Nt = Mtrace_iat2counts.mtrace_iat2counts(T, A, tUsed);
        Matrix Ninf = Mtrace_iat2counts.mtrace_iat2counts(T, A, tinfUsed);

        double[] btv = new double[m];
        double[] binfv = new double[m];
        double[] m3tv = new double[m];
        for (int i = 0; i < m; i++) {
            double[] col = column(Nt, i);
            double[] colInf = column(Ninf, i);
            // IDC = var(N)/(rate*t); MATLAB var normalizes by (N-1).
            btv[i] = sampleVariance(col) / (av[i] * tUsed);
            binfv[i] = sampleVariance(colInf) / (av[i] * tinfUsed);
            // third central moment of counts from raw moments (mean normalizes by N)
            double m1 = mean(col);
            double m2 = rawMoment(col, 2);
            double m3 = rawMoment(col, 3);
            m3tv[i] = m3 - 3.0 * m2 * m1 + 2.0 * m1 * m1 * m1;
        }

        Pair<MatrixCell, List<MatrixCell>> res =
                M3pp_superpos_fitc.m3pp_superpos_fitc(av, btv, binfv, m3tv, tUsed, tinfUsed);

        Matrix[] fit = toArray(res.getLeft());
        List<Matrix[]> comps = new ArrayList<Matrix[]>();
        for (MatrixCell comp : res.getRight()) comps.add(toArray(comp));
        return new Pair<Matrix[], List<Matrix[]>>(fit, comps);
    }

    public static Pair<Matrix[], List<Matrix[]>> m3pp_superpos_fitc_trace(double[] T, int[] A) {
        return m3pp_superpos_fitc_trace(T, A, null, null);
    }

    /** Matrix-argument overload (T and A as column matrices). */
    public static Pair<Matrix[], List<Matrix[]>> m3pp_superpos_fitc_trace(Matrix T, Matrix A,
                                                                          Double t, Double tinf) {
        double[] tArray = T.toArray1D();
        double[] aVals = A.toArray1D();
        int[] aArray = new int[aVals.length];
        for (int i = 0; i < aArray.length; i++) aArray[i] = (int) aVals[i];
        return m3pp_superpos_fitc_trace(tArray, aArray, t, tinf);
    }

    public static Pair<Matrix[], List<Matrix[]>> m3pp_superpos_fitc_trace(Matrix T, Matrix A) {
        return m3pp_superpos_fitc_trace(T, A, null, null);
    }

    private static Matrix[] toArray(MatrixCell cell) {
        Matrix[] out = new Matrix[cell.size()];
        for (int i = 0; i < cell.size(); i++) out[i] = cell.get(i);
        return out;
    }

    private static double[] column(Matrix M, int col) {
        int n = M.getNumRows();
        double[] out = new double[n];
        for (int r = 0; r < n; r++) out[r] = M.get(r, col);
        return out;
    }

    private static double sum(double[] arr) {
        double s = 0.0;
        for (double v : arr) s += v;
        return s;
    }

    private static double mean(double[] arr) {
        if (arr.length == 0) return 0.0;
        return sum(arr) / arr.length;
    }

    private static double rawMoment(double[] arr, int order) {
        if (arr.length == 0) return 0.0;
        double s = 0.0;
        for (double v : arr) s += Math.pow(v, order);
        return s / arr.length;
    }

    private static double sampleVariance(double[] arr) {
        int n = arr.length;
        if (n < 2) return 0.0;
        double mu = mean(arr);
        double s = 0.0;
        for (double v : arr) s += (v - mu) * (v - mu);
        return s / (n - 1);
    }
}
