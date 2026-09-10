/**
 * @file M3PP theoretical counting process fitting
 *
 * @since LINE 3.0
 */
package jline.api.mam.m3pp;

import jline.util.Pair;

import jline.api.mam.Map_count_mean;
import jline.api.mam.Map_count_moment;
import jline.api.mam.Map_count_var;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class M3pp2m_fitc_theoretical {
    private M3pp2m_fitc_theoretical() {}

    /**
     * Fits the theoretical characteristics of a MMAP(n,m) with a M3PP(2,m).
     */
    public static Matrix[] m3pp2m_fitc_theoretical(MatrixCell mmap, String method, double t, double tinf) {
        int m = mmap.size() - 2;

        if ("approx_cov".equals(method) && m > 2) {
            throw new IllegalArgumentException("Approximate covariance fitting only supported for two classes.");
        }

        double t1 = t;
        double t2 = t;
        double t3 = t;

        MatrixCell aggregateMap = new MatrixCell(2);
        aggregateMap.set(0, mmap.get(0));
        aggregateMap.set(1, mmap.get(1));

        double countMean = Map_count_mean.map_count_mean(aggregateMap, t1);
        double a = countMean / t1;

        double countVar1 = Map_count_var.map_count_var(aggregateMap, t1);
        double countVar2 = Map_count_var.map_count_var(aggregateMap, t2);
        double countVarInf = Map_count_var.map_count_var(aggregateMap, tinf);

        double bt1 = countVar1 / (a * t1);
        double bt2 = countVar2 / (a * t2);
        double binf = countVarInf / (a * tinf);

        double mt2 = Map_count_moment.map_count_moment(aggregateMap, t2, 3);
        double m1t2 = Map_count_mean.map_count_mean(aggregateMap, t2);
        double m2t2 = Map_count_var.map_count_var(aggregateMap, t2) + m1t2 * m1t2;
        double m3t2 = mt2 - 3.0 * m2t2 * m1t2 + 2.0 * m1t2 * m1t2 * m1t2;

        double[] ai = new double[m];
        for (int i = 0; i < m; i++) {
            MatrixCell singleClassMap = new MatrixCell(2);
            singleClassMap.set(0, mmap.get(0));
            singleClassMap.set(1, mmap.get(2 + i));
            ai[i] = Map_count_mean.map_count_mean(singleClassMap, 1.0);
        }

        String lower = method.toLowerCase();
        if ("exact_delta".equals(lower)) {
            double[] dvt3 = computeVarianceDifferenceTheoretical(mmap, m, t3);
            return M3pp2m_fitc.m3pp2m_fitc(a, bt1, bt2, binf, m3t2, t1, t2, ai, dvt3, t3);
        } else if ("approx_delta".equals(lower)) {
            double[] dvt3 = computeVarianceDifferenceTheoretical(mmap, m, t3);
            return M3pp2m_fitc_approx.m3pp2m_fitc_approx(a, bt1, bt2, binf, m3t2, t1, t2, ai, dvt3, t3);
        } else if ("approx_cov".equals(lower)) {
            Pair<double[], Double> covRes = computeCovarianceTheoretical(mmap, m, t3);
            double[] vi = covRes.getFirst();
            double s = covRes.getSecond();
            return M3pp22_fitc_approx_cov.m3pp22_fitc_approx_cov(a, bt1, bt2, binf, m3t2, t1, t2, ai, s, t3);
        } else if ("approx_ag".equals(lower)) {
            double[] gt3 = computeVarianceCovarianceTheoretical(mmap, m, t3);
            return M3pp2m_fitc_approx_ag.m3pp2m_fitc_approx_ag(a, bt1, bt2, binf, m3t2, t1, t2, ai, gt3, t3);
        } else {
            throw new IllegalArgumentException("Invalid method '" + method + "'");
        }
    }

    public static Matrix[] m3pp2m_fitc_theoretical(MatrixCell mmap, String method, double t) {
        return m3pp2m_fitc_theoretical(mmap, method, t, 1e4);
    }

    public static Matrix[] m3pp2m_fitc_theoretical(MatrixCell mmap, String method) {
        return m3pp2m_fitc_theoretical(mmap, method, 1.0, 1e4);
    }

    public static Matrix[] m3pp2m_fitc_theoretical(MatrixCell mmap) {
        return m3pp2m_fitc_theoretical(mmap, "approx_delta", 1.0, 1e4);
    }

    /**
     * Computes per-class variance difference for delta fitting methods.
     */
    private static double[] computeVarianceDifferenceTheoretical(MatrixCell mmap, int m, double t) {
        double[] dvt = new double[m];
        for (int i = 0; i < m; i++) {
            MatrixCell mmap2 = new MatrixCell(4);
            mmap2.set(0, mmap.get(0));
            mmap2.set(1, mmap.get(1));
            mmap2.set(2, mmap.get(2 + i));

            Matrix otherClasses = Matrix.zeros(mmap.get(0).getNumRows(), mmap.get(0).getNumCols());
            for (int j = 0; j < m; j++) {
                if (j != i) {
                    otherClasses.addEq(mmap.get(2 + j));
                }
            }
            mmap2.set(3, otherClasses);

            MatrixCell singleClass = new MatrixCell(2);
            singleClass.set(0, mmap2.get(0));
            singleClass.set(1, mmap2.get(2));
            double var1 = Map_count_var.map_count_var(singleClass, t);

            MatrixCell otherClass = new MatrixCell(2);
            otherClass.set(0, mmap2.get(0));
            otherClass.set(1, mmap2.get(3));
            double var2 = Map_count_var.map_count_var(otherClass, t);

            dvt[i] = var1 - var2;
        }
        return dvt;
    }

    /**
     * Computes covariance for 2-class fitting.
     */
    private static Pair<double[], Double> computeCovarianceTheoretical(MatrixCell mmap, int m, double t) {
        double[] vi = new double[m];

        MatrixCell aggregateMap = new MatrixCell(2);
        aggregateMap.set(0, mmap.get(0));
        aggregateMap.set(1, mmap.get(1));
        double totalVar = Map_count_var.map_count_var(aggregateMap, t);

        for (int i = 0; i < m; i++) {
            MatrixCell singleClass = new MatrixCell(2);
            singleClass.set(0, mmap.get(0));
            singleClass.set(1, mmap.get(2 + i));
            vi[i] = Map_count_var.map_count_var(singleClass, t);
        }

        double sumVi = 0.0;
        for (double v : vi) sumVi += v;
        double s = 0.5 * (totalVar - sumVi);
        return new Pair<double[], Double>(vi, s);
    }

    /**
     * Computes per-class variance + covariance for aggregate fitting.
     */
    private static double[] computeVarianceCovarianceTheoretical(MatrixCell mmap, int m, double t) {
        double[] gt = new double[m];
        for (int i = 0; i < m; i++) {
            MatrixCell mmap2 = new MatrixCell(4);
            mmap2.set(0, mmap.get(0));
            mmap2.set(1, mmap.get(1));
            mmap2.set(2, mmap.get(2 + i));

            Matrix otherClasses = Matrix.zeros(mmap.get(0).getNumRows(), mmap.get(0).getNumCols());
            for (int j = 0; j < m; j++) {
                if (j != i) {
                    otherClasses.addEq(mmap.get(2 + j));
                }
            }
            mmap2.set(3, otherClasses);

            MatrixCell singleClass = new MatrixCell(2);
            singleClass.set(0, mmap2.get(0));
            singleClass.set(1, mmap2.get(2));
            double varI = Map_count_var.map_count_var(singleClass, t);

            double covApprox = 0.0;

            gt[i] = varI + covApprox;
        }
        return gt;
    }
}
