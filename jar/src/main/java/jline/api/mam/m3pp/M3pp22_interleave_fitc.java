/**
 * @file M3PP(2,2) interleaved superposition fitting
 *
 * Implements lumped superposition of multiple M3PP(2,2) processes using interleaved
 * parameter fitting.
 *
 * @since LINE 3.0
 */
package jline.api.mam.m3pp;

import java.util.ArrayList;
import java.util.List;

import jline.util.Pair;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class M3pp22_interleave_fitc {
    private M3pp22_interleave_fitc() {}

    /**
     * Fits L pairs of classes into a single MMAP.
     *
     * @return Pair of (lumped superposition MMAP, list of component M3PPs)
     */
    public static Pair<MatrixCell, List<MatrixCell>> m3pp22_interleave_fitc(
            double[][] av, double[] btv, double[] binfv, double[] stv, double t) {
        int L = av.length;
        if (btv.length != L) throw new IllegalArgumentException("btv must have length L");
        if (binfv.length != L) throw new IllegalArgumentException("binfv must have length L");
        if (stv.length != L) throw new IllegalArgumentException("stv must have length L");
        for (double[] row : av) {
            if (row.length != 2) throw new IllegalArgumentException("av must be L x 2 matrix");
        }

        double[] upperBounds = new double[L];
        double[] lowerBounds = new double[L];
        double[] dValues = new double[L];

        for (int i = 0; i < L; i++) {
            double totalRate = av[i][0] + av[i][1];
            double d = computeD(btv[i], binfv[i], t);
            double z = (binfv[i] - 1.0) * d * d * d * totalRate;
            double u = d * z / (2 * totalRate * totalRate * d * d + z);

            upperBounds[i] = u;
            lowerBounds[i] = d;
            dValues[i] = d;
        }

        double[][] offDiagonalRates = solveOffDiagonalOptimization(upperBounds, lowerBounds, dValues, L);

        List<MatrixCell> m3pps = new ArrayList<MatrixCell>();

        for (int i = 0; i < L; i++) {
            double r1 = 0.0;
            for (int it = 0; it < L; it++) {
                if (it >= i) r1 += offDiagonalRates[0][it];
            }
            double r2 = 0.0;
            for (int it = 0; it < L; it++) {
                if (it <= i) r2 += offDiagonalRates[1][it];
            }

            MatrixCell mmpp = new MatrixCell(2);
            mmpp.set(0, new Matrix(2, 2));
            mmpp.set(1, new Matrix(2, 2));

            mmpp.get(0).set(0, 1, r1);
            mmpp.get(0).set(1, 0, r2);

            double totalRate = av[i][0] + av[i][1];
            double d = r1 + r2;
            double z = (binfv[i] - 1.0) * d * d * d * totalRate;
            double delta = Math.sqrt(z / (2 * r1 * r2));
            double lambda2 = totalRate - r2 / d * delta;
            double lambda1 = lambda2 + delta;

            mmpp.get(1).set(0, 0, lambda1);
            mmpp.get(1).set(1, 1, lambda2);

            mmpp.get(0).set(0, 0, -(mmpp.get(0).get(0, 1) + mmpp.get(1).get(0, 0)));
            mmpp.get(0).set(1, 1, -(mmpp.get(0).get(1, 0) + mmpp.get(1).get(1, 1)));

            double variance = computeMmppVariance(mmpp, t);
            System.out.println("MMPP " + i + " - Var(t): " + variance);

            MatrixCell m3pp = m3pp22_fitc_approx_cov_multiclass(mmpp, av[i], stv[i], t);
            m3pps.add(m3pp);
        }

        MatrixCell lumped = M3pp2m_interleave.m3pp2m_interleave(m3pps);

        return new Pair<MatrixCell, List<MatrixCell>>(lumped, m3pps);
    }

    private static double computeD(double bt1, double binf, double t1) {
        if (!(binf > bt1 && bt1 > 1.0)) {
            throw new IllegalArgumentException(
                    "No solution, infeasible IDC(t): IDC(" + t1 + ") = " + bt1 + ", IDC(inf) = " + binf);
        }
        double c = (binf - 1.0) / (binf - bt1);
        double z = -c * Math.exp(-c);
        double w = lambertW(z, 100, 1e-12);
        return (w + c) / t1;
    }

    private static double lambertW(double z, int maxIter, double tolerance) {
        double w = (z > -0.1) ? z : -1.0;

        for (int iter = 0; iter < maxIter; iter++) {
            double ew = Math.exp(w);
            double wew = w * ew;
            double f = wew - z;
            double df = ew * (w + 1.0);

            if (Math.abs(df) < tolerance) break;

            double delta = f / df;
            w -= delta;

            if (Math.abs(delta) < tolerance) break;
        }
        return w;
    }

    private static double[][] solveOffDiagonalOptimization(double[] upperBounds,
                                                           double[] lowerBounds,
                                                           double[] dValues, int L) {
        double[][] r = new double[2][L];
        for (int i = 0; i < L; i++) {
            double totalD = dValues[i];
            double minU = Math.min(upperBounds[i], totalD * 0.6);
            double minL = Math.min(lowerBounds[i] - minU, totalD * 0.4);

            r[0][i] = minU;
            r[1][i] = Math.max(0.0, minL);

            double sum = 0.0;
            for (int it = 0; it < L; it++) {
                if (it >= i) sum += r[0][it];
                if (it <= i) sum += r[1][it];
            }

            if (Math.abs(sum - totalD) > 1e-6) {
                double adjustment = (totalD - sum) / 2.0;
                r[0][i] += adjustment;
                r[1][i] += adjustment;
            }
        }
        return r;
    }

    private static double computeMmppVariance(MatrixCell mmpp, double t) {
        Matrix D0 = mmpp.get(0);
        Matrix D1 = mmpp.get(1);

        double lambda1 = D1.get(0, 0);
        double lambda2 = D1.get(1, 1);
        double r12 = D0.get(0, 1);
        double r21 = D0.get(1, 0);

        double meanRate = (lambda1 * r21 + lambda2 * r12) / (r12 + r21);
        return meanRate * t * (1.0 + Math.pow(lambda1 - lambda2, 2) / Math.pow(r12 + r21, 2));
    }

    private static MatrixCell m3pp22_fitc_approx_cov_multiclass(MatrixCell mmpp,
                                                                double[] classRates,
                                                                double covariance, double t) {
        MatrixCell m3pp = new MatrixCell(4);
        m3pp.set(0, mmpp.get(0));
        m3pp.set(1, mmpp.get(1));

        double totalRate = classRates[0] + classRates[1];
        if (totalRate > 0) {
            double prop1 = classRates[0] / totalRate;
            double prop2 = classRates[1] / totalRate;

            m3pp.set(2, scaleMatrix(mmpp.get(1), prop1));
            m3pp.set(3, scaleMatrix(mmpp.get(1), prop2));
        } else {
            m3pp.set(2, new Matrix(2, 2));
            m3pp.set(3, new Matrix(2, 2));
        }
        return m3pp;
    }

    private static Matrix scaleMatrix(Matrix m, double factor) {
        Matrix scaled = new Matrix(m.getNumRows(), m.getNumCols());
        for (int i = 0; i < m.getNumRows(); i++) {
            for (int j = 0; j < m.getNumCols(); j++) {
                scaled.set(i, j, m.get(i, j) * factor);
            }
        }
        return scaled;
    }
}
