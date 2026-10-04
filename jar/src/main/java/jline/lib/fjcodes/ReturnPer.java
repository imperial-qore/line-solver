/**
 * Extract percentiles from Phase-Type distribution
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
package jline.lib.fjcodes;

import java.util.ArrayList;
import java.util.List;

import jline.util.matrix.Matrix;

public final class ReturnPer {
    private ReturnPer() {}

    /**
     * Extract percentiles from a Phase-Type (PH) distribution
     */
    public static Matrix returnPer(Matrix vector, Matrix matrix, double[] pers) {
        int m = matrix.getNumCols();

        Matrix vectorDivMatrix = vector.rightMatrixDivide(matrix);
        double meanRT = -elementSumLocal(vectorDivMatrix);

        double c = 0.0;
        for (int i = 0; i < m; i++) {
            c = Math.max(c, -matrix.get(i, i));
        }

        Matrix P_res = matrix.scale(1.0 / c).add(1.0, Matrix.eye(m));

        Matrix ImMinusPres = Matrix.eye(m).add(-1.0, P_res);
        double M = elementSumLocal(vector.rightMatrixDivide(ImMinusPres));

        double a0 = elementSumLocal(vector);
        double sum_a = a0;
        int k = 0;
        Matrix vP = sumRowsLocal(P_res);
        List<Double> akList = new ArrayList<Double>();

        while (Math.abs(sum_a - M) >= 1e-10) {
            k++;
            double ak = vector.mult(vP).get(0, 0);
            akList.add(Double.valueOf(ak));
            sum_a += ak;
            vP = P_res.mult(vP);
        }
        int K1 = k;
        double[] ak = new double[akList.size()];
        for (int i = 0; i < akList.size(); i++) {
            ak[i] = akList.get(i).doubleValue();
        }

        Matrix percentileRTs = new Matrix(pers.length, 2);

        for (int p = 0; p < pers.length; p++) {
            percentileRTs.set(p, 0, pers[p]);

            if (pers[p] < 1.0 - elementSumLocal(vector)) {
                percentileRTs.set(p, 1, 0.0);
            } else {
                // Solve CDF(t) = pers[p] by bisection. The CDF (a uniformized
                // Poisson sum) is monotone increasing in t, so this converges in
                // ~60 evaluations instead of the previous O(maxTime/0.001) linear
                // scan (which made percentile analysis prohibitively slow).
                double hi = Math.max(3.0 * meanRT, 1e-9);
                int guard = 0;
                while (cdfAt(hi, c, a0, ak, K1) < pers[p] && guard++ < 200) {
                    hi *= 2.0;
                }
                double lo = 0.0;
                for (int it = 0; it < 60; it++) {
                    double mid = 0.5 * (lo + hi);
                    if (cdfAt(mid, c, a0, ak, K1) < pers[p]) {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }
                percentileRTs.set(p, 1, 0.5 * (lo + hi));
            }
        }

        return percentileRTs;
    }

    /**
     * CDF of the uniformized phase-type sojourn time at time t:
     * F(t) = 1 - sum_k a_k * Poisson(k; c*t). Monotone increasing in t.
     */
    private static double cdfAt(double t, double c, double a0, double[] ak, int K1) {
        double pM = Math.exp(-c * t);
        double F = pM * a0;
        for (int ki = 0; ki < K1; ki++) {
            pM = c * t * pM / (ki + 1);
            F += pM * ak[ki];
        }
        return 1.0 - F;
    }

    /**
     * Local element sum (private extension equivalent).
     */
    private static double elementSumLocal(Matrix m) {
        double sum = 0.0;
        for (int i = 0; i < m.getNumRows(); i++) {
            for (int j = 0; j < m.getNumCols(); j++) {
                sum += m.get(i, j);
            }
        }
        return sum;
    }

    /**
     * Local row sums as column vector (private extension equivalent).
     */
    private static Matrix sumRowsLocal(Matrix m) {
        Matrix result = new Matrix(m.getNumRows(), 1);
        for (int i = 0; i < m.getNumRows(); i++) {
            double sum = 0.0;
            for (int j = 0; j < m.getNumCols(); j++) {
                sum += m.get(i, j);
            }
            result.set(i, 0, sum);
        }
        return result;
    }
}
