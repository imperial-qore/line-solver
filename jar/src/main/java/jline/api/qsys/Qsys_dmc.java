/**
 * @file D/M/c queueing system analysis
 *
 * Computes time-average performance metrics for a D/M/c queue (deterministic
 * interarrivals, exponential service, c servers) by embedding at arrival
 * epochs and integrating over the inter-arrival cycle.
 */
package jline.api.qsys;

import jline.util.Maths;
import jline.util.matrix.Matrix;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.LUDecomposition;

public final class Qsys_dmc {
    private Qsys_dmc() {}

    public static DmcResult qsys_dmc(double lambdaArr, double mu, int c, int truncation, int quadSteps) {
        if (!(lambdaArr > 0.0)) throw new IllegalArgumentException("Arrival rate must be positive");
        if (!(mu > 0.0)) throw new IllegalArgumentException("Service rate must be positive");
        if (!(c >= 1)) throw new IllegalArgumentException("Number of servers must be >= 1");

        double rho = lambdaArr / (c * mu);
        if (!(rho < 1.0 - 1e-12)) throw new IllegalArgumentException("Load rho=" + rho + " must be strictly less than 1");

        double s = 1.0 / lambdaArr;

        // Cap the truncation so that the dense LU stays under ~10s. ρ=0.99 with
        // n=900 still yields >5 sig figs.
        int autoN = Math.max(200, Math.min(1000, (int) (8.0 / (1.0 - rho)) + 200));
        int nMax = (truncation > 0) ? truncation : autoN;
        int n = nMax + 1;

        // Death-only sub-generator A on [0, n)
        Matrix A = new Matrix(n, n);
        for (int m = 0; m < n; m++) {
            double rate = Math.min(m, c) * mu;
            A.set(m, m, -rate);
            if (m > 0) A.set(m, m - 1, rate);
        }

        Matrix expAs = Maths.matrixExp(A.scale(s));
        Matrix expAdt = Maths.matrixExp(A.scale(s / quadSteps));

        // Build (P^T - I) directly: P[X, e] = expAs[X+1, e] (cap X+1 at n-1)
        double[][] Marr = new double[n][n];
        int[] yIdx = new int[n];
        for (int i = 0; i < n; i++) yIdx[i] = Math.min(i + 1, n - 1);
        for (int X = 0; X < n; X++) {
            int Y = yIdx[X];
            for (int e = 0; e < n; e++) Marr[e][X] = expAs.get(Y, e);
        }
        for (int i = 0; i < n; i++) Marr[i][i] -= 1.0;
        for (int j = 0; j < n; j++) Marr[n - 1][j] = 1.0;
        double[] rhs = new double[n];
        rhs[n - 1] = 1.0;

        double[] piArr = new LUDecomposition(new Array2DRowRealMatrix(Marr, false)).getSolver()
                .solve(new ArrayRealVector(rhs)).toArray();

        // Build q[Y] = sum_{X: yIdx[X]=Y} pi_arr[X] so we can integrate q·expAt·w
        double[] q = new double[n];
        for (int X = 0; X < n; X++) q[yIdx[X]] += piArr[X];

        // Cache expAdt as a flat double[][] for fast vector·matrix products.
        double[][] ed = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                ed[i][j] = expAdt.get(i, j);
            }
        }

        double[] weightsLq = new double[n];
        for (int i = 0; i < n; i++) weightsLq[i] = (double) Math.max(0, i - c);
        double[] weightsN = new double[n];
        for (int i = 0; i < n; i++) weightsN[i] = (double) i;

        double[] v = q.clone();
        double lqSum = 0.5 * dotW(v, weightsLq);
        double nSum = 0.5 * dotW(v, weightsN);
        double[] newV = new double[n];
        for (int k = 1; k <= quadSteps; k++) {
            // newV[j] = sum_i v[i] * ed[i][j]
            for (int j = 0; j < n; j++) newV[j] = 0.0;
            for (int i = 0; i < n; i++) {
                double vi = v[i];
                if (vi == 0.0) continue;
                double[] row = ed[i];
                for (int j = 0; j < n; j++) newV[j] += vi * row[j];
            }
            for (int j = 0; j < n; j++) v[j] = newV[j];
            double wt = (k == quadSteps) ? 0.5 : 1.0;
            lqSum += wt * dotW(v, weightsLq);
            nSum += wt * dotW(v, weightsN);
        }
        double dt = s / quadSteps;
        double Lq = lqSum * dt / s;
        double L = nSum * dt / s;

        double Wq = Lq / lambdaArr;
        double W = Wq + 1.0 / mu;
        return new DmcResult(L, Lq, Wq, W, rho);
    }

    public static DmcResult qsys_dmc(double lambdaArr, double mu, int c, int truncation) {
        return qsys_dmc(lambdaArr, mu, c, truncation, 200);
    }

    public static DmcResult qsys_dmc(double lambdaArr, double mu, int c) {
        return qsys_dmc(lambdaArr, mu, c, -1, 200);
    }

    private static double dotW(double[] v, double[] w) {
        double sum = 0.0;
        for (int i = 0; i < v.length; i++) sum += v[i] * w[i];
        return sum;
    }
}
