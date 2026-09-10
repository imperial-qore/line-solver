/**
 * @file Exhaustive Polling System Analysis
 *
 * @since LINE 3.0
 */
package jline.api.polling;

import java.util.ArrayList;

import org.apache.commons.math3.util.FastMath;

import jline.api.mam.*;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Polling_qsys_exhaustive {
    private Polling_qsys_exhaustive() {}

    /**
     * Computes mean waiting times for an exhaustive polling system.
     */
    public static double[] polling_qsys_exhaustive(MatrixCell[] arvMAPs,
                                                   MatrixCell[] svcMAPs,
                                                   MatrixCell[] switchMAPs) {
        int n = arvMAPs.length;

        double[] lambda = new double[n];
        double[] b = new double[n];
        double[] b2 = new double[n];
        double[] rho1 = new double[n];
        double[] r1 = new double[n];
        double[] delta2 = new double[n];
        double rho = 0.0;
        double r = 0.0;

        for (int i = 0; i < n; i++) {
            lambda[i] = Map_lambda.map_lambda(arvMAPs[i]);
            b[i] = Map_mean.map_mean(svcMAPs[i]);
            b2[i] = Map_moment.map_moment(svcMAPs[i], 2);
            rho1[i] = lambda[i] * b[i];
            rho += rho1[i];
            r1[i] = Map_mean.map_mean(switchMAPs[i]);
            r += r1[i];
            delta2[i] = Map_var.map_var(switchMAPs[i]);
        }

        ArrayList<double[]> lst1 = new ArrayList<double[]>();
        ArrayList<Double> lst2 = new ArrayList<Double>();

        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= n; j++) {
                if (i > j) {
                    double[] t1 = new double[n * n];
                    for (int m = i + 1; m <= n; m++) {
                        t1[(j - 1) * n + m - 1]--;
                    }
                    for (int m = 1; m <= j - 1; m++) {
                        t1[(j - 1) * n + m - 1]--;
                    }
                    for (int m = j; m <= i - 1; m++) {
                        t1[(m - 1) * n + j - 1]--;
                    }
                    t1[(i - 1) * n + j - 1] += (1 - rho1[i - 1]) / rho1[i - 1];
                    lst1.add(t1);
                    lst2.add(0.0);
                } else if (j > i) {
                    double[] t1 = new double[n * n];
                    for (int m = i + 1; m <= j - 1; m++) {
                        t1[(j - 1) * n + m - 1]--;
                    }
                    for (int m = j; m <= n; m++) {
                        t1[(m - 1) * n + j - 1]--;
                    }
                    for (int m = 1; m <= i - 1; m++) {
                        t1[(m - 1) * n + j - 1]--;
                    }
                    t1[(i - 1) * n + j - 1] += (1 - rho1[i - 1]) / rho1[i - 1];
                    lst1.add(t1);
                    lst2.add(0.0);
                } else {
                    double[] t1 = new double[n * n];
                    t1[(i - 1) * n + i - 1]++;
                    for (int m = 1; m <= n; m++) {
                        if (i != m) {
                            t1[(i - 1) * n + m - 1] -= rho1[i - 1] / (1 - rho1[i - 1]);
                        }
                    }
                    lst1.add(t1);
                    double temp = 0.0;
                    if (i > 1) {
                        temp += delta2[i - 2] / FastMath.pow(1 - rho1[i - 1], 2);
                    } else {
                        temp += delta2[n - 1] / FastMath.pow(1 - rho1[i - 1], 2);
                    }
                    temp += lambda[i - 1] * b2[i - 1] * r * (1 - rho1[i - 1])
                            / ((1 - rho) * FastMath.pow(1 - rho1[i - 1], 3));
                    lst2.add(temp);
                }
            }
        }

        int size = lst1.size();
        double[][] A = new double[size][n * n];
        double[] B = new double[size];
        for (int k = 0; k < size; k++) {
            A[k] = lst1.get(k);
            B[k] = lst2.get(k);
        }

        Matrix rhs = new Matrix(B);
        Matrix x = Matrix.createLike(rhs);
        Matrix.solve(new Matrix(A), rhs, x);
        double[] finalSolution = x.toArray1D();
        double[] W = new double[n];

        for (int i = 1; i <= n; i++) {
            double temp = 0.0;
            temp += lambda[i - 1] * b2[i - 1] / (2 * (1 - rho1[i - 1]));
            temp += r * (1 - rho1[i - 1]) / (2 * (1 - rho));
            double sum = 0.0;
            for (int j = 1; j <= n; j++) {
                if (i != j) {
                    sum += finalSolution[(i - 1) * n + j - 1];
                }
            }
            sum *= (1 - rho1[i - 1]) / rho1[i - 1];
            if (i > 1) {
                sum += delta2[i - 2];
            } else {
                sum += delta2[n - 1];
            }
            sum /= r * (1 - rho1[i - 1]) * 2 / (1 - rho);
            temp += sum;
            W[i - 1] = temp;
        }

        return W;
    }
}
