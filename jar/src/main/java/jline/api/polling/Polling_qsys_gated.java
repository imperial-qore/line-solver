/**
 * @file Gated Polling System Analysis
 *
 * @since LINE 3.0
 */
package jline.api.polling;

import java.util.ArrayList;
import java.util.List;

import jline.api.mam.*;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Polling_qsys_gated {
    private Polling_qsys_gated() {}

    /**
     * Computes mean waiting times for a polling system with gated service discipline.
     */
    public static double[] polling_qsys_gated(MatrixCell[] arvMAPs,
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

        List<double[]> lst1 = new ArrayList<double[]>();
        List<Double> lst2 = new ArrayList<Double>();

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i > j) {
                    double[] t1 = new double[n * n];
                    for (int m = i; m < n; m++) {
                        t1[j * n + m] = -1.0;
                    }
                    for (int m = 0; m <= j - 1; m++) {
                        t1[j * n + m] = -1.0;
                    }
                    for (int m = j; m <= i - 1; m++) {
                        t1[m * n + j] = -1.0;
                    }
                    t1[i * n + j] = 1 / rho1[i];
                    lst1.add(t1);
                    lst2.add(0.0);
                } else if (j > i) {
                    double[] t1 = new double[n * n];
                    for (int m = i; m <= j - 1; m++) {
                        t1[j * n + m] = -1.0;
                    }
                    for (int m = j; m < n; m++) {
                        t1[m * n + j] = -1.0;
                    }
                    for (int m = 0; m <= i - 1; m++) {
                        t1[m * n + j] = -1.0;
                    }
                    t1[i * n + j] = 1 / rho1[i];
                    lst1.add(t1);
                    lst2.add(0.0);
                } else {
                    double[] t1 = new double[n * n];
                    t1[i * n + i] = t1[i * n + i] + 1;
                    for (int m = 0; m < n; m++) {
                        if (i != m) {
                            t1[i * n + m] = -rho1[i];
                        }
                    }
                    for (int m = 0; m < n; m++) {
                        t1[m * n + i] = t1[m * n + i] - (rho1[i] * rho1[i]);
                    }
                    lst1.add(t1);
                    double temp = delta2[i] + lambda[i] * b2[i] * r / (1 - rho);
                    lst2.add(temp);
                }
            }
        }

        int n1 = lst1.size();
        Matrix lhs = new Matrix(n1, n1, n1 * n1);
        for (int i = 0; i < n1; i++) {
            for (int j = 0; j < n1; j++) {
                lhs.set(i, j, lst1.get(i)[j]);
            }
        }
        double[] lst2arr = new double[lst2.size()];
        for (int i = 0; i < lst2.size(); i++) lst2arr[i] = lst2.get(i);
        Matrix rhs = new Matrix(lst2arr);
        Matrix x = Matrix.createLike(rhs);
        Matrix.solve(lhs, rhs, x);
        double[] finalSolution = x.toArray1D();
        double[] W = new double[n];

        for (int i = 0; i < n; i++) {
            double temp = (1 + rho1[i]) * r / (2 * (1 - rho));
            double sum = 0.0;
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    sum += finalSolution[i * n + j];
                }
            }
            sum = sum * (1 / rho1[i]);
            for (int j = 0; j < n; j++) {
                sum += finalSolution[j * n + i];
            }
            temp = temp + (1 - rho) * (1 + rho1[i]) * sum / (2 * r);
            W[i] = temp;
        }
        return W;
    }
}
