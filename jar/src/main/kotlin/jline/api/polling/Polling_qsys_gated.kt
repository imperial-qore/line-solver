/**
 * @file Gated Polling System Analysis
 * 
 * Implements analysis algorithms for gated polling systems where the server
 * serves all customers present at the beginning of a visit period. Provides
 * mean-value analysis for multi-queue polling systems with gated service.
 * 
 * @since LINE 3.0
 */
package jline.api.polling

import jline.api.mam.map_lambda
import jline.api.mam.map_mean
import jline.api.mam.map_moment
import jline.api.mam.map_var
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Gated polling system analysis algorithms.
 * 
 * Provides methods for analyzing polling systems with gated service discipline.
 * In gated polling, each queue is served until all customers present at the
 * beginning of the service period are served, then the server moves to the next queue.
 *
 * @since LINE 3.0
 */

/**
 * Computes the exact mean waiting time solution for a polling system with open arrivals,
 * where all queues use gated service discipline. The calculation is based on the equations
 * provided by Takagi in ACM Computing Surveys, Vol. 20, No. 1, March 1988, eq (20).
 *
 * @param arvMAPs    an array of MatrixCell objects representing the arrival process MAPs.
 * @param svcMAPs    an array of MatrixCell objects representing the service process MAPs.
 * @param switchMAPs an array of MatrixCell objects representing the switching times MAPs.
 * @return a double array containing the mean waiting times for each queue in the system.
 *
 *
 * Example usage:
 * <pre>
 * `MatrixCell[] A = new MatrixCell[2];
 * MatrixCell[] S = new MatrixCell[2];
 * MatrixCell[] C = new MatrixCell[2];
 * A[0] = map_exponential(1 / 0.6);
 * A[1] = map_exponential(1 / 0.2);
 * S[0] = map_exponential(1.0);
 * S[1] = map_exponential(1.0);
 * C[0] = map_exponential(1.0);
 * C[1] = map_exponential(1.0);
 * double[] W = polling_qsys_gated(A, S, C);
 * new Matrix(W).print();
` *
</pre> *
 */
fun polling_qsys_gated(arvMAPs: Array<MatrixCell>,
                       svcMAPs: Array<MatrixCell>,
                       switchMAPs: Array<MatrixCell>): DoubleArray {
    run {
        val n = arvMAPs.size
        val lambda = DoubleArray(n)
        val b = DoubleArray(n)
        val b2 = DoubleArray(n)
        val rho1 = DoubleArray(n)
        val r1 = DoubleArray(n)
        val delta2 = DoubleArray(n)
        var rho = 0.0
        var r = 0.0

        for (i in 0..<n) {
            lambda[i] = map_lambda(arvMAPs[i])
            b[i] = map_mean(svcMAPs[i])
            b2[i] = map_moment(svcMAPs[i], 2)
            rho1[i] = lambda[i] * b[i]
            rho += rho1[i]
            r1[i] = map_mean(switchMAPs[i])
            r += r1[i]
            delta2[i] = map_var(switchMAPs[i])
        }

        val lst1: MutableList<DoubleArray> = ArrayList()
        val lst2: MutableList<Double> = ArrayList()

        for (i in 0..<n) {
            for (j in 0..<n) {
                if (i > j) {
                    val t1 = DoubleArray(n * n)
                    for (m in i..<n) {
                        t1[j * n + m] = -1.0
                    }
                    for (m in 0..j - 1) {
                        t1[j * n + m] = -1.0
                    }
                    for (m in j..i - 1) {
                        t1[m * n + j] = -1.0
                    }
                    t1[i * n + j] = 1 / rho1[i]
                    lst1.add(t1)
                    lst2.add(0.0)
                } else if (j > i) {
                    val t1 = DoubleArray(n * n)
                    for (m in i..j - 1) {
                        t1[j * n + m] = -1.0
                    }
                    for (m in j..<n) {
                        t1[m * n + j] = -1.0
                    }
                    for (m in 0..i - 1) {
                        t1[m * n + j] = -1.0
                    }
                    t1[i * n + j] = 1 / rho1[i]
                    lst1.add(t1)
                    lst2.add(0.0)
                } else {
                    val t1 = DoubleArray(n * n)
                    t1[i * n + i] = t1[i * n + i] + 1
                    for (m in 0..<n) {
                        if (i != m) {
                            t1[i * n + m] = -rho1[i]
                        }
                    }
                    for (m in 0..<n) {
                        t1[m * n + i] = t1[m * n + i] - (rho1[i] * rho1[i])
                    }
                    lst1.add(t1)
                    val temp = delta2[i] + lambda[i] * b2[i] * r / (1 - rho)
                    lst2.add(temp)
                }
            }
        }

        val n1 = lst1.size
        val lhs = Matrix(n1, n1, n1 * n1)
        for (i in 0..<n1) {
            for (j in 0..<n1) {
                lhs[i, j] = lst1[i][j]
            }
        }
        val rhs = Matrix(lst2)
        val x = Matrix.createLike(rhs)
        Matrix.solve(lhs, rhs, x)
        val finalSolution = x.toArray1D()
        val W = DoubleArray(n)

        for (i in 0..<n) {
            var temp = (1 + rho1[i]) * r / (2 * (1 - rho))
            var sum = 0.0
            for (j in 0..<n) {
                if (i != j) {
                    sum += finalSolution[i * n + j]
                }
            }
            sum = sum * (1 / rho1[i])
            for (j in 0..<n) {
                sum += finalSolution[j * n + i]
            }
            temp = temp + (1 - rho) * (1 + rho1[i]) * sum / (2 * r)
            W[i] = temp
        }
        return W
    }
}
/**
 * Polling Qsys Gated algorithms
 */
@Suppress("unused")
class PollingQsysGated {
    companion object {
        // Class documentation marker for Dokka
    }
}