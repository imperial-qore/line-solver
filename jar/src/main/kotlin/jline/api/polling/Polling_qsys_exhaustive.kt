/**
 * @file Exhaustive Polling System Analysis
 * 
 * Implements analysis algorithms for exhaustive polling systems where the
 * server continues serving until a queue becomes empty. Provides performance
 * evaluation for multi-queue polling systems with exhaustive service discipline.
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
import org.apache.commons.math3.util.FastMath

/**
 * Exhaustive polling system analysis algorithms.
 * 
 * Provides methods for analyzing polling systems with exhaustive service discipline.
 * In exhaustive polling, the server continues to serve a queue until it becomes empty
 * before moving to the next queue in the polling cycle.
 *
 * @since LINE 3.0
 */

/**
 * Computes the exact mean waiting time solution for a polling system with open arrivals,
 * where all queues are served exhaustively. The calculation is based on the equations
 * provided by Takagi in ACM Computing Surveys, Vol. 20, No. 1, March 1988, eq (15).
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
 * S[0] = map_exponential(0.1);
 * S[1] = map_exponential(1.0);
 * C[0] = map_exponential(1.0);
 * C[1] = map_hyperexp(1, 2, 0.9);
 * double[] W = polling_qsys_exhaustive(A, S, C);
 * new Matrix(W).print();
` *
</pre> *
 */
fun polling_qsys_exhaustive(arvMAPs: Array<MatrixCell>,
                            svcMAPs: Array<MatrixCell>,
                            switchMAPs: Array<MatrixCell>): DoubleArray {
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

    val lst1 = ArrayList<DoubleArray>()
    val lst2 = ArrayList<Double>()

    for (i in 1..n) {
        for (j in 1..n) {
            if (i > j) {
                val t1 = DoubleArray(n * n)
                for (m in i + 1..n) {
                    t1[(j - 1) * n + m - 1]--
                }
                for (m in 1..j - 1) {
                    t1[(j - 1) * n + m - 1]--
                }
                for (m in j..i - 1) {
                    t1[(m - 1) * n + j - 1]--
                }
                t1[(i - 1) * n + j - 1] += (1 - rho1[i - 1]) / rho1[i - 1]
                lst1.add(t1)
                lst2.add(0.0)
            } else if (j > i) {
                val t1 = DoubleArray(n * n)
                for (m in i + 1..j - 1) {
                    t1[(j - 1) * n + m - 1]--
                }
                for (m in j..n) {
                    t1[(m - 1) * n + j - 1]--
                }
                for (m in 1..i - 1) {
                    t1[(m - 1) * n + j - 1]--
                }
                t1[(i - 1) * n + j - 1] += (1 - rho1[i - 1]) / rho1[i - 1]
                lst1.add(t1)
                lst2.add(0.0)
            } else {
                val t1 = DoubleArray(n * n)
                t1[(i - 1) * n + i - 1]++
                for (m in 1..n) {
                    if (i != m) {
                        t1[(i - 1) * n + m - 1] -= rho1[i - 1] / (1 - rho1[i - 1])
                    }
                }
                lst1.add(t1)
                var temp = 0.0
                temp += if (i > 1) {
                    delta2[i - 2] / FastMath.pow(1 - rho1[i - 1], 2)
                } else {
                    delta2[n - 1] / FastMath.pow(1 - rho1[i - 1], 2)
                }
                temp += lambda[i - 1] * b2[i - 1] * r * (1 - rho1[i - 1]) / ((1 - rho) * FastMath.pow(1 - rho1[i - 1],
                    3))
                lst2.add(temp)
            }
        }
    }

    val size = lst1.size
    val A = Array(size) { DoubleArray(n * n) }
    val B = DoubleArray(size)
    for (k in 0..<size) {
        A[k] = lst1[k]
        B[k] = lst2[k]
    }

    val rhs = Matrix(B)
    val x = Matrix.createLike(rhs)
    Matrix.solve(Matrix(A), rhs, x)
    val finalSolution = x.toArray1D()
    val W = DoubleArray(n)

    for (i in 1..n) {
        var temp = 0.0
        temp += lambda[i - 1] * b2[i - 1] / (2 * (1 - rho1[i - 1]))
        temp += r * (1 - rho1[i - 1]) / (2 * (1 - rho))
        var sum = 0.0
        for (j in 1..n) {
            if (i != j) {
                sum += finalSolution[(i - 1) * n + j - 1]
            }
        }
        sum *= (1 - rho1[i - 1]) / rho1[i - 1]
        sum += if (i > 1) {
            delta2[i - 2]
        } else {
            delta2[n - 1]
        }
        sum /= r * (1 - rho1[i - 1]) * 2 / (1 - rho)
        temp += sum
        W[i - 1] = temp
    }

    // % station time method Ferguson and Aminetzah 1985 as reported by Takagi
    return W
}
/**
 * Polling Qsys Exhaustive algorithms
 */
@Suppress("unused")
class PollingQsysExhaustive {
    companion object {
        // Class documentation marker for Dokka
    }
}