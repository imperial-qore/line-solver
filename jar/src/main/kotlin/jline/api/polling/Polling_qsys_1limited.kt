/**
 * @file 1-Limited Polling System Analysis
 * 
 * Implements analysis algorithms for 1-limited polling systems where the
 * server serves at most one customer per visit to each queue. Provides
 * performance evaluation for multi-queue polling systems with limited service.
 * 
 * @since LINE 3.0
 */
package jline.api.polling

import jline.api.mam.map_lambda
import jline.api.mam.map_mean
import jline.api.mam.map_moment
import jline.api.mam.map_var
import jline.util.matrix.MatrixCell

/**
 * 1-limited polling system analysis algorithms.
 * 
 * Provides methods for analyzing polling systems with 1-limited service discipline.
 * In 1-limited polling, the server serves at most one customer from each queue before
 * moving to the next queue in the polling cycle, regardless of queue length.
 *
 * @since LINE 3.0
 */

/**
 * Computes the exact mean waiting time solution for a polling system with open arrivals.
 * The system assumes that all queues use gated service discipline.
 * The calculation is based on the equations provided by Takagi in ACM Computing Surveys,
 * Vol. 20, No. 1, March 1988, eq (20).
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
 * A[0] = map_exponential(1/0.6);
 * A[1] = map_exponential(1/0.2);
 * S[0] = map_exponential(1.0);
 * S[1] = map_exponential(1.0);
 * C[0] = map_exponential(1.0);
 * C[1] = map_exponential(1.0);
 * double[] W = polling_qsys_1limited(A, S, C);
 * new Matrix(W).print();
` *
</pre> *
 */
fun polling_qsys_1limited(arvMAPs: Array<MatrixCell>,
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
    var d = 0.0

    for (i in 0..<n) {
        lambda[i] = map_lambda(arvMAPs[i])
        b[i] = map_mean(svcMAPs[i])
        b2[i] = map_moment(svcMAPs[i], 2)
        rho1[i] = lambda[i] * b[i]
        rho += rho1[i]
        r1[i] = map_mean(switchMAPs[i])
        r += r1[i]
        delta2[i] = map_var(switchMAPs[i])
        d += delta2[i]
    }

    val W = DoubleArray(n)
    var sumOfSquares = 0.0
    for (i in rho1.indices) {
        sumOfSquares += rho1[i] * rho1[i]
    }
    var sumProduct = 0.0
    for (i in lambda.indices) {
        sumProduct += lambda[i] * b2[i]
    }

    for (i in 0..<n) {
        W[i] = (1 - rho + rho1[i]) / (1 - rho - lambda[i] * r)

        W[i] = W[i] * (1 - rho) / ((1 - rho) * rho + sumOfSquares)
        W[i] =
            W[i] * (rho / (2 * (1 - rho)) * sumProduct + rho * d / 2 / r + r / (2 * (1 - rho)) * rho1[i] * (1 + rho1[i]))
    }
    return W
}
/**
 * Polling Qsys 1Limited algorithms
 */
@Suppress("unused")
class PollingQsys1limited {
    companion object {
        // Class documentation marker for Dokka
    }
}