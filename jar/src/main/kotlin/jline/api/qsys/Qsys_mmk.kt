/**
 * @file M/M/k queueing system analysis
 * 
 * Implements exact analytical solutions for M/M/k queues with Poisson arrivals,
 * exponential service times, and k parallel servers. Uses the Erlang-C formula
 * to compute blocking probabilities and performance metrics for multi-server systems.
 *
 * @since LINE 3.0
 */
package jline.api.qsys

import jline.io.Ret.qsys
import org.apache.commons.math3.util.FastMath

/**
 * Analyzes an M/M/k queueing system.
 *
 * @param lambda Arrival rate.
 * @param mu     Service rate.
 * @param k      Number of servers.
 * @return qsysReturn containing average waiting time (W) and utilization (rho).
 */
fun qsys_mmk(lambda: Double, mu: Double, k: Int): qsys {
    val rho = lambda / mu / k
    val Q = rho / (1 - rho) * ErlangC(rho, k) + k * rho
    val W = Q / lambda
    return qsys(W, rho)
}


/**
 * Calculates the probability that an arriving customer is forced to join the queue (i.e., all servers are occupied) in an M/M/k system.
 *
 * @param nu Utilization.
 * @param C  The number of servers.
 * @return Probability that an arriving customer is forced to join the queue.
 */
fun ErlangC(nu: Double, C: Int): Double {
    var S = 0.0
    var factj: Int
    var factj_1 = 1
    for (j in 0..C - 1) {
        if (j == 0) {
            S += FastMath.pow(C * nu, j)
        } else {
            factj = j * factj_1
            S += FastMath.pow(C * nu, j) / factj
            factj_1 = factj
        }
    }
    return 1.0 / (1 + (1 - nu) * (C * factj_1) / FastMath.pow(C * nu, C) * S)
}
/**
 * Queueing system mmk algorithms
 */
@Suppress("unused")
class QsysMmkAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}