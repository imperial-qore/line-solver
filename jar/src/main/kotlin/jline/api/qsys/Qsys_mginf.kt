package jline.api.qsys

import jline.GlobalConstants.Inf
import java.util.HashMap
import kotlin.math.exp
import kotlin.math.pow

/**
 * M/G/∞ queue analysis (infinite servers)
 * 
 * Computes exact performance measures for the M/G/∞ queueing system
 * with unlimited servers. This represents a pure delay system where
 * customers never wait in queue.
 * 
 * @param lambda arrival rate
 * @param mu service rate per server
 * @param cv2 squared coefficient of variation of service time (for consistency, not used in calculations)
 * @return HashMap containing:
 *         - L: average number of customers in system
 *         - Lq: average number of customers in queue (always 0)
 *         - W: average time in system (= service time)
 *         - Wq: average waiting time in queue (always 0)
 *         - p0: probability of empty system
 */
@JvmName("qsys_mginf")
fun qsys_mginf(lambda: Double, mu: Double, cv2: Double): HashMap<String, Any> {
    val result = HashMap<String, Any>()
    
    val rho = lambda / mu
    
    // Exact results for M/G/∞
    val L = rho  // Average number of customers = average busy servers
    val Lq = 0.0  // No queue since infinite servers
    val W = 1.0 / mu  // Average time = service time
    val Wq = 0.0  // No waiting time
    val p0 = exp(-rho)  // Poisson distribution
    
    result["L"] = L
    result["Lq"] = Lq
    result["W"] = W
    result["Wq"] = Wq
    result["p0"] = p0
    
    return result
}

/**
 * M/G/∞ queue analysis with state probability
 * 
 * @param lambda arrival rate
 * @param mu service rate per server
 * @param cv2 squared coefficient of variation of service time
 * @param k specific state for probability computation
 * @return HashMap containing all metrics plus pk (probability of k customers)
 */
@JvmName("qsys_mginf")
fun qsys_mginf(lambda: Double, mu: Double, cv2: Double, k: Int): HashMap<String, Any> {
    val result = qsys_mginf(lambda, mu, cv2)
    
    val rho = lambda / mu
    
    // Poisson distribution: P(k) = e^(-ρ) * ρ^k / k!
    val pk = exp(-rho) * rho.pow(k.toDouble()) / factorial(k)
    
    result["pk"] = pk
    
    return result
}

/**
 * Factorial function with overflow protection
 */
private fun factorial(n: Int): Double {
    if (n == 0 || n == 1) return 1.0
    var result = 1.0
    for (i in 2..n) {
        result *= i
        if (result.isInfinite()) {
            return Inf
        }
    }
    return result
}
/**
 * Queueing system mginf algorithms
 */
@Suppress("unused")
class QsysMginfAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}