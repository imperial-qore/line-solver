package jline.api.qsys

import jline.io.Ret.qsys
import java.util.HashMap
import kotlin.math.pow
import kotlin.math.sqrt

/**
 * G/G/k queue approximation using Whitt's method
 * 
 * This approximation uses the Queue Network Analyzer (QNA) methodology
 * to provide accurate estimates for multi-server queues with general distributions.
 * 
 * @param lambda arrival rate
 * @param mu service rate per server
 * @param ca2 squared coefficient of variation of inter-arrival time
 * @param cs2 squared coefficient of variation of service time
 * @param k number of servers
 * @return HashMap containing:
 *         - L: average number of customers in system
 *         - W: average time in system
 *         - Q: average queue length
 *         - U: server utilization
 */
@JvmName("qsys_gigk_approx_whitt")
fun qsys_gigk_approx_whitt(lambda: Double, mu: Double, ca2: Double, cs2: Double, k: Int): HashMap<String, Any> {
    val result = HashMap<String, Any>()
    
    val rho = lambda / (k * mu)
    
    // First get M/M/k results as baseline
    qsys_mmk(lambda, mu, k)
    val W_mmk = qsys.W
    val Wq_mmk = W_mmk - 1.0 / mu
    
    // Whitt's approximation factors
    
    // Heavy traffic factor
    val htFactor = if (rho > 0.7) {
        // Enhanced accuracy in heavy traffic
        1.0 + (4.0 * (rho - 0.7)).pow(2.0)
    } else {
        1.0
    }
    
    // Variability correction factor
    val g = if (ca2 >= 1.0 && cs2 >= 1.0) {
        1.0
    } else if (ca2 <= 1.0 && cs2 <= 1.0) {
        // Low variability correction
        val phi = (1.0 - ca2) * (1.0 - cs2) / (1.0 + cs2)
        1.0 - phi * sqrt(k.toDouble())
    } else {
        // Mixed variability
        (ca2 + cs2) / 2.0
    }
    
    // QNA formula
    val Wq = Wq_mmk * ((ca2 + cs2) / 2.0) * g * htFactor
    
    // Ensure non-negative waiting time
    val Wq_final = if (Wq > 0) Wq else 0.0
    
    val W = Wq_final + 1.0 / mu
    val L = lambda * W
    val Q = lambda * Wq_final
    val U = rho
    
    result["L"] = L
    result["W"] = W
    result["Q"] = Q
    result["U"] = U
    
    return result
}
/**
 * Queueing system gigk approx whitt algorithms
 */
@Suppress("unused")
class QsysGigkApproxWhittAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}