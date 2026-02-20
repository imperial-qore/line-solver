/**
 * @file G/G/1 queueing system analysis
 * 
 * Provides comprehensive analysis of G/G/1 queues with general arrival and service processes.
 * Uses exact methods for special cases (M/M/1, M/G/1, G/M/1) and high-quality approximations
 * like Allen-Cunneen for the general case. Essential for modeling real-world systems.
 *
 * @since LINE 3.0
 */
package jline.api.qsys

import jline.io.Ret.qsys
import java.util.HashMap
import kotlin.math.abs

/**
 * G/G/1 queue analysis (general arrivals and service)
 * 
 * Computes performance measures for the G/G/1 queueing system using
 * exact methods when available (special cases) or high-quality approximations.
 * 
 * @param lambda arrival rate
 * @param mu service rate
 * @param ca2 squared coefficient of variation of inter-arrival time
 * @param cs2 squared coefficient of variation of service time
 * @return HashMap containing:
 *         - L: average number of customers in system
 *         - Lq: average number of customers in queue
 *         - W: average time in system
 *         - Wq: average waiting time in queue
 *         - p0: probability of empty system
 */
@JvmName("qsys_gg1")
fun qsys_gg1(lambda: Double, mu: Double, ca2: Double, cs2: Double): HashMap<String, Any> {
    val result = HashMap<String, Any>()
    
    val rho = lambda / mu
    val tolerance = 1e-8
    
    // Check special cases
    when {
        // M/M/1 case (ca2 = cs2 = 1)
        abs(ca2 - 1.0) < tolerance && abs(cs2 - 1.0) < tolerance -> {
            qsys_mm1(lambda, mu)
            val W = qsys.W
            val Wq = W - 1.0 / mu
            val L = lambda * W
            val Lq = lambda * Wq
            val p0 = 1.0 - rho
            
            result["L"] = L
            result["Lq"] = Lq
            result["W"] = W
            result["Wq"] = Wq
            result["p0"] = p0
            return result
        }
        
        // M/G/1 case (ca2 = 1)
        abs(ca2 - 1.0) < tolerance -> {
            qsys_mg1(lambda, mu, cs2)
            val W = qsys.W
            val Wq = W - 1.0 / mu
            val L = lambda * W
            val Lq = lambda * Wq
            val p0 = 1.0 - rho
            
            result["L"] = L
            result["Lq"] = Lq
            result["W"] = W
            result["Wq"] = Wq
            result["p0"] = p0
            return result
        }
        
        // G/M/1 case (cs2 = 1)
        abs(cs2 - 1.0) < tolerance -> {
            // G/M/1 expects (rho, mu) not (lambda, mu)
            qsys_gm1(rho, mu)
            val W = qsys.W
            val Wq = W - 1.0 / mu
            val L = lambda * W
            val Lq = lambda * Wq
            val p0 = 1.0 - rho
            
            result["L"] = L
            result["Lq"] = Lq
            result["W"] = W
            result["Wq"] = Wq
            result["p0"] = p0
            return result
        }
        
        // General G/G/1 case - use Allen-Cunneen approximation as default
        else -> {
            // Allen-Cunneen is generally considered one of the best approximations
            qsys_gig1_approx_allencunneen(lambda, mu, ca2, cs2)
            
            // Extract W from approximation result
            val W = qsys.W
            val Wq = W - 1.0 / mu
            val L = lambda * W
            val Lq = lambda * Wq
            val p0 = 1.0 - rho
            
            result["L"] = L
            result["Lq"] = Lq
            result["W"] = W
            result["Wq"] = Wq
            result["p0"] = p0
        }
    }
    
    return result
}

/**
 * G/G/1 queue analysis with state probability
 * 
 * @param lambda arrival rate
 * @param mu service rate
 * @param ca2 squared coefficient of variation of inter-arrival time
 * @param cs2 squared coefficient of variation of service time
 * @param k specific state for probability computation
 * @return HashMap containing all metrics plus pk (probability of k customers)
 */
@JvmName("qsys_gg1")
fun qsys_gg1(lambda: Double, mu: Double, ca2: Double, cs2: Double, k: Int): HashMap<String, Any> {
    val result = qsys_gg1(lambda, mu, ca2, cs2)
    
    // For general G/G/1, exact state probabilities are not available
    // Use approximation based on geometric distribution with modified parameter
    val rho = lambda / mu
    val L = result["L"] as Double
    
    // Approximate state probabilities using modified geometric distribution
    val pk = if (k == 0) {
        result["p0"] as Double
    } else {
        // Use a geometric-like approximation
        val p0 = result["p0"] as Double
        val r = (L - rho) / L  // Modification factor
        p0 * Math.pow(rho * (1.0 + r * (ca2 + cs2 - 2.0) / 2.0), k.toDouble())
    }
    
    result["pk"] = pk
    
    return result
}
/**
 * Queueing system gg1 algorithms
 */
@Suppress("unused")
class QsysGg1Algo {
    companion object {
        // Class documentation marker for Dokka
    }
}