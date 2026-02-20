package jline.api.qsys

import java.util.HashMap
import kotlin.math.pow
import kotlin.math.sqrt

/**
 * G/G/1 queue approximation using enhanced Myskja's method
 * 
 * Enhanced third moment-based approximation for G/G/1 queues
 * 
 * @param lambda arrival rate
 * @param mu service rate
 * @param ca squared coefficient of variation of inter-arrival time
 * @param cs squared coefficient of variation of service time
 * @param q0 lowest value of the relative third moment for a given mean and SCV
 * @param qa third relative moment E[X^3]/6/E[X]^3, X=inter-arrival time r.v.
 * @return HashMap containing W (mean waiting time in system)
 */
@JvmName("qsys_gig1_approx_myskja2")
fun qsys_gig1_approx_myskja2(lambda: Double, mu: Double, ca: Double, cs: Double, q0: Double, qa: Double): HashMap<String, Any> {
    val result = HashMap<String, Any>()
    
    val ra = (1.0 + ca) / 2.0
    val rs = (1.0 + cs) / 2.0
    val rho = lambda / mu
    
    val theta = (rho * (qa - ra) - (qa - ra * ra)) / (2.0 * rho * (ra - 1.0))
    val d = (1.0 + 1.0 / ra) * (1.0 - rs) * (1.0 - (q0 / qa).pow(3.0)) * (1.0 - rho.pow(3.0))
    val D = (rs - theta).pow(2.0) + (2.0 * rs - 1.0 + d) * (ra - 1.0)
    
    val W = (rho / (1.0 - rho)) / lambda * (rs + (1.0 / rho) * (sqrt(D) - (rs - theta)))
    
    result["W"] = W
    
    return result
}
/**
 * Queueing system gig1 approx myskja2 algorithms
 */
@Suppress("unused")
class QsysGig1ApproxMyskja2Algo {
    companion object {
        // Class documentation marker for Dokka
    }
}