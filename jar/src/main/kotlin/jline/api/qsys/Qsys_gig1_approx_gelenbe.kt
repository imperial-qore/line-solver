package jline.api.qsys

import java.util.HashMap

/**
 * G/G/1 queue approximation using Gelenbe's method
 * 
 * Simple two-moment approximation for G/G/1 queues
 * 
 * @param lambda arrival rate
 * @param mu service rate
 * @param ca squared coefficient of variation of inter-arrival time
 * @param cs squared coefficient of variation of service time
 * @return HashMap containing W (mean waiting time in system)
 */
@JvmName("qsys_gig1_approx_gelenbe")
fun qsys_gig1_approx_gelenbe(lambda: Double, mu: Double, ca: Double, cs: Double): HashMap<String, Any> {
    val result = HashMap<String, Any>()
    
    val rho = lambda / mu
    val W = (rho * ca + cs) / 2.0 / (1.0 - rho) / lambda
    
    result["W"] = W
    
    return result
}
/**
 * Queueing system gig1 approx gelenbe algorithms
 */
@Suppress("unused")
class QsysGig1ApproxGelenbeAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}