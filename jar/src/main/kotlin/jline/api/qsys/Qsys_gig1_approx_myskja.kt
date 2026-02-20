package jline.api.qsys

import java.util.HashMap
import kotlin.math.pow

/**
 * G/G/1 queue approximation using Myskja's method
 * 
 * Third moment-based approximation for G/G/1 queues
 * 
 * @param lambda arrival rate
 * @param mu service rate
 * @param ca squared coefficient of variation of inter-arrival time
 * @param cs squared coefficient of variation of service time
 * @param q0 lowest value of the relative third moment for a given mean and SCV
 * @param qa third relative moment E[X^3]/6/E[X]^3, X=inter-arrival time r.v.
 * @return HashMap containing W (mean waiting time in system)
 */
@JvmName("qsys_gig1_approx_myskja")
fun qsys_gig1_approx_myskja(lambda: Double, mu: Double, ca: Double, cs: Double, q0: Double, qa: Double): HashMap<String, Any> {
    val result = HashMap<String, Any>()
    
    val rho = lambda / mu
    val q = qa  // Assuming q is the same as qa based on the MATLAB code
    
    val W = rho / (2.0 * mu * (1.0 - rho)) * 
            ((1.0 + cs) + (q0 / q).pow(1.0 / rho - rho) * (1.0 / rho) * (ca - 1.0))
    
    result["W"] = W
    
    return result
}
/**
 * Queueing system gig1 approx myskja algorithms
 */
@Suppress("unused")
class QsysGig1ApproxMyskjaAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}