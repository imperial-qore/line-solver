package jline.api.qsys

import java.util.HashMap

/**
 * M/G/1/K loss probability using MacGregor Smith approximation
 * 
 * Calculates the loss probability for an M/G/1/K queue using the
 * MacGregor Smith approximation method.
 * Reference: J. MacGregor Smith - Optimal Design and Performance Modelling 
 * of M/G/1/K Queueing Systems
 * 
 * @param lambda arrival rate
 * @param mu mean service rate (1/mean service time)
 * @param mu_scv squared coefficient of variation of service time
 * @param K system capacity (maximum number of customers)
 * @return HashMap containing lossprob (loss probability) and rho (utilization)
 */
@JvmName("qsys_mg1k_loss_mgs")
fun qsys_mg1k_loss_mgs(lambda: Double, mu: Double, mu_scv: Double, K: Int): HashMap<String, Any> {
    val result = HashMap<String, Any>()
    
    val rho = lambda / mu
    val s = Math.sqrt(mu_scv)
    val sqrt_rho = Math.sqrt(rho)
    
    val exponent1 = (sqrt_rho * s * s - sqrt_rho + 2.0 * K) / (2.0 + sqrt_rho * s * s - sqrt_rho)
    val exponent2 = 2.0 * (1.0 + sqrt_rho * s * s - sqrt_rho + K) / (2.0 + sqrt_rho * s * s - sqrt_rho)
    
    val lossprob_num = Math.pow(rho, exponent1) * (rho - 1.0)
    val lossprob_den = Math.pow(rho, exponent2) - 1.0
    
    val lossprob = lossprob_num / lossprob_den
    
    result["lossprob"] = lossprob
    result["rho"] = rho
    
    return result
}
/**
 * Queueing system mg1k loss mgs algorithms
 */
@Suppress("unused")
class QsysMg1kLossMgsAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}