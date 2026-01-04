package jline.api.qsys

import java.util.HashMap

/**
 * M/M/1/K loss probability calculation
 * 
 * Calculates the loss probability for an M/M/1/K queue (finite capacity K)
 * Based on: Niu-Cooper, Transform-Free Analysis of M/G/1/K and Related Queues,
 * Mathematics of Operations Research Vol. 18, No. 2 (May, 1993), pp. 486-510
 * 
 * @param lambda arrival rate
 * @param mu service rate
 * @param K system capacity (maximum number of customers)
 * @return HashMap containing lossprob (loss probability) and rho (utilization)
 */
@JvmName("qsys_mm1k_loss")
fun qsys_mm1k_loss(lambda: Double, mu: Double, K: Int): HashMap<String, Any> {
    val result = HashMap<String, Any>()
    
    val rho = lambda / mu
    
    val lossprob = if (Math.abs(rho - 1.0) < 1e-10) {
        // Special case when rho = 1
        1.0 / (K + 1.0)
    } else {
        // General case
        (1.0 - rho) / (1.0 - Math.pow(rho, K + 1.0)) * Math.pow(rho, K.toDouble())
    }
    
    result["lossprob"] = lossprob
    result["rho"] = rho
    
    return result
}
/**
 * Queueing system mm1k loss algorithms
 */
@Suppress("unused")
class QsysMm1kLossAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}