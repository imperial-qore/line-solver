package jline.api.qsys

import jline.api.mc.dtmc_makestochastic
import jline.api.mc.dtmc_solve
import jline.util.matrix.Matrix
import java.util.HashMap
import kotlin.math.exp
import kotlin.math.pow

/**
 * M/G/1/K loss probability calculation
 * 
 * Calculates the loss probability for an M/G/1/K queue using the
 * Niu-Cooper embedded chain analysis method.
 * Reference: Niu-Cooper, Transform-Free Analysis of M/G/1/K and Related Queues,
 * Mathematics of Operations Research Vol. 18, No. 2 (May, 1993), pp. 486-510
 * 
 * @param lambda arrival rate
 * @param svc_density service time density function (t -> density at t)
 * @param K system capacity (maximum number of customers)
 * @return HashMap containing lossprob (loss probability) and rho (utilization)
 */
@JvmName("qsys_mg1k_loss")
fun qsys_mg1k_loss(lambda: Double, svc_density: (Double) -> Double, K: Int): HashMap<String, Any> {
    val result = HashMap<String, Any>()
    
    val tmax = 10.0 / lambda  // Integration upper limit
    
    // Calculate mean service time
    val mu = 1.0 / integrate({ t -> t * svc_density(t) }, 0.0, tmax)
    
    // Calculate arrival probabilities a_j
    val a = DoubleArray(K)
    for (j in 0 until K) {
        val factj = factorial(j)
        a[j] = integrate({ t -> 
            exp(-lambda * t) * (lambda * t).pow(j.toDouble()) * svc_density(t) 
        }, 0.0, tmax) / factj
        
        if (a[j] < 1e-12) {
            break
        }
    }
    
    // Build transition matrix P
    val P = Matrix(K - 1, K - 1)
    
    // Traditional embedded process at departure
    for (i in 0 until K - 1) {
        P[0, i] = a[i]
        if (K - 1 > 1) {
            P[1, i] = a[i]
        }
    }
    P[0, K - 2] = 1.0 - P.getRow(0).elementSum() + P[0, K - 2]
    if (K - 1 > 1) {
        P[1, K - 2] = 1.0 - P.getRow(1).elementSum() + P[1, K - 2]
    }
    
    for (j in 2 until K - 1) {
        for (i in j - 1 until K - 1) {
            P[j, i] = a[i - j + 1]
        }
        P[j, K - 2] = 1.0 - P.getRow(j).elementSum() + P[j, K - 2]
    }
    
    // Make stochastic and solve
    val Pstoch = dtmc_makestochastic(P)
    val sigma = dtmc_solve(Pstoch)
    
    val rho = lambda / mu
    val lossprob = 1.0 - 1.0 / (sigma[0] * a[0] + rho)
    
    result["lossprob"] = lossprob
    result["rho"] = rho
    
    return result
}

/**
 * Simple numerical integration using trapezoidal rule
 */
private fun integrate(f: (Double) -> Double, a: Double, b: Double, n: Int = 1000): Double {
    val h = (b - a) / n
    var sum = 0.5 * (f(a) + f(b))
    for (i in 1 until n) {
        sum += f(a + i * h)
    }
    return sum * h
}

/**
 * Factorial function with memoization for efficiency
 */
private val factorialCache = mutableMapOf<Int, Double>()
private fun factorial(n: Int): Double {
    if (n == 0 || n == 1) return 1.0
    return factorialCache.getOrPut(n) { n * factorial(n - 1) }
}
/**
 * Queueing system mg1k loss algorithms
 */
@Suppress("unused")
class QsysMg1kLossAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}