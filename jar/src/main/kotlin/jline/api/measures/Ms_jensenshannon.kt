/**
 * @file Jensen-Shannon divergence metric
 * 
 * Implements Jensen-Shannon divergence, a symmetric and bounded version of 
 * Kullback-Leibler divergence. Measures the similarity between probability 
 * distributions with values between 0 and 1, commonly used in phylogenetics and text analysis.
 * 
 * @since LINE 3.0
 */
package jline.api.measures

import jline.util.matrix.Matrix
import kotlin.math.ln

/**
 * Jensen-Shannon divergence between two probability distributions.
 * Part of Shannon's entropy family. Symmetric version of KL divergence.
 * 
 * @param P first probability distribution
 * @param Q second probability distribution
 * @return Jensen-Shannon divergence
 */
fun ms_jensenshannon(P: Matrix, Q: Matrix): Double {
    require(P.length() == Q.length()) { "Distributions must have the same size" }
    
    var sum1 = 0.0
    var sum2 = 0.0
    
    for (i in 0 until P.length()) {
        val pi = P.get(i)
        val qi = Q.get(i)
        val m = (pi + qi) / 2.0
        
        if (pi > 0 && m > 0) {
            sum1 += pi * ln(pi / m)
        }
        
        if (qi > 0 && m > 0) {
            sum2 += qi * ln(qi / m)
        }
    }
    
    return 0.5 * (sum1 + sum2)
}
/**
 * Jensenshannon metric algorithms
 */
@Suppress("unused")
class MsJensenshannonAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}