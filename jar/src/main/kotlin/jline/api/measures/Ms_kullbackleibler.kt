/**
 * @file Kullback-Leibler divergence metric
 * 
 * Implements the Kullback-Leibler divergence D_KL(P||Q) = ∑pᵢlog(pᵢ/qᵢ) measuring 
 * the relative entropy between two probability distributions. A fundamental measure 
 * in information theory for distribution comparison and model selection.
 * 
 * @since LINE 3.0
 */
package jline.api.measures

import jline.GlobalConstants.Inf

import jline.util.matrix.Matrix
import kotlin.math.ln

/**
 * Kullback-Leibler divergence between two probability distributions.
 * Part of Shannon's entropy family.
 * Measures information lost when Q is used to approximate P.
 * 
 * @param P exact probability distribution
 * @param Q model probability distribution
 * @return KL divergence
 */
fun ms_kullbackleibler(P: Matrix, Q: Matrix): Double {
    require(P.length() == Q.length()) { "Distributions must have the same size" }
    
    var kl = 0.0
    for (i in 0 until P.length()) {
        val pi = P.get(i)
        val qi = Q.get(i)
        if (pi > 0) {
            if (qi == 0.0) {
                return Inf
            }
            kl += pi * ln(pi / qi)
        }
    }
    
    return kl
}
/**
 * Kullbackleibler metric algorithms
 */
@Suppress("unused")
class MsKullbackleiblerAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}