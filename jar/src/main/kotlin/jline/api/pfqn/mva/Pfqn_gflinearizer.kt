/**
 * @file General-form linearizer approximate MVA
 * 
 * Implements the general-form linearizer approximation for closed queueing networks with
 * configurable linearization parameters. Provides wrapper for the extended general-form
 * linearizer algorithm with scalar linearization factor.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.mva

import jline.io.Ret
import jline.lang.constant.SchedStrategy
import jline.util.matrix.Matrix


/* General-form linearizer algorithm */
fun pfqn_gflinearizer(L: Matrix,
                      N: Matrix,
                      Z: Matrix,
                      type: Array<SchedStrategy>,
                      tol: Double,
                      maxiter: Int,
                      alpha: Double): Ret.pfqnAMVA {
    val alphaM = Matrix(1, N.numCols)
    alphaM.fill(alpha)
    return pfqn_egflinearizer(L, N, Z, type, tol, maxiter, alphaM)
}
/**
 * PFQN gflinearizer algorithms
 */
@Suppress("unused")
class PfqnGflinearizerAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}