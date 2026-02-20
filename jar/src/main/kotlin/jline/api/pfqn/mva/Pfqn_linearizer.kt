/**
 * Linearizer Approximate MVA for Product-Form Networks
 * 
 * Implements the linearizer approximate MVA method for large closed queueing networks
 * where exact MVA becomes computationally prohibitive. Provides near-exact accuracy
 * with significantly reduced computational complexity for multi-class systems.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.mva

import jline.io.Ret
import jline.lang.constant.SchedStrategy
import jline.util.matrix.Matrix

/**
 * Linearizer approximate mean value analysis algorithm
 */
fun pfqn_linearizer(L: Matrix,
                    N: Matrix,
                    Z: Matrix,
                    type: Array<SchedStrategy>,
                    tol: Double = 1.0e-8,
                    maxiter: Int = 1000): Ret.pfqnAMVA {
    return pfqn_gflinearizer(L, N, Z, type, tol, maxiter, 1.0)
}
/**
 * PFQN linearizer algorithms
 */
@Suppress("unused")
class PfqnLinearizerAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}