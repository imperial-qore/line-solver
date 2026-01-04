package jline.api.mam

import jline.util.matrix.Matrix

/**
 * QBD R-matrix computation algorithms.
 * 
 * Provides methods for computing the R-matrix in Quasi-Birth-Death (QBD) processes using
 * matrix analytic methods. The R-matrix is fundamental to analyzing QBD processes and
 * represents the conditional probability that the level increases before decreasing.
 *
 * @since LINE 3.0
 */

/**
 * Compute R matrix using successive substitutions method for QBD processes
 *
 * @param B Backward matrix
 * @param L Local matrix
 * @param F Forward matrix
 * @param iterMax Maximum number of iterations (default: 100000)
 * @return R matrix
 */
fun qbd_R(B: Matrix, L: Matrix, F: Matrix, iterMax: Int = 100000): Matrix {
    val LInv = L.inv()
    val Fil = F.mult(LInv)
    val BiL = B.mult(LInv)
    
    var R = Fil.scale(-1.0)
    var Rprime = Fil.scale(-1.0).sub(R.mult(R).mult(BiL))
    
    for (iter in 1..iterMax) {
        R = Rprime.copy()
        Rprime = Fil.scale(-1.0).sub(R.mult(R).mult(BiL))
        
        val diff = R.sub(Rprime)
        if (diff.norm() <= 1e-12) {
            break
        }
    }
    
    return Rprime
}
/**
 * QBD R algorithms
 */
@Suppress("unused")
class QbdRAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}