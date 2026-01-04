/**
 * @file Quasi-Birth-Death process R-matrix logarithmic reduction
 * 
 * Computes QBD R-matrix using logarithmic reduction method for numerical stability.
 * Advanced algorithm for solving matrix quadratic equations in QBD analysis.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix

/**
 * Compute R matrix using logarithmic reduction method for QBD processes
 *
 * @param B Backward matrix
 * @param L Local matrix
 * @param F Forward matrix
 * @param iterMax Maximum number of iterations (default: 100000)
 * @return R matrix
 */
fun qbd_R_logred(B: Matrix, L: Matrix, F: Matrix, iterMax: Int = 100000): Matrix {
    val r = L.numRows
    val LInv = L.inv()
    
    var iLF = LInv.mult(F).scale(-1.0)
    var iLB = LInv.mult(B).scale(-1.0)
    
    val T = iLF.copy()
    var S = iLB.copy()
    
    for (iter in 1..iterMax) {
        val D = iLF.mult(iLB).add(iLB.mult(iLF))
        val eyeMinusD = Matrix.eye(r).sub(D)
        val eyeMinusDInv = eyeMinusD.inv()
        
        iLF = eyeMinusDInv.mult(iLF.mult(iLF))
        iLB = eyeMinusDInv.mult(iLB.mult(iLB))
        
        S = S.add(T.mult(iLB))
        T.multEq(iLF)
        
        val ones = Matrix.ones(r, 1)
        val convergenceTest = ones.sub(S.mult(ones))
        if (convergenceTest.norm() <= 1e-12) {
            break
        }
    }
    
    val U = L.add(F.mult(S))
    return F.scale(-1.0).mult(U.inv())
}
/**
 * QBD R logred algorithms
 */
@Suppress("unused")
class QbdRLogredAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}