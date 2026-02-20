package jline.api.mc

import jline.solvers.ctmc.SolverCTMC
import jline.util.matrix.Matrix
import kotlin.math.ceil

fun ctmc_pseudostochcomp(Q: Matrix, I_list: MutableList<Double?>): SolverCTMC.StochCompResult {
    var I = Matrix(I_list)
    if (I_list.isEmpty()) {
        // Default: use first half of states
        val halfSize = ceil(Q.getNumCols().toDouble() / 2).toInt()
        val defaultIndices = mutableListOf<Double?>()
        for (i in 0 until halfSize) {
            defaultIndices.add(i.toDouble())
        }
        I = Matrix(defaultIndices)
    }
    
    // Compute Ic as complement of I
    val diff_values: MutableList<Double?> = ArrayList<Double?>()
    for (checkValue in 0..<Q.getNumCols()) {
        var check = false
        for (I_row in 0..<I.getNumRows()) {
            if (I.get(I_row) == checkValue.toDouble()) {
                check = true
                break
            }
        }
        if (!check) {
            diff_values.add(checkValue.toDouble())
        }
    }
    val Ic = Matrix(diff_values)
    
    // Extract submatrices
    val Q11 = Q.getSubMatrix(I, I)
    val Q12 = Q.getSubMatrix(I, Ic)
    val Q21 = Q.getSubMatrix(Ic, I)
    val Q22 = Q.getSubMatrix(Ic, Ic)
    
    // Solve for stationary distribution
    val abar = ctmc_solve(Q)
    
    // Extract abar(Ic) - stationary probabilities for complement states
    val abar_Ic = abar.getSubMatrix(Ic, Matrix(mutableListOf(0.0)))
    
    // Y = abar(Ic) * Q21
    val Y = abar_Ic.transpose().mult(Q21)
    
    // Create ones vector of size Q12.columns
    val ones = Matrix.ones(Q12.getNumCols(), 1)
    
    // S = Q11 + Q12 * ones * Y / sum(Y)
    val sumY = Y.sumCols().get(0)
    val scaledY = Y.scale(1.0 / sumY)
    val Q12_ones = Q12.mult(ones)
    val addTerm = Q12_ones.mult(scaledY)
    val S = Q11.add(1.0, addTerm)
    
    // For compatibility with existing StochCompResult, we can set T = Q12 * ones * Y / sum(Y)
    val T = addTerm
    
    return SolverCTMC.StochCompResult(S, Q11, Q12, Q21, Q22, T)
}
/**
 * CTMC pseudostochcomp algorithms
 */
@Suppress("unused")
class CtmcPseudostochcompAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}