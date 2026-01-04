/**
 * @file CTMC stochastic complementarity analysis
 * 
 * Implements stochastic complementarity analysis for CTMCs to identify strongly 
 * connected components and analyze their communication structure. Essential for 
 * decomposing large Markov chains and understanding long-run behavior.
 * 
 * @since LINE 3.0
 */
package jline.api.mc

import jline.solvers.ctmc.SolverCTMC
import jline.util.matrix.Matrix
import kotlin.math.ceil

fun ctmc_stochcomp(Q: Matrix, I_list: MutableList<Double?>): SolverCTMC.StochCompResult {
    
    var I = Matrix(I_list)
    if (I_list.isEmpty()) {
        I = Matrix(1, ceil(Q.getNumCols().toDouble() / 2).toInt())
    }
    val diff_values: MutableList<Double?> = ArrayList<Double?>()
    for (checkValue in 0..<Q.getNumCols()) {
        var check = false
        for (I_row in 0..<I.getNumRows()) {
            if (I.get(I_row) == checkValue.toDouble()) {
                check = true
            }
        }
        if (!check) {
            diff_values.add(checkValue.toDouble())
        }
    }
    val Ic = Matrix(diff_values)
    val Q11: Matrix
    val Q12: Matrix
    val Q21: Matrix?
    val Q22: Matrix

    Q11 = Q.getSubMatrix(I, I)
    Q12 = Q.getSubMatrix(I, Ic)
    Q21 = Q.getSubMatrix(Ic, I)
    Q22 = Q.getSubMatrix(Ic, Ic)
    Q22.mulByMinusOne()

    var T = Q22.inv().mult(Q21)
    T = Q12.mult(T)
    val S = Q11.add(1.0, T)
    return SolverCTMC.StochCompResult(S, Q11, Q12, Q21, Q22, T)
}
/**
 * CTMC stochcomp algorithms
 */
@Suppress("unused")
class CtmcStochcompAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}