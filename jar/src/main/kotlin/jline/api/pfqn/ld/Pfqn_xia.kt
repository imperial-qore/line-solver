/**
 * @file Xia's asymptotic approximation for load-dependent normalizing constants
 * 
 * Implements Xia's asymptotic approximation method for computing normalizing constants
 * in load-dependent closed queueing networks. Provides efficient approximation for
 * high-population regimes with state-dependent service rates and bottleneck identification.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.ld

import jline.solvers.SolverOptions
import jline.util.Maths
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.CombinatoricsUtils
import org.apache.commons.math3.util.FastMath

private fun pfqn_xia_F(u: Double, k: Double): Double {
    var ret = 0.0
    var j = 0
    while (j < k) {
        ret += FastMath.pow(u, j) / CombinatoricsUtils.factorial(j)
        j++
    }
    ret += FastMath.pow(u, k) / CombinatoricsUtils.factorial(k.toInt()) / (1 - u / k)
    return ret
}

fun pfqn_xia(L: Matrix, N: Int, s: Matrix, options: SolverOptions?): Double {
    val L_local = L.copy()
    val rho = Matrix(L_local.numRows, L_local.numCols)
    for (i in 0..<s.numRows) {
        for (j in 0..<s.numCols) {
            rho[i, j] = L_local[i, j] / s[i, j]
        }
    }
    val scalefactor = 1 / rho.elementMax()
    L_local.scaleEq(scalefactor)
    val M = L_local.numRows
    rho.scaleEq(scalefactor)
    val bnkset: MutableList<Int> = ArrayList()
    val nbnkset: MutableList<Int> = ArrayList()
    for (i in 0..<rho.numElements) {
        if (rho[i] == rho.elementMax()) {
            bnkset.add(i)
        }
    }
    for (i in 0..<M) {
        if (!bnkset.contains(i)) {
            nbnkset.add(i)
        }
    }
    val B = bnkset.size
    var logGasy = -Maths.factln(B - 1) - N * FastMath.log(scalefactor)
    for (b in bnkset) {
        logGasy = logGasy + s[b] * FastMath.log(L_local[b]) - Maths.factln(s[b])
    }
    for (k in nbnkset) {
        val f = pfqn_xia_F(L_local[k], s[k])
        if (java.lang.Double.isFinite(f)) {
            logGasy += FastMath.log(f)
        }
    }
    return logGasy
}
/**
 * PFQN xia algorithms
 */
@Suppress("unused")
class PfqnXiaAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}