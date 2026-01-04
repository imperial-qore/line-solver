/**
 * @file Multi-class repairman model integration using Simpson's rule
 * 
 * Implements numerical integration for computing normalizing constants in multi-class
 * repairman models using Simpson's rule integration. Specialized for single-station
 * networks with multiple job classes and think times.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.nc

import jline.io.Ret
import jline.util.matrix.Matrix
import org.apache.commons.math3.analysis.UnivariateFunction
import org.apache.commons.math3.analysis.integration.BaseAbstractUnivariateIntegrator
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator
import org.apache.commons.math3.util.FastMath

fun pfqn_mmint2(L: Matrix, N: Matrix, Z: Matrix): Ret.pfqnNc {
    val nnzClasses = N.findNonNegative().toIntArray1D()
    val order = 12

    val f = { u: Double ->
        val expTerm = FastMath.exp(-u)
        var prodTerm = 1.0
        for (j in nnzClasses) {
            val term = Z.get(j) + L.get(j) * u
            prodTerm *= FastMath.pow(term, N.get(j))
        }
        expTerm * prodTerm
    }

    val p = 1 - FastMath.pow(10.0, -order.toDouble())
    val exp1prctile = -FastMath.log(1 - p)

    val func = UnivariateFunction { t -> f(t) }
    val integrator: BaseAbstractUnivariateIntegrator = SimpsonIntegrator(1e-12, 1e-8, 3, 64)

    val integralValue: Double = try {
        integrator.integrate(Int.MAX_VALUE, func, 0.0, exp1prctile)
    } catch (e: Exception) {
        throw RuntimeException("Integration failed", e)
    }

    val Nmat = Matrix(N)
    val lG = FastMath.log(integralValue) - Nmat.factln().elementSum()
    val G = FastMath.exp(lG)

    return Ret.pfqnNc(G, lG)
}
/**
 * PFQN mmint2 algorithms
 */
@Suppress("unused")
class PfqnMmint2Algo {
    companion object {
        // Class documentation marker for Dokka
    }
}