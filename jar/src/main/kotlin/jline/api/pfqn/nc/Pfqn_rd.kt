/**
 * @file Random Discretization (RD) method for load-dependent normalizing constants
 * 
 * Implements the Random Discretization approach for computing normalizing constants in
 * load-dependent closed queueing networks. Approximates continuous load-dependent service
 * rates through discrete randomization techniques with gamma distribution sampling.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.nc

import jline.GlobalConstants.Inf
import jline.GlobalConstants.NegInf
import jline.api.pfqn.ld.pfqn_gldsingle
import jline.api.pfqn.mva.pfqn_mva
import jline.io.Ret
import jline.solvers.SolverOptions
import jline.solvers.nc.SolverNC
import jline.util.Maths
import jline.util.Utils
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath
import java.util.*
import kotlin.math.abs

fun pfqn_rd(L: Matrix, N: Matrix, Z: Matrix, mu: Matrix, options: SolverOptions?): Ret.pfqnRd {
    var options = options
    if (options == null) {
        options = SolverNC.defaultOptions()
    }
    val method = options.method
    val M = L.numRows
    val R = L.numCols
    val lambda = Matrix(1, R)
    lambda.zero()
    if (N.elementSum() < 0) {
        return Ret.pfqnRd(NegInf)
    }
    var isLi = true
    for (i in 0..<M) {
        for (j in 1..<R) {
            if (mu[i, j] != mu[i, 0]) {
                isLi = false
                break
            }
        }
        if (isLi) {
            for (j in 1..<R) {
                L[i, j] = L[i, j] / mu[i, 0]
            }
        }
    }
    if (N.elementSum() == 0.0) {
        return Ret.pfqnRd(0.0)
    }

    val gamma = Matrix(M, FastMath.ceil(N.elementSum()).toInt())
    gamma.ones()

    var mu_new: Matrix
    if (N.elementSum().toInt() >= mu.numCols) {
        mu_new = mu.copy()
    } else {
        mu_new = Matrix(mu.numRows, 0)
        var i = 0
        while (i < N.elementSum()) {
            val mu_col_i = Matrix(mu.numRows, 1)
            Matrix.extract(mu, 0, mu.numRows, i, i + 1, mu_col_i, 0, 0)
            mu_new = Matrix.concatColumns(mu_new, mu_col_i, null)
            i++
        }
    }

    for (i in 0..<mu_new.numRows) {
        for (j in 0..<mu_new.numCols) {
            if (java.lang.Double.isNaN(mu_new[i, j])) {
                mu_new[i, j] = Inf
            }
        }
    }
    val s = Matrix(M, 1)
    for (i in 0..<M) {
        if (mu[i, mu.numCols - 1] != Inf) {
            val mu_i = Matrix(1, mu.numCols)
            for (j in 0..<mu.numCols) {
                val mu_ij = FastMath.abs(mu[i, j] - mu[i, mu.numCols - 1])
                mu_i[j] = (if (mu_ij < options.tol) 1 else 0).toDouble()
            }
            s[i] = mu_i.find().elementMin()
            //					if (s.get(i) == -1) {
//						s.set(i, N.elementSum() - 1);
//					}
        } else {
            s[i] = N.elementSum() - 1
        }
    }
    val isDelay = Matrix(ArrayList(Collections.nCopies(M, 0.0)))
    Matrix(ArrayList(Collections.nCopies(M, 0.0)))
    val y = L.copy()

    for (i in 0..<M) {
        if (Utils.isInf(mu[i, s[i].toInt()])) {
            for (j in mu.numCols - 1 downTo 0) {
                if (java.lang.Double.isFinite(mu[i, j])) {
                    s[i] = j.toDouble()
                    break
                }
            }
        }
        for (j in 0..<y.numCols) {
            y[i, j] = y[i, j] / mu[i, s[i].toInt()]
        }
    }
    for (i in 0..<M) {
        for (j in 0..<gamma.numCols) {
            gamma[i, j] = mu[i, j] / mu[i, s[i].toInt()]
        }
        var max = Double.NaN
        var Ncum = 0.0
        for (j in 0..<mu.numCols) {
            Ncum += 1.0
            val x = FastMath.abs(mu[i, j] - (Ncum))
            if (java.lang.Double.isNaN(max) || x > max) {
                max = x
            }
        }
        if (max < options.tol) {
            isDelay[i] = 1.0
        }
    }
    //% eliminating the delays seems to produce problems
    //% Z = sum([Z; L(isDelay,:)],1);
    //% L(isDelay,:)=[];
    //% mu(isDelay,:)=[];
    //% gamma(isDelay,:)=[];
    //% y(isDelay,:)=[];
    //% isLI(isDelay) = [];
    //% M = M - sum(isDelay);
    val beta = Matrix.ones(M, FastMath.ceil(N.elementSum()).toInt())
    for (i in 0..<M) {
        var beta_ij = gamma[i, 0] / (1 - gamma[i, 0])
        beta[i, 0] = if (java.lang.Double.isNaN(beta_ij)) Inf else beta_ij
        var j = 1
        while (j < FastMath.ceil(N.elementSum())) {
            beta_ij = (1 - gamma[i, j - 1]) * (gamma[i, j] / (1 - gamma[i, j]))
            beta[i, j] = if (java.lang.Double.isNaN(beta_ij)) Inf else beta_ij
            j++
        }
    }
    var isInf = true
    var idx = 0
    while (isInf && idx < beta.numElements) {
        if (java.lang.Double.isFinite(beta[idx])) {
            isInf = false
        }
        idx++
    }
    if (isInf) {
        options.method = "default"
        val lG = pfqn_nc(lambda, L, N, Z, options).lG
        options.method = method
        return Ret.pfqnRd(lG)
    }
    var Cgamma = 0.0
    val sld: MutableList<Double> = ArrayList()
    for (i in 0..<s.numElements) {
        if (s[i] > 0) {
            sld.add(s[i])
        }
    }
    val sldM = Matrix(sld)
    val vmax = FastMath.min(sldM.elementSum().toInt(), FastMath.ceil(N.elementSum()).toInt())
    val Y = pfqn_mva(y, N, Matrix.scaleMult(N, 0.0), Matrix.ones(1, M)).X
    val rhoN = y.mult(Y.transpose())
    val lEN = Matrix(1, vmax + 1)
    lEN.zero()
    for (vtot in 0..<vmax) {
        lEN[vtot + 1] =
            FastMath.log(abs(pfqn_gldsingle(rhoN, Matrix.singleton((vtot + 1).toDouble()), beta, options).G))
    }
    var EN: Double
    for (vtot in 0..<lEN.numElements) {
        EN = FastMath.exp(lEN[vtot])
        Cgamma += ((N.elementSum() - Maths.max(0.0, (vtot - 1).toDouble())) / N.elementSum()) * EN
    }
    options.method = "default"
    var lGN = pfqn_nc(lambda, y, N, Z, options).lG
    options.method = method
    lGN += FastMath.log(Cgamma)
    return Ret.pfqnRd(lGN, Cgamma)
}
/**
 * PFQN rd algorithms
 */
@Suppress("unused")
class PfqnRdAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}