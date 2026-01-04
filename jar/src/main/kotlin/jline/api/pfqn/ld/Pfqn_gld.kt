/**
 * @file Generalized load-dependent normalizing constant computation
 * 
 * Implements the generalized convolution algorithm for computing normalizing constants in
 * load-dependent closed queueing networks. Handles multi-class systems with state-dependent
 * service rates, providing exact normalization for product-form networks with complex dependencies.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.ld

import jline.io.Ret
import jline.GlobalConstants
import jline.GlobalConstants.NegInf
import jline.api.pfqn.nc.pfqn_nc
import jline.solvers.SolverOptions
import jline.solvers.nc.SolverNC
import jline.util.Maths
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath
import kotlin.math.ln

/**
 * Compute the normalizing constant of a single-class load-dependent closed queueing network model
 *
 * @param L       - demands at all stations
 * @param N       - number of jobs for each class
 * @param mu      - load-depedent scalings
 * @param options - solver options
 * @return normalizing constant and its logarithm
 */

fun pfqn_gld(L: Matrix, N: Matrix, mu: Matrix?, options: SolverOptions?): Ret.pfqnNc {
    val M = L.numRows
    val R = L.numCols
    val lambda = Matrix(1, R)
    var G: Double?
    val lG: Double

    if (M == 1) {
        var N_tmp = Matrix(1, 0)
        var L_tmp = Matrix(1, 0)
        for (i in 0..<R) {
            if (L[i] > GlobalConstants.FineTol) {
                val N_tmp2 = Matrix(1, 1)
                N_tmp2.fill(N[i])
                val L_tmp2 = Matrix(1, 1)
                L_tmp2.fill(ln(L[0, i]))
                N_tmp = Matrix.concatColumns(N_tmp, N_tmp2, null)
                L_tmp = Matrix.concatColumns(L_tmp, L_tmp2, null)
            }
        }
        var mu_new: Matrix
        if (N.elementSum().toInt() >= mu!!.numCols) {
            mu_new = Matrix.extractRows(mu, 0, 1, null)
        } else {
            mu_new = Matrix(1, 0)
            var i = 0
            while (i < N.elementSum()) {
                val mu_col_i = Matrix(1, 1)
                Matrix.extract(mu, 0, 1, i, i + 1, mu_col_i, 0, 0)
                mu_new = Matrix.concatColumns(mu_new, mu_col_i, null)
                i++
            }
        }

        for (i in 0..<mu_new.length()) {
            mu_new[i] = FastMath.log(mu_new[i])
        }

        lG = (Maths.factln(N.elementSum()) - Matrix.factln(N)
            .elementSum() + N_tmp.mult(L_tmp.transpose())[0] - mu_new.elementSum())
        G = FastMath.exp(lG)
        return Ret.pfqnNc(G, lG)
    }

    if (R == 1) {
        val ret = pfqn_gldsingle(L, N, mu!!, null)
        lG = ret.lG
        G = ret.G
        return Ret.pfqnNc(G, lG)
    }

    if (L.isEmpty) {
        G = 0.0
        lG = NegInf
        return Ret.pfqnNc(G, lG)
    }

    val mu_new: Matrix
    if (mu == null) {
        mu_new = Matrix(M, N.elementSum().toInt())
        mu_new.fill(1.0)
    } else {
        mu_new = mu.copy()
    }
    val options_new = options ?: SolverNC.defaultOptions()

    var isLoadDep = false
    val isInfServer = BooleanArray(M)
    for (i in 0..<M) {
        val mu_row_i = Matrix(1, N.elementSum().toInt())
        Matrix.extract(mu_new, i, i + 1, 0, N.elementSum().toInt(), mu_row_i, 0, 0)
        var flag = true
        for (j in 0..<mu_row_i.numCols) {
            if (FastMath.abs(mu_row_i[j] - (j + 1)) > GlobalConstants.FineTol) {
                flag = false
                break
            }
        }

        if (FastMath.abs(mu_row_i.elementMin() - 1) < GlobalConstants.FineTol && FastMath.abs(mu_row_i.elementMax() - 1) < GlobalConstants.FineTol) {
            isInfServer[i] = false
        } else if (flag) {
            isInfServer[i] = true
        } else {
            isInfServer[i] = false
            isLoadDep = true
        }
    }

    if (!isLoadDep) {
        var Lli = Matrix(0, L.numCols)
        var Zli = Matrix(0, L.numCols)
        for (i in 0..<M) {
            val L_row_i = Matrix.extractRows(L, i, i + 1, null)
            if (isInfServer[i]) {
                Zli = Matrix.concatRows(Zli, L_row_i, null)
            } else {
                Lli = Matrix.concatRows(Lli, L_row_i, null)
            }
        }
        if (Lli.isEmpty) {
            Lli = N.copy()
            Lli.fill(0.0)
        }
        if (Zli.isEmpty) {
            Zli = N.copy()
            Zli.fill(0.0)
        }
        options_new.method = "exact"
        lG = pfqn_nc(lambda, Lli, N, Zli.sumCols(), options_new).lG
        G = FastMath.exp(lG)
        return Ret.pfqnNc(G, lG)
    }

    G = 0.0
    if (M == 0) {
        lG = FastMath.log(G!!)
        return Ret.pfqnNc(G, lG)
    }

    if (FastMath.abs(N.elementMax()) < GlobalConstants.FineTol && FastMath.abs(N.elementMin()) < GlobalConstants.FineTol) {
        G = 1.0
        lG = FastMath.log(G!!)
        return Ret.pfqnNc(G, lG)
    }

    if (R == 1) {
        G = pfqn_gldsingle(L, N, mu_new, null).G
        lG = FastMath.log(G!!)
        return Ret.pfqnNc(G, lG)
    }

    G = pfqn_gld(Matrix.extractRows(L, 0, M - 1, null), N, Matrix.extractRows(mu_new, 0, M - 1, null), options_new).G

    for (r in 0..<R) {
        if (N[r] > GlobalConstants.FineTol) {
            val N_1 = N.copy()
            if (R > 1) {
                N_1[r] = N_1[r] - 1
            } else {
                for (i in 0..<N_1.length()) {
                    N_1[i] = N_1[i] - 1
                }
            }
            G += (L[M - 1, r] / mu_new[M - 1, 0] * pfqn_gld(L, N_1, pfqn_mushift(mu!!, M - 1), options_new).G)
        }
    }
    lG = FastMath.log(G!!)
    return Ret.pfqnNc(G, lG)
}
/**
 * PFQN gld algorithms
 */
@Suppress("unused")
class PfqnGldAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}