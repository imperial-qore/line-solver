/**
 * @file Complex-valued generalized load-dependent normalizing constant computation
 * 
 * Extends the generalized load-dependent convolution algorithm to handle complex-valued
 * service demands and parameters. Supports analysis of queueing networks with complex
 * service characteristics while maintaining load-dependent service rate capabilities.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.ld

import jline.io.Ret
import jline.GlobalConstants
import jline.GlobalConstants.NegInf
import jline.solvers.SolverOptions
import jline.solvers.nc.SolverNC
import jline.util.Maths
import jline.util.matrix.ComplexMatrix
import jline.util.matrix.Matrix
import org.apache.commons.math3.complex.Complex
import org.apache.commons.math3.util.FastMath
import kotlin.math.atan
import kotlin.math.ln
import kotlin.math.pow

/**
 * Compute the normalizing constant of a single-class load-dependent closed queueing network model with complex
 * demands.
 *
 * @param L       - demands at all stations
 * @param N       - number of jobs for each class
 * @param mu      - load-depedent scalings
 * @param options - solver options
 * @return normalizing constant and its logarithm
 */

fun pfqn_gld_complex(L: ComplexMatrix, N: Matrix, mu: Matrix?, options: SolverOptions?): Ret.pfqnNcComplex {
    val M = L.numRows
    val R = L.numCols
    Matrix(1, R)
    var G: Complex
    val lG: Complex

    if (M == 1) {
        var N_tmp = Matrix(1, 0)
        val L_tmp = ComplexMatrix(1, 0)
        for (i in 0..<R) {
            if (L.real[i] > GlobalConstants.FineTol || L.im[i] > GlobalConstants.FineTol) {
                val N_tmp2 = Matrix(1, 1)
                N_tmp2.fill(N[i])
                val L_tmp2 = ComplexMatrix(1, 1)
                L_tmp2.real.fill(0.5 * FastMath.log(L.real[0, i].pow(2.0) + FastMath.pow(L.im[0, i], 2)))
                L_tmp2.im.fill(ln(atan(L.im[0, i] / L.real[0, i])))
                N_tmp = Matrix.concatColumns(N_tmp, N_tmp2, null)
                L_tmp.real = Matrix.concatColumns(L_tmp.real, L_tmp2.real, null)
                L_tmp.im = Matrix.concatColumns(L_tmp.im, L_tmp2.im, null)
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
        lG = Complex(N_tmp.mult(L_tmp.real.transpose())[0],
            N_tmp.mult(L_tmp.im.transpose())[0]).add(Maths.factln(N.elementSum()) - Matrix.factln(N).elementSum())
            .subtract(mu_new.elementSum())
        G = lG.exp()
        return Ret.pfqnNcComplex(G, lG)
    }

    if (R == 1) {
        val ret = pfqn_gldsingle_complex(L, N, mu!!, null)
        lG = ret.lG
        G = ret.G
        return Ret.pfqnNcComplex(G, lG)
    }

    if (L.isEmpty) {
        G = Complex(0.0)
        lG = Complex(NegInf)
        return Ret.pfqnNcComplex(G, lG)
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
        throw RuntimeException("pfqn_gld_complex is only implemented for load dependent models.")
    }

    G = Complex(0.0)
    if (M == 0) {
        lG = G.log()
        return Ret.pfqnNcComplex(G, lG)
    }

    if (FastMath.abs(N.elementMax()) < GlobalConstants.FineTol && FastMath.abs(N.elementMin()) < GlobalConstants.FineTol) {
        G = Complex(1.0)
        lG = G.log()
        return Ret.pfqnNcComplex(G, lG)
    }

    if (R == 1) {
        return pfqn_gldsingle_complex(L, N, mu_new, null)
    }

    G = pfqn_gld_complex(ComplexMatrix.extractRows(L, 0, M - 1, null),
        N,
        Matrix.extractRows(mu_new, 0, M - 1, null),
        options_new).G

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
            G = G.add(L[M - 1, r].divide(mu_new[M - 1, 0])
                .multiply(pfqn_gld_complex(L, N_1, pfqn_mushift(mu!!, M - 1), options_new).G))
        }
    }
    lG = G.log()
    return Ret.pfqnNcComplex(G, lG)
}
/**
 * PFQN gld complex algorithms
 */
@Suppress("unused")
class PfqnGldComplexAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}