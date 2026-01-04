/**
 * @file Normal distribution Random Permutation (NRP) method for normalizing constants
 * 
 * Implements the Normal Random Permutation approach for computing normalizing constants
 * using cumulative distribution function transformations. Combines normal distribution
 * sampling with Laplace approximation for enhanced accuracy in closed networks.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.nc

import jline.api.pfqn.ld.pfqn_gld
import jline.api.pfqn.ld.pfqn_gld_complex
import jline.solvers.SolverOptions
import jline.util.Maths
import jline.util.SerializableFunction
import jline.util.matrix.ComplexMatrix
import jline.util.matrix.Matrix
import org.apache.commons.math3.complex.Complex
import org.apache.commons.math3.distribution.NormalDistribution
import org.apache.commons.math3.util.FastMath

fun pfqn_nrp(L: Matrix, N: Matrix, Z: Matrix, alpha: Matrix, options: SolverOptions?): Double {
    var L = L
    var alpha = alpha
    val Nt = N.elementSum()
    if (Z.elementSum() > 0) {
        L = Matrix.concatRows(L, Z, null)
        val alpha_tmp = Matrix(1, Nt.toInt())
        var i = 0
        while (i < Nt + 1) {
            alpha_tmp[1, i] = i + 1
            i++
        }
        alpha = Matrix.concatRows(alpha, alpha_tmp, null)
    }
    val M = L.numRows
    val R = L.numCols
    if (M == 1 && Z.elementSum() == 0.0) {
        return pfqn_gld(L, N, alpha, options).lG
    }
    val Lmax = L.elementMax()
    val Lmax_mat = Matrix.singleton(Lmax).repmat(1, R)
    val L_final = L.elementDivide(Lmax_mat.repmat(M, 1))
    val x0 = Matrix(1, R)
    x0.zero()
    val lG = Maths.laplaceapprox_h_complex(x0, infradius_hnorm(L_final, N, alpha)).logI.real
    val Lmax_log_mat = Lmax_mat.copy()
    for (i in 0..<Lmax_mat.numElements) {
        Lmax_log_mat[i] = FastMath.log(Lmax_mat[i])
    }
    return lG + N.mult(Lmax_log_mat.transpose())[0]
}

fun infradius_hnorm(L: Matrix, N: Matrix, alpha: Matrix): SerializableFunction<Matrix, ComplexMatrix> {
    return object : SerializableFunction<Matrix, ComplexMatrix> {
        override fun apply(x: Matrix): ComplexMatrix {
            val MU = 0.0
            val SIGMA = 1.0
            val Nt = N.elementSum()
            val beta = Matrix(N.numRows, N.numCols)
            N.divide(Nt, beta, true)
            val Z = NormalDistribution(MU, SIGMA)

            val t = x.copy()
            for (i in 0 until t.numRows) {
                for (j in 0 until t.numCols) {
                    t[i, j] = Z.cumulativeProbability(x[i, j])
                }
            }

            var tb = 0.0
            for (i in 0 until t.numRows) {
                for (j in 0 until t.numCols) {
                    tb += beta[i, j] * t[i, j]
                }
            }

            val tsubtb = t.copy()
            for (i in 0 until t.numElements) {
                tsubtb[i] = tsubtb[i] - tb
            }

            val y = ComplexMatrix(Matrix(x.numRows, 1), Matrix(x.numRows, 1))
            val h = nrp_h(L, tsubtb, Nt, alpha)
            for (i in 0 until x.numRows) {
                val xi = Matrix.extractRows(x, i, x.numRows, null)
                y[i] = h(xi)
            }

            return ComplexMatrix(y.real)
        }
    }
}

fun nrp_h(L: Matrix, tsubtb: Matrix, Nt: Double, alpha: Matrix): (Matrix) -> Complex {
    return fun(x: Matrix): Complex {
        val M = L.numRows
        val c = ComplexMatrix(tsubtb.numRows, tsubtb.numCols)

        for (i in 0 until c.numElements) {
            c[i] = Complex(0.0, 2 * Math.PI * tsubtb[i]).exp()
        }

        val LComplex = ComplexMatrix(L.elementMult(c.real.repmat(M, 1), null), L.elementMult(c.im.repmat(M, 1), null))

        val Z = NormalDistribution()
        val normpdfX = x.copy()
        for (i in 0 until x.numElements) {
            normpdfX[i] = Z.density(x[i])
        }

        return pfqn_gld_complex(
            LComplex.sumRows(),
            Matrix.singleton(Nt),
            alpha,
            null
        ).G.multiply(normpdfX.elementMult())
    }
}
/**
 * PFQN nrp algorithms
 */
@Suppress("unused")
class PfqnNrpAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}