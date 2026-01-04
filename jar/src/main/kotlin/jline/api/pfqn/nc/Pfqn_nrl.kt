/**
 * @file Normal distribution-based Random Lattice (NRL) method for normalizing constants
 * 
 * Implements the Normal Random Lattice approach for computing normalizing constants
 * using Laplace approximation with complex-valued service demands. Provides high-accuracy
 * approximation for closed networks through sophisticated lattice-based integration.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.nc

import jline.GlobalConstants.NegInf
import jline.api.pfqn.ld.pfqn_gld
import jline.api.pfqn.ld.pfqn_gld_complex
import jline.solvers.SolverOptions
import jline.util.Maths
import jline.util.SerializableFunction
import jline.util.matrix.ComplexMatrix
import jline.util.matrix.Matrix
import org.apache.commons.math3.complex.Complex
import org.apache.commons.math3.util.FastMath
import kotlin.math.exp

fun pfqn_nrl(L: Matrix, N: Matrix, Z: Matrix, alpha: Matrix, options: SolverOptions?): Double {
    var L = L
    var alpha = alpha
    val Nt = N.elementSum()
    if (Nt < 0) {
        return NegInf
    }
    if (Nt == 0.0) {
        return 0.0
    }
    val M = L.numRows
    val R = L.numCols
    if (Z.elementSum() < 0) {
        L = Matrix.concatRows(L, Z, null)
        val Ntrange: MutableList<Double> = ArrayList()
        var i = 1.0
        while (i < Nt) {
            Ntrange.add(i)
            i++
        }
        alpha = Matrix.concatRows(alpha, Matrix(Ntrange), null)
    }
    if (M == 1 && Z.elementSum() == 0.0) {
        return pfqn_gld(L, N, alpha, options).lG
    }

    val Lmax = L.elementMax()
    val Lmax_mat = Matrix.singleton(Lmax).repmat(1, R)
    val L_final = L.elementDivide(Lmax_mat.repmat(M, 1))
    val x0 = Matrix(1, R)
    x0.zero()
    val lG = Maths.laplaceapprox_h_complex(x0, infradius_h(L_final, N, alpha)).logI.real
    val Lmax_log_mat = Lmax_mat.copy()
    for (i in 0..<Lmax_mat.numElements) {
        Lmax_log_mat[i] = FastMath.log(Lmax_mat[i])
    }
    return lG + N.mult(Lmax_log_mat.transpose())[0]
}

fun infradius_h(L: Matrix, N: Matrix, alpha: Matrix): SerializableFunction<Matrix, ComplexMatrix> {
    return object : SerializableFunction<Matrix, ComplexMatrix> {
        override fun apply(x: Matrix): ComplexMatrix {
            L.numRows
            val Nt = N.elementSum()
            val beta = Matrix(N.numRows, N.numCols)
            N.divide(Nt, beta, true)
            val t = Matrix(x.numRows, x.numCols)
            val tbMat = Matrix(x.numRows, 1)

            for (i in 0 until x.numRows) {
                var tbRow = 0.0
                for (j in 0 until x.numCols) {
                    val expVal = exp(x[i, j])
                    val tVal = expVal / (1.0 + expVal)
                    t[i, j] = tVal
                    tbRow += beta[i, j] * tVal
                }
                tbMat[i, 0] = tbRow
            }

            val tb = tbMat.elementSum()
            val tsubtb = t.copy()
            for (i in 0 until t.numElements) {
                tsubtb[i] = tsubtb[i] - tb
            }

            val h = nrl_h(L, tsubtb, Nt, alpha)
            val y = ComplexMatrix(Matrix(x.numRows, 1), Matrix(x.numRows, 1))
            for (i in 0 until x.numRows) {
                val xi = Matrix.extractRows(x, i, x.numRows, null)
                y[i] = h(xi)
            }

            return y
        }
    }
}

fun nrl_h(L: Matrix, tsubtb: Matrix, Nt: Double, alpha: Matrix): (Matrix) -> Complex {
    return fun(x: Matrix): Complex {
        val M = L.numRows
        val c = ComplexMatrix(tsubtb.numRows, tsubtb.numCols)

        for (i in 0 until c.numElements) {
            c[i] = Complex(0.0, 2 * Math.PI * tsubtb[i]).exp()
        }

        val LComplex = ComplexMatrix(L.elementMult(c.real.repmat(M, 1), null), L.elementMult(c.im.repmat(M, 1), null))

        val expX = x.copy()
        for (i in 0 until x.numElements) {
            val expVal = exp(x[i])
            expX[i] = expVal / ((1 + expVal) * (1 + expVal))
        }

        return pfqn_gld_complex(LComplex.sumRows(), Matrix.singleton(Nt), alpha, null).G.multiply(expX.elementMult())
    }
}
/**
 * PFQN nrl algorithms
 */
@Suppress("unused")
class PfqnNrlAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}