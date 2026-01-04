/**
 * @file Logistic sampling method for normalizing constant computation
 * 
 * Implements the logistic sampling approach for computing normalizing constants in
 * closed product-form queueing networks. Uses importance sampling with multivariate
 * normal distributions centered at the Laguerre expansion maximum.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.nc

import jline.io.Ret
import jline.GlobalConstants
import jline.util.Maths
import jline.util.matrix.Matrix
import org.apache.commons.math3.distribution.MultivariateNormalDistribution
import org.apache.commons.math3.util.FastMath
import kotlin.math.exp

/**
 * Logistic sampling method to compute the normalizing constant
 *
 * @param L - demands at all stations
 * @param N - number of jobs for each class
 * @param Z - think time for each class
 * @param I - number of samples
 * @return normalizing constant and its logarithm
 */
fun pfqn_ls(L: Matrix, N: Matrix, Z: Matrix = Matrix(1, L.numCols), I: Long, seed: Long = 23000): Ret.pfqnNc {
    var M = L.numRows
    var R = L.numCols
    val Lsum = Matrix(M, 1)
    for (i in 0..<M) {
        Lsum[i] = L.sumRows(i)
    }
    var L_new = Matrix(0, R)
    for (i in 0..<M) {
        val L_row_i = Matrix(1, R)
        Matrix.extract(L, i, i + 1, 0, R, L_row_i, 0, 0)
        if (Lsum[i] > GlobalConstants.CoarseTol) {
            L_new = Matrix.concatRows(L_new, L_row_i, null)
        }
    }
    M = L_new.numRows
    R = L_new.numCols
    var sample: Array<DoubleArray>? = null
    val lGn: Double

    if (L_new.isEmpty || N.isEmpty || N.elementSum() == 0.0 || L_new.elementSum() < GlobalConstants.CoarseTol) {
        val tmp = Matrix(1, Z!!.numCols)
        for (i in 0..<tmp.length()) {
            tmp[i] = FastMath.log(Z.sumCols(i))
        }
        lGn = (-Matrix.factln(N).elementSum() + N.elementMult(tmp, null).elementSum())
    } else if (Z.isEmpty) {
        val ret = pfqn_le_fpi(L_new, N)
        val umax = ret.u
        val A = pfqn_le_hessian(L_new, N, umax.transpose())
        val A_t = A.transpose()
        for (i in 0..<A.numRows) {
            for (j in 0..<A.numCols) {
                A[i, j] = (A[i, j] + A_t[i, j]) / 2.0
            }
        }
        val iA = A.inv()
        val x0 = Matrix(1, M - 1)
        for (i in 0..<M - 1) {
            x0[i] = FastMath.log(umax[i] / umax[M - 1])
        }

        val x0_array = DoubleArray(M - 1)
        for (i in 0..<M - 1) {
            x0_array[i] = x0[i]
        }
        val iA_array = arrayOfNulls<DoubleArray>(iA.numRows)
        for (i in 0..<iA.numRows) {
            val tmp_row = DoubleArray(iA.numCols)
            for (j in 0..<iA.numCols) {
                tmp_row[j] = iA[i, j]
            }
            iA_array[i] = tmp_row
        }
        val mvd = MultivariateNormalDistribution(x0_array, iA_array)
        mvd.reseedRandomGenerator(seed)

        if (sample == null) {
            sample = mvd.sample(I.toInt())
        }

        val T = Matrix(I.toInt(), 1)
        val dpdf = Matrix(1, I.toInt())
        for (i in 0..<I) {
            val sample_i = sample!![i.toInt()]
            T[i.toInt()] = Maths.simplex_fun(sample_i, L_new, N)
            dpdf[i.toInt()] = mvd.density(sample_i)
        }

        val tmp = Matrix(1, N.length() + 1)
        Matrix.extract(N, 0, 1, 0, N.length(), tmp, 0, 0)
        tmp[N.length()] = (M - 1).toDouble()
        var div_sum = 0.0
        for (i in 0..<I) {
            div_sum += T[i.toInt()] / dpdf[i.toInt()]
        }
        lGn = Maths.multinomialln(tmp) + Maths.factln(M - 1) + FastMath.log(div_sum / I)
    } else {
        val ret = pfqn_le_fpiZ(L_new, N, Z)
        val umax = ret.u
        val vmax = ret.v
        val A = pfqn_le_hessianZ(L_new, N, Z, umax.transpose(), vmax)
        val A_t = A.transpose()
        for (i in 0..<A.numRows) {
            for (j in 0..<A.numCols) {
                A[i, j] = (A[i, j] + A_t[i, j]) / 2.0
            }
        }
        val iA = A.inv()
        val x0 = Matrix(1, M)
        for (i in 0..<M - 1) {
            x0[i] = FastMath.log(umax[i] / umax[M - 1])
        }
        x0[M - 1] = FastMath.log(vmax)

        val x0_array = DoubleArray(M)
        for (i in 0..<M) {
            x0_array[i] = x0[i]
        }
        val iA_array = arrayOfNulls<DoubleArray>(iA.numRows)
        for (i in 0..<iA.numRows) {
            val tmp_row = DoubleArray(iA.numCols)
            for (j in 0..<iA.numCols) {
                tmp_row[j] = iA[i, j]
            }
            iA_array[i] = tmp_row
        }
        val mvd = MultivariateNormalDistribution(x0_array, iA_array)
        mvd.reseedRandomGenerator(seed)

        if (sample == null) {
            sample = mvd.sample(I.toInt())
        }

        val T = Matrix(I.toInt(), 1)
        val epsilon = 1e-10
        val eN = epsilon * N.elementSum()
        val eta = N.elementSum() + M * (1 + eN)
        val K = M
        val dpdf = Matrix(1, I.toInt())

        for (i in 0..<I) {
            val sample_i = sample!![i.toInt()]
            T[i.toInt()] = pfqn_ls_helper(sample_i, K, M, eta, eN, L_new, N, Z)
            dpdf[i.toInt()] = mvd.density(sample_i)
        }

        var div_sum = 0.0
        for (i in 0..<I) {
            div_sum += T[i.toInt()] / dpdf[i.toInt()]
        }
        lGn = FastMath.log(exp(-Matrix.factln(N).elementSum()) * div_sum / I)
    }
    val Gn = FastMath.exp(lGn)
    return Ret.pfqnNc(Gn, lGn)
}


/**
 * Auxiliary function used in the logistic sampling method
 */
internal fun pfqn_ls_helper(x: DoubleArray,
                            K: Int,
                            M: Int,
                            eta: Double,
                            eN: Double,
                            L: Matrix,
                            N: Matrix,
                            Z: Matrix): Double {
    var res = -exp(x[K - 1]) + K * (1 + eN) * x[M - 1]
    var tmp1 = 0.0
    var tmp2 = 0.0
    for (i in 0..<K - 1) {
        tmp1 += x[i]
        tmp2 += FastMath.exp(x[i])
    }
    tmp2 = FastMath.log(tmp2 + 1)
    tmp2 *= -eta
    res += (tmp1 + tmp2)

    val L_row_K = Matrix(1, L.numCols)
    val L_first_K_minus_one_rows = Matrix(K - 1, L.numCols)
    var x_first_K_minus_one_elements = Matrix(1, K - 1)
    Matrix.extract(L, K - 1, K, 0, L.numCols, L_row_K, 0, 0)
    Matrix.extract(L, 0, K - 1, 0, L.numCols, L_first_K_minus_one_rows, 0, 0)
    for (i in 0..<K - 1) {
        x_first_K_minus_one_elements[i] = FastMath.exp(x[i])
    }
    for (i in 0..<K - 1) {
        for (j in 0..<L.numCols) {
            L_first_K_minus_one_rows[i, j] = L_first_K_minus_one_rows[i, j] * FastMath.exp(x[K - 1]) + Z[j]
        }
    }
    x_first_K_minus_one_elements = x_first_K_minus_one_elements.mult(L_first_K_minus_one_rows)
    for (i in 0..<L_row_K.length()) {
        L_row_K[i] = FastMath.log(L_row_K[i] * FastMath.exp(x[K - 1]) + Z[i] + x_first_K_minus_one_elements[i])
    }
    res += N.mult(L_row_K.transpose()).elementSum()

    res = FastMath.exp(res)
    return res
}
/**
 * PFQN ls algorithms
 */
@Suppress("unused")
class PfqnLsAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}