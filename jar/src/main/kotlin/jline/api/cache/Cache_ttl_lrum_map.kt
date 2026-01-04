/**
 * @file TTL-based LRUM-MAP Cache Analysis
 * 
 * Combines TTL approximation with LRUM cache policies and Markovian Arrival
 * Process (MAP) input streams. Provides sophisticated analysis for multi-server
 * cache systems with correlated arrival patterns and TTL-based eviction.
 * 
 * @since LINE 3.0
 */
package jline.api.cache

import jline.api.mam.map_lambda
import jline.api.mam.map_piq
import jline.util.Maths
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

// currently D0k = -pk
fun cache_ttl_lrum_map(D0Matrix: Array<MatrixCell>, D1Matrix: Array<MatrixCell>, m: Matrix): Double {
    val n = D0Matrix.size

    val h = D0Matrix[0].size()
    val d = D0Matrix[0][0].numRows


    val t = cache_t_lrum_map(D0Matrix, D1Matrix, m)


    val probh = Matrix(n, h) // steady state 1,...,h
    val prob0 = Matrix(n, 1) // steady state 0
    val expD0 = arrayOfNulls<MatrixCell>(n)
    val trans = arrayOfNulls<MatrixCell>(n)
    val pih = arrayOfNulls<MatrixCell>(n)
    val pi0 = arrayOfNulls<MatrixCell>(n)
    val R = arrayOfNulls<MatrixCell>(n)
    val Nh = arrayOfNulls<MatrixCell>(n)
    val N0 = arrayOfNulls<MatrixCell>(n)
    arrayOfNulls<MatrixCell>(n)
    val Q = arrayOfNulls<MatrixCell>(n)
    val lambda = Matrix(n, 1)
    var lambda_sum = 0.0



    for (i in 0..<n) {
        expD0[i] = MatrixCell(h)
        trans[i] = MatrixCell(h)
        pih[i] = MatrixCell(h)
        pi0[i] = MatrixCell(1)
        R[i] = MatrixCell(h)
        Nh[i] = MatrixCell(h)
        N0[i] = MatrixCell(1)
        Q[i] = MatrixCell(h)
        pi0[i]!![0] = map_piq(D0Matrix[i][0], D1Matrix[i][0])
        val current_lambda = (pi0[i]!![0].mult(D1Matrix[0][0])).mult(Matrix.ones(d, 1))
        lambda[i] = current_lambda[0]
        lambda_sum += lambda[i]
    }


    for (k in 0..<n) {
        for (l in 0..<h) {
            val index = Matrix.scaleMult(D0Matrix[k][l], t[l])
            val expD0kl = Maths.matrixExp(index)
            expD0[k]!![l] = expD0kl
            val temp = (Matrix.negative(D0Matrix[k][l])).inv()
            Nh[k]!![l] = Matrix.oneMinusMatrix(expD0kl).mult(temp)
            trans[k]!![l] = Matrix.oneMinusMatrix(expD0kl).mult(temp).mult(D1Matrix[k][l])
        }

        for (l in h - 1 downTo 0) {
            if (l == h - 1) {
                R[k]!![l] = trans[k]!![l - 1].mult(Matrix.oneMinusMatrix(trans[k]!![l]).inv())
            } else if (l == 0) {
                //R[k].set(l, Matrix.oneMinusMatrix(R[k].get(l+1).mult(expD0[k].get(l+1))).inv());
                val temp = (Matrix.negative(D0Matrix[k][l])).inv().mult(D1Matrix[k][l]) //(A0)
                R[k]!![l] = temp.mult(Matrix.oneMinusMatrix(R[k]!![l + 1].mult(expD0[k]!![l + 1])).inv())
            } else {
                R[k]!![l] = trans[k]!![l - 1].mult(Matrix.oneMinusMatrix(R[k]!![l + 1].mult(expD0[k]!![l + 1])).inv())
            }
            N0[k]!![0] = (Matrix.negative(D0Matrix[k][0])).inv()
        }
    }

    val prod_R = arrayOfNulls<MatrixCell>(n)
    val sum_prod_R = MatrixCell()

    for (k in 0..<n) {
        sum_prod_R[k] = R[k]!![1]
        prod_R[k] = MatrixCell(h + 1)
        prod_R[k]!![0] = R[k]!![1]
        for (l in 1..<h) {
            prod_R[k]!![l] = prod_R[k]!![l - 1].mult(R[k]!![l])
            sum_prod_R[k] = sum_prod_R[k].add(prod_R[k]!![l])
        }
    }


    for (k in 0..<n) {
        pih[k]!![0] = pi0[k]!![0].mult(R[k]!![0])
        for (l in 1..<h) {
            pih[k]!![l] = pih[k]!![l - 1].mult(R[k]!![l])
        }
    }


    val e = Matrix(d, 1)
    e.ones()
    val lambdaMatrix = Array(n) { DoubleArray(h) }
    val lambdaMatrixsum = DoubleArray(h)
    for (l in 0..<h) {
        lambdaMatrix[0][l] = (pih[0]!![l].mult(D1Matrix[0][l]).mult(e))[0, 0]
        lambdaMatrixsum[l] = lambdaMatrix[0][l]
        for (k in 1..<n) {
            lambdaMatrix[k][l] = (pih[k]!![l].mult(D1Matrix[k][l]).mult(e))[0, 0]
            lambdaMatrixsum[l] = lambdaMatrixsum[l] + lambdaMatrix[k][l]
        }
    }


    val production = Array(n) { DoubleArray(h + 1) }
    val productionD = Array(n) { DoubleArray(h + 1) }
    for (k in 0..<n) {
        production[k][0] = (pi0[k]!![0]).mult(N0[k]!![0]).mult(e)[0, 0]
        productionD[k][0] = (pi0[k]!![0]).mult(N0[k]!![0]).mult(e)[0, 0]
    }
    for (l in 1..<h + 1) {
        for (k in 0..<n) {
            production[k][l] = (pih[k]!![l - 1]).mult(Nh[k]!![l - 1]).mult(e)[0, 0]
            productionD[k][l] = (pih[k]!![l - 1]).mult(Nh[k]!![l - 1]).mult(D1Matrix[k][l - 1]).mult(e)[0, 0]
        }
    }

    val m_matrixD = DoubleArray(h + 1)
    for (l in 1..<h + 1) {
        for (k in 0..<n) {
            var down_sum = 0.0
            for (j in 0..<h + 1) {
                down_sum += production[k][j]
            }
            probh[k, l - 1] = productionD[k][l] / down_sum
            m_matrixD[l] += productionD[k][l] / down_sum
        }
    }

    val hitrate = DoubleArray(h)
    for (k in 0..<n) {
        val probhK = Matrix(1, h)
        for (i in 0..<h) {
            probhK[0, i] = probh[k, i]
            hitrate[i] += probhK[0, i]
        }
        prob0[k, 0] = 1 - probhK.elementSum()
    }

    val arrival = DoubleArray(n)
    var sum_arrival = 0.0
    for (i in 0..<n) {
        arrival[i] = map_lambda(D0Matrix[i][0], D1Matrix[i][0])
        sum_arrival += arrival[i]
    }


    var avg_hitrate = 0.0
    for (i in 0..<h) {
        avg_hitrate += hitrate[i]
    }

    avg_hitrate = avg_hitrate / sum_arrival


    return avg_hitrate
}
/**
 * Cache ttl lrum map algorithms
 */
@Suppress("unused")
class CacheTtlLrumMapAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}