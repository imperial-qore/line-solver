package jline.api.cache

import jline.io.Ret
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath
import kotlin.math.sqrt

/**
 * Approximate the normalizing constant of the cache steady state distribution using the ray method.
 *
 * @param gamma - Matrix representing the cache access factors.
 * @param m     - Matrix representing the cache capacity vector.
 * @return cacheRayIntReturn - the approximated normalizing constant (Z, and its logarithm lE) and the
 * xi terms.
 */

fun cache_rayint(gamma: Matrix, m: Matrix): Ret.cacheSpm {
    var gamma = gamma
    val rowsToKeep = BooleanArray(gamma.numRows)
    val colsToKeep = BooleanArray(gamma.numCols)
    for (i in 0..<gamma.numCols) {
        colsToKeep[i] = true
    }
    for (i in 0..<gamma.numRows) {
        if (gamma.getRow(i).elementSum() > 0) {
            rowsToKeep[i] = true
        } else {
            rowsToKeep[i] = false
        }
    }
    gamma = gamma.getSlice(rowsToKeep, colsToKeep)

    val h = m.numElements
    val n = gamma.numRows
    val mt = m.elementSum()

    if (n.toDouble() == mt) {
        System.err.println("The number of items equals the cache capacity.")
    }

    val xi = cache_xi_bvh(gamma, m)

    val S = Matrix(n, 1)
    for (k in 0..<n) {
        var Sk = 0.0
        for (l in 0..<h) {
            Sk += gamma[k, l] * xi[l]
        }
        S[k, 0] = Sk
    }

    // phi
    var phi = 0.0
    for (k in 0..<n) {
        phi += FastMath.log(1 + S[k])
    }
    phi -= xi.copy().log().mult(m.copy().transpose()).elementSum()

    // A
    val delta = Matrix.eye(h)
    val C = Matrix(h, h)
    for (j in 0..<h) {
        for (l in 0..<h) {
            var C1 = 0.0
            for (k in 0..<n) {
                C1 += gamma[k, j] / (1 + S[k])
            }
            var C2 = 0.0
            for (k in 0..<n) {
                C2 += gamma[k, j] * gamma[k, l] / FastMath.pow(1 + S[k], 2)
            }
            C[j, l] = delta[j, l] * C1 - xi[j] * C2
        }
    }

    // Z
    val Z = FastMath.exp(phi) * FastMath.pow(sqrt(2 * FastMath.PI), -h) * m.fact().elementMult() / xi.sqrt()
        .elementMult() / FastMath.sqrt(C.det())
    val lZ = -h * FastMath.log(sqrt(2 * FastMath.PI)) + phi + m.factln().elementSum() - xi.sqrt().log()
        .elementSum() - FastMath.log(sqrt(C.det()))

    return Ret.cacheSpm(Z, lZ, xi)
}