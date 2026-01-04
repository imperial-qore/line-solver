/**
 * @file Gauss-Legendre integration for multi-class repairman model normalizing constants
 * 
 * Implements Gauss-Legendre quadrature integration for computing normalizing constants
 * in multi-class repairman models. Uses precomputed Gauss-Legendre nodes and weights
 * for high-precision numerical integration over infinite intervals.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.nc

import jline.io.Ret
import jline.util.Maths
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath
import java.io.BufferedReader
import java.io.FileNotFoundException
import java.io.IOException
import java.io.InputStreamReader

/**
 * Compute the normalizing constant of a repairmen model using Gauss-Legendre integration
 *
 * @param L - demands at all stations
 * @param N - number of jobs for each class
 * @param Z - think times
 * @return normalizing constant and its logarithm
 */
fun pfqn_mmint2_gausslegendre(L: Matrix, N: Matrix, Z: Matrix, m: Int?): Ret.pfqnNc {
    var m = m
    if (m == null) {
        m = 1
    }

    val gausslegendreNodes: MutableList<Double> = ArrayList()
    val gausslegendreWeights: MutableList<Double> = ArrayList()

    try {
        val nodeStream = object {}.javaClass.getResourceAsStream("/gausslegendre-nodes.txt")
            ?: throw FileNotFoundException("Resource gausslegendre-nodes.txt not found.")
        val nodeReader = BufferedReader(InputStreamReader(nodeStream))
        var line: String
        while ((nodeReader.readLine().also { line = it }) != null) {
            gausslegendreNodes.add(line.toDouble())
        }

        val weightStream = object {}.javaClass.getResourceAsStream("/gausslegendre-weights.txt")
        val weightReader = BufferedReader(InputStreamReader(weightStream))
        while ((weightReader.readLine().also { line = it }) != null) {
            gausslegendreWeights.add(line.toDouble())
        }
    } catch (e1: FileNotFoundException) {
        e1.printStackTrace()
    } catch (e: IOException) {
        throw RuntimeException(e)
    }

    val n = FastMath.max(300.0,
        FastMath.min(gausslegendreNodes.size.toDouble(), 2 * (N.sumRows().sumCols()[0] + m - 1) - 1)).toInt()
    val y = Matrix(1, n)
    y.fill(0.0)

    if (!(Z.numRows == L.numRows && Z.numCols == L.numCols)) {
        throw RuntimeException("The dimensions of Z and L are not the same.")
    }
    for (i in 0..<n) {
        val tmp = L.copy()

        for (j in 0..<tmp.numRows) {
            for (k in 0..<tmp.numCols) {
                tmp[j, k] = FastMath.log(Z[j, k] + gausslegendreNodes[i] * tmp[j, k])
            }
        }
        y[i] = (N.mult(tmp.transpose())).value()
    }

    val g = y.copy()
    val nodes = Matrix(1, n)
    val logNodes = Matrix(1, n)
    val logWeights = Matrix(1, n)
    for (i in 0..<n) {
        nodes[i] = gausslegendreNodes[i]
        logNodes[i] = FastMath.log(gausslegendreNodes[i])
        logWeights[i] = FastMath.log(gausslegendreWeights[i])
    }

    for (i in 0..<n) {
        g[i] = g[i] + logWeights[i] - nodes[i]
    }

    var coeff = 0.0
    for (i in 0..<N.length()) {
        coeff -= Maths.factln(N[i])
    }
    coeff -= Maths.factln(m - 1)
    coeff += (m - 1) * logNodes.elementSum()

    var lG = 0.0
    for (i in 0..<g.length()) {
        lG += FastMath.exp(g[i])
    }
    lG = FastMath.log(lG) + coeff
    if (!java.lang.Double.isFinite(lG)) {
        lG = Matrix.logsumexp(g) + coeff
    }
    val G = FastMath.exp(lG)
    return Ret.pfqnNc(G, lG)
}
/**
 * PFQN mmint2 gausslegendre algorithms
 */
@Suppress("unused")
class PfqnMmint2GausslegendreAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}