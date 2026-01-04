/**
 * @file Cache Xi Terms via Gast-van Houdt Method
 * 
 * Computes cache xi terms using the iterative method from Gast-van Houdt
 * (SIGMETRICS 2015). Assumes monotone access factors with list index for
 * analyzing item distribution patterns in hierarchical cache systems.
 * 
 * @since LINE 3.0
 */
package jline.api.cache

import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath

/**
 * Computes the cache xi terms using the iterative method used in Gast-van Houdt, SIGMETRICS 2015.
 * This method calculates the xi values, which are important for understanding the distribution of items
 * in the cache. The script assumes (like the paper) that the access factors are monotone with the list index, so
 * it may not work with arbitrary (non-monotone) access costs.
 *
 * @param gamma - Matrix representing the cache access factors.
 * @param m     - Matrix representing the cache capacity vector.
 * @return Matrix - A matrix containing the computed xi terms.
 */

fun cache_xi_iter(gamma: Matrix, m: Matrix): Matrix {
    val n = gamma.numRows
    val f = m.scale(1.0 / n)
    val h = f.numCols

    val pp = Matrix(h + 1, n)
    pp.setRow(0, Matrix.ones(1, n))
    for (i in 0..<h) {
        pp.setRow(i + 1, gamma.getColumn(i).transpose())
    }

    val zOld = Matrix.zeros(1, h + 1)
    val z = Matrix.ones(1, h + 1)

    while (z.sub(zOld).elementMaxAbs() > FastMath.pow(10.0, -12) * zOld.elementMaxAbs()) {
        zOld.setTo(z)
        val temp = z.mult(pp).scale(n.toDouble())
        for (i in 0..<h) {
            val a = temp.sub(pp.getRow(i + 1).scale(z.scale(n.toDouble())[0, i + 1]))
            val Fi = pp.getRow(i + 1).elementDiv(pp.getRow(i + 1).scale(n.toDouble()).add(a)).elementSum()

            var ziMin: Double
            var ziMax: Double
            if (Fi > f[0, i]) {
                ziMin = 0.0
                ziMax = 1.0
            } else {
                ziMin = 1.0
                ziMax = 2.0
                while (pp.getRow(i + 1).scale(ziMax).div(pp.getRow(i + 1).scale(n.toDouble()).scale(ziMax).add(a))
                        .elementSum() < f[0, i]) {
                    ziMin = ziMax
                    ziMax = ziMax * 2
                }
            }

            for (x in 0..49) {
                val zi = (ziMin + ziMax) / 2
                z[0, i + 1] = zi
                if (pp.getRow(i + 1).scale(zi).div(pp.getRow(i + 1).scale(zi).scale(n.toDouble()).add(a))
                        .elementSum() < f[0, i]) {
                    ziMin = zi
                } else {
                    ziMax = zi
                }
            }
        }
    }

    z.removeCols(setOf(0))
    return z
}
/**
 * Cache xi iter algorithms
 */
@Suppress("unused")
class CacheXiIterAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}