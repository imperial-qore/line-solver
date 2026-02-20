/**
 * @file Limited load-dependent function evaluation with spline interpolation
 * 
 * Evaluates limited load-dependent (LLD) scaling functions using spline interpolation for
 * continuous queue-length values. Supports multi-server stations with load-dependent service
 * rates and provides smooth interpolation between discrete scaling points.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.ld

import jline.util.Maths
import jline.util.Utils
import jline.util.matrix.Matrix
import org.apache.commons.math3.analysis.interpolation.SplineInterpolator
import org.apache.commons.math3.util.FastMath
import java.lang.Double

/**
 * Evaluate limited-load dependent (LLD) function
 *
 * @param n          Queue-length values. The values can be continuous.
 * @param lldscaling If not null, then the LLD function uses lldscaling to interpolate continuous values of n
 * @param nservers   If not null, then the LLD function is set to be for a multi-server with nserver stations
 * @return Interpolated LLD function values
 */

fun pfqn_lldfun(n: Matrix, lldscaling: Matrix, nservers: Matrix?): Matrix {
    val M = n.length()
    val r = Matrix(M, 1)
    r.fill(1.0)
    val smax = lldscaling.numCols
    val alpha = 20.0

    for (i in 0..<M) {
        if (!(nservers == null || nservers.isEmpty)) {
            if (Utils.isInf(nservers[i])) {
                r[i, 0] = 1
            } else {
                val softminValue = r[i, 0] / Maths.softmin(n[i], nservers[i], alpha)
                if (Double.isNaN(softminValue)) r[i, 0] = 1.0 / FastMath.min(n[i], nservers[i])
                else r[i, 0] = softminValue
            }
        }

        if (!lldscaling.isEmpty) {
            val lldscaling_i = Matrix(1, smax)
            Matrix.extract(lldscaling, i, i + 1, 0, smax, lldscaling_i, 0, 0)
            if (lldscaling_i.elementMax() != lldscaling_i.elementMin()) {
                val X = DoubleArray(smax)
                val V = DoubleArray(smax)
                for (j in 0..<smax) {
                    X[j] = (j + 1).toDouble()
                    V[j] = lldscaling[i, j]
                }
                r[i, 0] = r[i, 0] / (SplineInterpolator().interpolate(X, V)).value(n[i])
            }
        }
    }
    return r
}
/**
 * PFQN lldfun algorithms
 */
@Suppress("unused")
class PfqnLldfunAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}