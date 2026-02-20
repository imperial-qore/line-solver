/**
 * @file Absorbing Phase-type distribution comprehensive fitting
 * 
 * Fits multiple APH(2) distributions to match given moments with exhaustive parameter search.
 * Provides comprehensive fitting solutions for phase-type distribution approximation.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import org.apache.commons.math3.util.FastMath
import kotlin.math.min

/**
 * Fits a set of acyclic phase-type (APH) distributions with two phases (APH(2)) to match the given moments of a random variable.
 *
 *
 * This method approximates a random variable with multiple APH(2) distributions using the first three moments:
 * the mean (M1), the second moment (M2), and the third moment (M3). The method checks the feasibility of these moments
 * and computes the parameters for possible APH(2) distributions. If the moments cannot be exactly matched, an approximate
 * APH(2) distribution is returned.
 *
 *
 * The method first calculates the squared coefficient of variation (SCV) and the lower bound for M3 (M3lb) based on the moments.
 * It then determines the number of solutions (n) and computes the parameters (h1, h2, r1) for the APH(2) distributions.
 * Valid distributions are stored in the map `APHS`, with each entry keyed by an integer and containing a MatrixCell
 * representing the transition matrices `D0` and `D1`.
 *
 * @param M1 the first moment (mean) of the distribution
 * @param M2 the second moment of the distribution
 * @param M3 the third moment of the distribution
 * @return a map containing possible APH(2) distributions, each represented by a MatrixCell
 */
fun aph2_fitall(M1: Double, M2: Double, M3: Double): Map<Int, MatrixCell> {
    val APHS: MutableMap<Int, MatrixCell> = HashMap()
    val degentol = 1e-8
    val SCV = (M2 - FastMath.pow(M1, 2)) / FastMath.pow(M1, 2)
    val M3lb = 3 * FastMath.pow(M1, 3) * (3 * SCV - 1 + FastMath.sqrt(2.0) * FastMath.pow(1 - SCV, 1.5))
    val tmp0: Double
    if (SCV <= 1 && FastMath.abs(M3 - M3lb) < degentol) {
        tmp0 = 0.0
    } else {
        tmp0 = FastMath.pow(M3, 2.0 / 9) + ((8 * FastMath.pow(M1, 3)) / 3.0 - 2 * M2 * M1) * M3 - 3 * FastMath.pow(M1,
            2) * FastMath.pow(M2, 2) + 2 * FastMath.pow(M2, 3)
        if (tmp0 < 0) {
            APHS[0] = aph_fit(M1, M2, M3, 2)
            return APHS
        }
    }

    val tmp1 = 3 * FastMath.sqrt(tmp0)
    val tmp2 = M3 - 3 * M1 * M2
    val tmp3 = (6 * M2 - 12 * FastMath.pow(M1, 2))
    val n = if (tmp0 == 0.0) {
        1
    } else {
        2
    }

    val h1v = Matrix(n, 1, n)
    val h2v = Matrix(n, 1, n)
    val r1v = Matrix(n, 1, n)

    if (n == 1) {
        h2v[0, 0] = tmp2 / tmp3
        h1v[0, 0] = tmp2 / tmp3
    } else {
        h2v[0, 0] = (tmp2 + tmp1) / tmp3
        h2v[1, 0] = (tmp2 - tmp1) / tmp3
        h1v[1, 0] = h2v[0, 0]
        h1v[0, 0] = h2v[1, 0]
    }

    for (j in 0..<n) {
        val h1 = h1v[j]
        val h2 = h2v[j]
        r1v[j] = (M1 - h1) / h2
    }

    var idx = 0
    for (j in 0..<n) {
        val h1 = h1v[j]
        val h2 = h2v[j]
        var r1 = r1v[j]
        if (h1 > 0 && h2 > 0 && r1 >= -degentol && r1 <= (1 + degentol)) {
            r1 = FastMath.max(min(r1, 1.0), 0.0)
            APHS[idx] = aph2_assemble(h1, h2, r1)
            idx++
        }
    }

    if (APHS.isEmpty()) {
        APHS[0] = aph_fit(M1, M2, M3, 2)
    }

    return APHS
}
/**
 * APH 2 fitall algorithms
 */
@Suppress("unused")
class Aph2Fitall {
    companion object {
        // Class documentation marker for Dokka
    }
}