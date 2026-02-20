/**
 * @file Absorbing Phase-type distribution general fitting algorithms
 * 
 * Fits APH distributions to specified moments using optimization and approximation techniques.
 * Fundamental algorithms for phase-type distribution parameter estimation and model construction.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.GlobalConstants.Inf

import jline.util.Utils
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import org.apache.commons.math3.util.FastMath
import kotlin.math.pow

/**
 * Fits an acyclic phase-type (APH) distribution to the given moments of a random variable.
 *
 *
 * This method approximates a random variable with an APH distribution, characterized by the first three moments:
 * the mean (e1), the second moment (e2), and the third moment (e3). The APH distribution is represented as a MAP
 * with two matrices: D0 and D1. If the second or third moment is infinite, the method defaults to fitting an exponential MAP.
 *
 *
 * The method iteratively tries to find an appropriate order `n` for the APH distribution up to a maximum order `nmax`.
 * The order determines the number of phases in the distribution. If the moments cannot be matched exactly within the
 * given `nmax`, an approximate solution is returned.
 *
 * @param e1   the first moment (mean) of the distribution
 * @param e2   the second moment of the distribution
 * @param e3   the third moment of the distribution
 * @param nmax the maximum number of phases to consider in the APH distribution
 * @return a MatrixCell containing the transition matrices `D0` and `D1` of the fitted APH distribution
 */
fun aph_fit(e1: Double, e2: Double, e3: Double, nmax: Int): MatrixCell {
    var APH = MatrixCell()
    if (Utils.isInf(e2) || Utils.isInf(e3)) {
        return map_exponential(e1)
    }

    var n2 = e2 / FastMath.pow(e1, 2)
    var n3 = e3 / e1 / e2

    var n2_feas = false
    var n3_ubfeas = false
    var n3_lbfeas = false
    var n = 1.0
    var un = 0.0
    val un_1 = un
    while ((!n2_feas || !n3_lbfeas || !n3_ubfeas) && n < nmax) {
        n++
        val pn =
            ((n + 1) * (n2 - 2) / (3 * n2 * (n - 1))) * (-2 * FastMath.sqrt(n + 1) / FastMath.sqrt(4 * (n + 1) - 3 * n * n2) - 1)
        val an = (n2 - 2) / (pn * (1 - n2) + FastMath.sqrt(pn.pow(2.0) + pn * n * (n2 - 2) / (n - 1)))
        val ln =
            ((3 + an) * (n - 1) + 2 * an) / ((n - 1) * (1 + an * pn)) - (2 * an * (n + 1)) / (2 * (n - 1) + an * pn * (n * an + 2 * n - 2))

        un =
            (1 / (n.pow(2.0) * n2)) * (2 * (n - 2) * (n * n2 - n - 1) * FastMath.sqrt(1 + n * (n2 - 2) / (n - 1)) + (n + 2) * (3 * n * n2 - 2 * n - 2))

        if (n2 >= (n + 1) / n && n2 <= (n + 4) / (n + 1)) {
            n2_feas = true
            if (n3 >= ln) {
                n3_lbfeas = true
            }
        } else if (n2 >= (n + 4) / (n + 1)) {
            n2_feas = true
            if (n3 >= n2 * (n + 1) / n) {
                n3_lbfeas = true
            }
        }

        if (n2 >= (n + 1) / n && n2 <= n / (n - 1)) {
            n2_feas = true
            if (n3 <= un) {
                n3_ubfeas = true
            }
        } else if (n2 >= n / (n - 1)) {
            n2_feas = true
            if (n3 < Inf) {
                n3_ubfeas = true
            }
        }
    }

    if ((!n2_feas || !n3_lbfeas || !n3_ubfeas) || (n == nmax.toDouble())) {
        print("'cannot match moment set exactly'")
        n2 = (n + 1) / n
        n3 = 2 * n2 - 1
    }

    if (n2 <= n / (n - 1) || n3 <= 2 * n2 - 1) {
        val b =
            2 * (4 - n * (3 * n2 - 4)) / (n2 * (4 + n - n * n3) + FastMath.sqrt(n * n2) * FastMath.sqrt(12 * FastMath.pow(
                n2,
                2) * (n + 1) + 16 * n3 * (n + 1) + n2 * (n * (n3 - 15) * (n3 + 1) - 8 * (n3 + 3))))
        val a = (b * n2 - 2) * (n - 1) * b / ((b - 1) * n)
        val p = (b - 1) / a

        val lambda = 1.0
        val mu = lambda * (n - 1) / a
        val alpha = Matrix(1, n.toInt(), n.toInt())
        alpha[0, 0] = p
        alpha[0, n.toInt() - 1] = 1 - p
        val T = Matrix(n.toInt(), n.toInt(), FastMath.pow(n, 2).toInt())
        run {
            var i = 0
            while (i < n - 1) {
                T[i, i] = -mu
                T[i, i + 1] = mu

                i++
            }
        }
        T[n.toInt() - 1, n.toInt() - 1] = -lambda
        val neg_T = T.copy()
        neg_T.scaleEq(-1.0)
        val one = Matrix(1, n.toInt() - 1, n.toInt() - 1)
        var i = 0
        while (i < n - 1) {
            one[0, i] = 1
            i++
        }
        APH = map_scale(map_normalize(T, neg_T.mult(one).mult(alpha))[0],
            map_normalize(T, neg_T.mult(one).mult(alpha))[1],
            e1)
    } else if (n2 > n / (n - 1) && n3 > un_1) {
        val K1 = n - 1
        val K2 = n - 2
        val K3 = 3 * n2 - 2 * n3
        val K4 = n3 - 3
        val K5 = n - n2
        val K6 = 1 + n2 - n3
        val K7 = n + n2 - n * n2
        val K8 = 3 + 3 * FastMath.pow(n2, 2) + n3 - 3 * n2 * n3
        val K9 = 108 * FastMath.pow(K1, 2) * (4 * FastMath.pow(K2, 2) * K3 * FastMath.pow(n, 2) * n2 + FastMath.pow(K1,
            2) * K2 * FastMath.pow(K4, 2) * n * FastMath.pow(n2,
            2) + 4 * K1 * K5 * (K5.pow(2.0) - 3 * K2 * K6 * n * n2) + FastMath.sqrt(-16 * FastMath.pow(K1,
            2) * FastMath.pow(K7, 6) + FastMath.pow(4 * K1 * FastMath.pow(K5, 3) + FastMath.pow(K1,
            2) * K2 * FastMath.pow(K4, 2) * n * FastMath.pow(n2, 2) + 4 * K2 * n * n2 * (K4 * FastMath.pow(n,
            2) - 3 * K6 * n2 + K8 * n), 2)))
        val K10 = FastMath.pow(K4, 2) / (4 * FastMath.pow(K3, 2)) - K5 / (K1 * K3 * n2)
        val K11 = FastMath.pow(2.0, 1.0 / 3) * (3 * FastMath.pow(K5,
            2) + K2 * (K3 + 2 * K4) * n * n2) / (K3 * FastMath.pow(K9, 1.0 / 3) * n2)
        val K12 = FastMath.pow(K9, 1.0 / 3) / (3 * FastMath.pow(2.0, 7.0 / 3) * FastMath.pow(K1, 2) * K3 * n2)
        val K13 = FastMath.sqrt(K10 + K11 + K12)
        val K14 = (6 * K1 * K3 * K4 * K5 + 4 * K2 * FastMath.pow(K3, 2) * n - FastMath.pow(K1, 2) * FastMath.pow(K4,
            2) * n2) / (4 * FastMath.pow(K1, 2) * FastMath.pow(K3, 3) * K13 * n2)
        val K15 = -K4 / (2 * K3)
        val K16 = FastMath.sqrt(2 * K10 - K11 - K12 - K14)
        val K17 = FastMath.sqrt(2 * K10 - K11 - K12 + K14)
        val K18 = 36 * FastMath.pow(K5, 3) + 36 * K2 * K4 * K5 * n * n2 + 9 * K1 * K2 * FastMath.pow(K4,
            2) * n * FastMath.pow(n2, 2) - FastMath.sqrt(81 * FastMath.pow(4 * FastMath.pow(K5,
            2) + 4 * K2 * K4 * K5 * n * n2 + K1 * K2 * FastMath.pow(K4, 2) * n * FastMath.pow(n2, 2),
            2) - 48 * FastMath.pow(3 * FastMath.pow(K5, 2) + 2 * K2 * K4 * n * n2, 3))
        val K19 = -K5 / (K1 * K4 * n2) - FastMath.pow(2.0, 2.0 / 3) * (3 * FastMath.pow(K5,
            2) + 2 * K2 * K4 * n * n2) / (3.0.pow(1.0 / 3) * K1 * K4 * n2 * FastMath.pow(K18, 1.0 / 3)) - FastMath.pow(
            K18,
            1.0 / 3) / (6.0.pow(2.0 / 3) * K1 * K4 * n2)
        val K20 =
            6 * K1 * K3 * K4 * K5 + 4 * K2 * FastMath.pow(K3, 2) * n - FastMath.pow(K1, 2) * FastMath.pow(K4, 2) * n2
        val K21 = K11 + K12 + K5 / (2 * n * K1 * K3)
        val K22 =
            FastMath.sqrt(3 * FastMath.pow(K4, 2) / (4 * FastMath.pow(K3, 2)) - 3 * K5 / (K1 * K3 * n2) + FastMath.sqrt(
                4 * FastMath.pow(K21, 2) - n * K2 / (n2 * FastMath.pow(K1, 2) * K3)))
        var f = 0.0
        if (n3 > un_1 && n3 < 3 * n2 / 2) {
            f = K13 + K15 - K17
        } else if (n3 == 2 * n2 / 2) {
            f = K19
        } else if (n3 > 3 * n2 / 2 && K20 > 0) {
            f = -K13 + K15 + K16
        } else if (K20 == 0.0) {
            f = K15 + K22
        } else if (K20 < 0) {
            f = K13 + K15 + K17
        }
        val a = 2 * (f - 1) * (n - 1) / ((n - 1) * (n2 * FastMath.pow(f, 2) - 2 * f + 2) - n)
        val p = (f - 1) * a
        val lambda = 1.0
        val mu = lambda * (n - 1) / a
        val alpha = Matrix(1, n.toInt(), n.toInt())
        alpha[0, 0] = p
        alpha[0, 1] = 1 - p
        val T = Matrix(n.toInt(), n.toInt(), FastMath.pow(n, 2).toInt())
        run {
            var i = 0
            while (i < n - 1) {
                T[i, i] = -mu
                T[i, i + 1] = mu
                i++
            }
        }
        T[n.toInt() - 1, n.toInt() - 1] = -mu
        T[0, 0] = -lambda
        T[0, 1] = lambda
        val neg_T = T.copy()
        neg_T.scaleEq(-1.0)
        val one = Matrix(1, n.toInt(), n.toInt())
        var i = 0
        while (i < n) {
            one[0, i] = 1
            i++
        }
        APH = map_scale(map_normalize(T, neg_T.mult(one).mult(alpha))[0],
            map_normalize(T, neg_T.mult(one).mult(alpha))[1],
            e1)
    } else {
        print("moment set cannot be matched with an APH distribution")
    }

    return APH
}

/**
 * Fits an acyclic phase-type (APH) distribution to the given moments of a random variable.
 *
 *
 * This method approximates a random variable with an APH distribution, characterized by the first three moments:
 * the mean (e1), the second moment (e2), and the third moment (e3). The APH distribution is represented as a MAP
 * with two matrices: D0 and D1. If the second or third moment is infinite, the method defaults to fitting an exponential MAP.
 *
 *
 * This version of the method uses a default maximum order of 10 for the APH distribution.
 *
 * @param e1 the first moment (mean) of the distribution
 * @param e2 the second moment of the distribution
 * @param e3 the third moment of the distribution
 * @return a MatrixCell containing the transition matrices `D0` and `D1` of the fitted APH distribution
 */
fun aph_fit(e1: Double, e2: Double, e3: Double): MatrixCell {
    return aph_fit(e1, e2, e3, 10)
}
/**
 * APH  fit algorithms
 */
@Suppress("unused")
class AphFitAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}