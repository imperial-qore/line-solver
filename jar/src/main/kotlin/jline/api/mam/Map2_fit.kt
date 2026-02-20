/**
 * @file Markovian Arrival Process two-phase fitting algorithms
 * 
 * Fits MAP(2) processes to match specified moments and autocorrelation decay rates.
 * Fundamental algorithm for constructing realistic arrival process models from statistical data.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.io.Ret
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import org.apache.commons.math3.util.FastMath
import java.util.*
import kotlin.math.pow

/**
 * Fits a 2-phase Markovian Arrival Process (MAP2) to match the given moments and decay rate of the autocorrelation function.
 *
 *
 * This method calculates the parameters for a 2-phase MAP based on the first three moments (e1, e2, e3) and the
 * decay rate (g2) of the autocorrelation function.
 *
 * @param e1 the first moment (mean)
 * @param e2 the second moment
 * @param e3 the third moment
 * @param g2 the decay rate of the autocorrelation function
 * @return a Map_fit_return_type object containing the fitted MAP transition matrices and additional fitting information
 */
fun map2_fit(e1: Double, e2: Double, e3: Double, g2: Double): Ret.mamMAPFitReturn {
    var e3 = e3
    val result = Ret.mamMAPFitReturn()
    val r1 = e1
    val r2 = e2 / 2
    val h2 = (r2 - FastMath.pow(r1, 2)) / FastMath.pow(r1, 2)
    if (e3 == -1.0) {
        val scv = (e2 - FastMath.pow(e1, 2)) / FastMath.pow(e1, 2)
        if (1 <= scv && scv < 3) {
            if (g2 < 0) {
                val h3 = h2 - FastMath.pow(h2, 2)
                e3 = 12 * FastMath.pow(e1, 3) * h2 + 6 * FastMath.pow(e1, 3) * h3 + 6 * FastMath.pow(e1,
                    3) * (1 + FastMath.pow(h2, 2))
            } else {
                e3 = 1.501 * FastMath.pow(e2, 2) / e1
            }
        } else if (3 <= scv) {
            e3 = 1.501 * FastMath.pow(e2, 2) / e1
        } else if (0 < scv && scv < 1) {
            e3 = (1 + 1e-10) * (12 * FastMath.pow(e1, 3) * h2 + 6 * FastMath.pow(e1,
                3) * (h2 * (1 - h2 - 2 * FastMath.sqrt(-h2))) + 6 * FastMath.pow(e1, 3) * (1 + FastMath.pow(h2, 2)))
        }
    }

    if (e3 == -2.0) {
        val scv = (e2 - FastMath.pow(e1, 2)) / FastMath.pow(e1, 2)

        if (scv >= 1) {
            e3 = (3.0 / 2 + 1e-6) * FastMath.pow(e2, 2) / e1
        } else if (0 < scv && scv < 1) {
            val h3 = h2 * (1 - h2 - 2 * FastMath.sqrt(-h2))
            e3 = 6 * FastMath.pow(e1, 3) * (h2.pow(2.0) + h3)
        }
    }

    if (e3 == -3.0) {
        val scv = (e2 - FastMath.pow(e1, 2)) / FastMath.pow(e1, 2)

        if (scv >= 1) {
            e3 = FastMath.pow(10.0, 6)
        } else if (0 < scv && scv < 1) {
            val h3 = FastMath.pow(-h2, 2)
            e3 = 6 * FastMath.pow(e1, 3) * (h2.pow(2.0) + h3)
        }
    }

    if (e3 == -4.0) {
        val scv = (e2 - FastMath.pow(e1, 2)) / FastMath.pow(e1, 2)
        val r = Random().nextDouble()

        if (scv >= 1) {
            e3 = r * (3.0 / 2 + 1e-6) * FastMath.pow(e2, 2) / e1 + (1 - r) * FastMath.pow(10.0, 6)
        } else if (0 < scv && scv < 1) {
            val h3 = r * FastMath.pow(-h2, 2) + (1 - r) * h2 * (1 - h2 - 2 * FastMath.sqrt(-h2))
            e3 = 6 * FastMath.pow(e1, 3) * (h2.pow(2.0) + h3)
        }
    }

    if (e3 > -1 && e3 < 0) {
        val scv = (e2 - FastMath.pow(e1, 2)) / FastMath.pow(e1, 2)
        val r = FastMath.abs(e3)

        if (scv >= 1) {
            e3 = r * (3.0 / 2 + 1e-6) * FastMath.pow(e2, 2) / e1 + (1 - r) * FastMath.pow(10.0, 6)
        } else if (0 < scv && scv < 1) {
            val h3 = r * h2 * (1 - h2 - 2 * FastMath.sqrt(-h2)) + (1 - r) * FastMath.pow(-h2, 2)
            e3 = 6 * FastMath.pow(e1, 3) * (h2.pow(2.0) + h3)
        }
    }

    val r3 = e3 / 6
    val h3 = (r3 * r1 - FastMath.pow(r2, 2)) / FastMath.pow(r1, 4)
    val b = h3 + FastMath.pow(r1, 2) - r1
    var c = FastMath.sqrt(b.pow(2.0) + 4 * FastMath.pow(r1, 3))

    if (r1 <= 0) {
        result.error = 10.0
        return result
    }

    if (h2 == 0.0) {
        if (h3 == 0.0 && g2 == 0.0) {
            result.MAP = map_exponential(e1)
        } else {
            result.error = 20.0
            return result
        }
    }

    if (h2 > 0 && h3 > 0) {
        if (b >= 0) {
            if ((b - c) / (b + c) <= g2) {
                val MAP = MatrixCell()
                val D0 = Matrix(2, 2, 4)
                D0[0, 0] = -(2 * h2 + b - c)
                D0[0, 1] = 0
                D0[1, 0] = 0
                D0[1, 1] = -(2 * h2 + b + c)
                D0.scaleEq(1.0 / (2 * r1 * h3))
                MAP[0] = D0

                val D1 = Matrix(2, 2, 4)
                D1[0, 0] = (2 * h2 + b - c) * (1 - b / c + g2 * (1 + b / c))
                D1[0, 1] = (2 * h2 + b - c) * (1 + b / c) * (1 - g2)
                D1[1, 0] = (2 * h2 + b + c) * (1 - b / c) * (1 - g2)
                D1[1, 1] = (2 * h2 + b + c) * (1 + b / c + g2 * (1 - b / c))
                D1.scaleEq(1 / (4 * r1 * h3))
                MAP[1] = D1

                result.MAP = MAP
            } else {
                result.error = 51.0
            }
        } else if (b < 0) {
            if (0 <= g2 && g2 < 1) {
                val MAP = MatrixCell()
                val D0 = Matrix(2, 2, 4)
                D0[0, 0] = -(2 * h2 + b - c)
                D0[0, 1] = 0
                D0[1, 0] = 0
                D0[1, 1] = -(2 * h2 + b + c)
                D0.scaleEq(1 / (2 * r1 * h3))
                MAP[0] = D0

                val D1 = Matrix(2, 2, 4)
                D1[0, 0] = (2 * h2 + b - c) * (1 - b / c + g2 * (1 + b / c))
                D1[0, 1] = (2 * h2 + b - c) * (1 + b / c) * (1 - g2)
                D1[1, 0] = (2 * h2 + b + c) * (1 - b / c) * (1 - g2)
                D1[1, 1] = (2 * h2 + b + c) * (1 + b / c + g2 * (1 - b / c))
                D1.scaleEq(1 / (4 * r1 * h3))
                MAP[1] = D1

                result.MAP = MAP
            } else if (-(h3 + FastMath.pow(h2, 2)) / h2 <= g2 && g2 < 0) {
                val a = (h3 + FastMath.pow(h2, 2)) / h2
                val d1 =
                    ((1 - a) * (2 * h2 * g2 + b - c) + g2 * (b + c) - (b - c)) / ((1 - a) * (2 * h2 + b - c) + 2 * c)
                val d2 = ((g2 - 1) * (b - c)) / ((1 - a) * (2 * h2 + b - c) + 2 * c)
                val MAP = MatrixCell()
                val D0 = Matrix(2, 2, 4)
                D0[0, 0] = -(2 * h2 + b - c)
                D0[0, 1] = (2 * h2 + b - c) * (1 - a)
                D0[1, 0] = 0
                D0[1, 1] = -(2 * h2 + b + c)
                D0.scaleEq(1 / (2 * r1 * h3))
                MAP[0] = D0

                val D1 = Matrix(2, 2, 4)
                D1[0, 0] = (2 * h2 + b - c) * d1
                D1[0, 1] = (2 * h2 + b - c) * (a - d1)
                D1[1, 0] = (2 * h2 + b + c) * d2
                D1[1, 1] = (2 * h2 + b + c) * (1 - d2)
                D1.scaleEq(1 / (2 * r1 * h3))
                MAP[1] = D1

                result.MAP = MAP
            } else {
                result.error = 52.0
            }
        }
        if (result.MAP.size() > 0 && !map_isfeasible(result.MAP)) {
            result.error = -1.0
        }
        return result
    } else if (-0.25 <= h2 && h2 < 0 && h2 * (1 - h2 - 2 * FastMath.sqrt(-h2)) <= h3 && h3 <= -h2.pow(2.0)) {
        if (g2 >= 0) {
            if (g2 <= -(h2 + FastMath.sqrt(-h3)).pow(2.0) / h2) {
                val a = (2 * h2 + b - c) * (h2 + FastMath.sqrt(-h3)) / (2 * h2 * FastMath.sqrt(-h3))
                c = -c
                val d1 =
                    ((1 - a) * (2 * h2 * g2 + b - c) + g2 * (b + c) - (b - c)) / ((1 - a) * (2 * h2 + b - c) + 2 * c)
                val d2 = ((g2 - 1) * (b - c)) / ((1 - a) * (2 * h2 + b - c) + 2 * c)
                val MAP = MatrixCell()
                val D0 = Matrix(2, 2, 4)
                D0[0, 0] = -(2 * h2 + b - c)
                D0[0, 1] = (2 * h2 + b - c) * (1 - a)
                D0[1, 0] = 0
                D0[1, 1] = -(2 * h2 + b + c)
                D0.scaleEq(1 / (2 * r1 * h3))
                MAP[0] = D0

                val D1 = Matrix(2, 2, 4)
                D1[0, 0] = (2 * h2 + b - c) * d1
                D1[0, 1] = (2 * h2 + b - c) * (a - d1)
                D1[1, 0] = (2 * h2 + b + c) * d2
                D1[1, 1] = (2 * h2 + b + c) * (1 - d2)
                D1.scaleEq(1 / (2 * r1 * h3))
                MAP[1] = D1

                result.MAP = MAP
            } else {
                result.error = 53.0
            }
        } else if (g2 < 0) {
            if (g2 >= -(h3 + FastMath.pow(h2, 2)) / h2) {
                val a = (h3 + FastMath.pow(h2, 2)) / h2
                c = -c
                val d1 =
                    ((1 - a) * (2 * h2 * g2 + b - c) + g2 * (b + c) - (b - c)) / ((1 - a) * (2 * h2 + b - c) + 2 * c)
                val d2 = ((g2 - 1) * (b - c)) / ((1 - a) * (2 * h2 + b - c) + 2 * c)
                val MAP = MatrixCell()
                val D0 = Matrix(2, 2, 4)
                D0[0, 0] = -(2 * h2 + b - c)
                D0[0, 1] = (2 * h2 + b - c) * (1 - a)
                D0[1, 0] = 0
                D0[1, 1] = -(2 * h2 + b + c)
                D0.scaleEq(1 / (2 * r1 * h3))
                MAP[0] = D0

                val D1 = Matrix(2, 2, 4)
                D1[0, 0] = (2 * h2 + b - c) * d1
                D1[0, 1] = (2 * h2 + b - c) * (a - d1)
                D1[1, 0] = (2 * h2 + b + c) * d2
                D1[1, 1] = (2 * h2 + b + c) * (1 - d2)
                D1.scaleEq(1 / (2 * r1 * h3))
                MAP[1] = D1

                result.MAP = MAP
            } else {
                result.error = 54.0
            }
        }

        if (result.MAP.size() > 0 && !map_isfeasible(result.MAP)) {
            result.error = -1.0
        }
        return result
    } else {
        if (!(-0.25 <= h2 && h2 < 0 && h2 * (1 - h2 - 2 * FastMath.sqrt(-h2)) <= h3 && h3 <= -h2.pow(2.0))) {
            result.error = 30.0
        } else if ((h2 > 0 && h3 < 0) || h2 * (1 - h2 - 2 * FastMath.sqrt(-h2)) > h3 || h3 <= FastMath.pow(h2, 2)) {
            result.error = 40.0
        } else {
            result.error = 60.0 // Infeasible moment set
        }
    }

    return result
}

/**
 * Fits a 2-phase Markovian Arrival Process (MAP2) using the first three moments.
 *
 *
 * This method provides a simplified fitting interface when `g2` is not available or specified. It uses a default
 * or inferred value for `g2` based on the provided moments.
 *
 * @param e1 the first moment (mean)
 * @param e2 the second moment
 * @param e3 the third moment
 * @return a mapFitReturn object containing the fitted MAP transition matrices and additional fitting information
 */
fun map2_fit(e1: Double, e2: Double, e3: Double): Ret.mamMAPFitReturn {
    return map2_fit(e1, e2, -1.0, e3)
}
/**
 * Map2 Fit algorithms
 */
@Suppress("unused")
class Map2FitAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}