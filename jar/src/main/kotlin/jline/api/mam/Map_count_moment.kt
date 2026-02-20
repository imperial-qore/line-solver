/**
 * @file Markovian Arrival Process counting process moment analysis
 * 
 * Computes power moments of MAP counting processes using moment generating functions and
 * numerical differentiation. Essential for advanced statistical analysis of arrival patterns.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.Maths
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Computes power moments of counts at resolution t for a Markovian Arrival Process (MAP).
 *
 * This function calculates the power moments of the counting process N(t) for a MAP
 * at a given time resolution t. The moments are computed using the moment generating
 * function approach with numerical differentiation.
 *
 * @param MAP The Markovian Arrival Process stored in a MatrixCell, containing the (D0, D1) matrices
 * @param t The time resolution at which to compute the moments
 * @param order The order of the moment to compute (e.g., 1 for mean, 2 for second moment)
 * @return The power moment of the specified order
 */
fun map_count_moment(MAP: MatrixCell, t: Double, order: Int): Double {
    val n = MAP[0].getNumRows()
    val theta = map_prob(MAP)
    val e = Matrix.ones(n, 1)
    
    // Use numerical differentiation for moment computation
    // We compute the moment generating function M(z) = theta * exp((D0 + D1*exp(z))*t) * e
    // The k-th moment is the k-th derivative of M(z) evaluated at z=0
    
    return computeHigherOrderMoment(MAP, t, theta, e, order)
}

/**
 * Computes multiple power moments of counts at resolution t for a MAP.
 *
 * @param MAP The Markovian Arrival Process stored in a MatrixCell
 * @param t The time resolution at which to compute the moments
 * @param orders Array of moment orders to compute
 * @return Array of power moments corresponding to the specified orders
 */
fun map_count_moment(MAP: MatrixCell, t: Double, orders: IntArray): DoubleArray {
    val results = DoubleArray(orders.size)
    for (i in orders.indices) {
        results[i] = map_count_moment(MAP, t, orders[i])
    }
    return results
}

/**
 * Moment generating function M(z) = theta * exp((D0 + D1*exp(z))*t) * e
 */
private fun mgfunc(z: Double, MAP: MatrixCell, t: Double, theta: Matrix, e: Matrix): Double {
    val D0 = MAP[0]
    val D1 = MAP[1]
    val n = D0.getNumRows()
    
    // Compute D0 + D1*exp(z)
    val expZ = kotlin.math.exp(z)
    val D1_scaled = D1.scale(expZ)
    val generator = D0.add(1.0, D1_scaled)
    
    // Compute exp((D0 + D1*exp(z))*t)
    val generatorScaled = generator.scale(t)
    val expMatrix = Maths.matrixExp(generatorScaled)
    
    // Compute theta * exp(...) * e
    val temp = theta.mult(expMatrix)
    return temp.mult(e).get(0, 0)
}

/**
 * Computes higher-order moments using finite differences.
 */
private fun computeHigherOrderMoment(MAP: MatrixCell, t: Double, theta: Matrix, e: Matrix, order: Int): Double {
    val h = 1e-4
    val points = order + 1
    val weights = computeFiniteDifferenceWeights(order, points)
    
    var result = 0.0
    val offset = points / 2
    
    for (i in 0 until points) {
        val z = (i - offset) * h
        val mgf_val = mgfunc(z, MAP, t, theta, e)
        result += weights[i] * mgf_val
    }
    
    return result / Math.pow(h, order.toDouble())
}

/**
 * Computes finite difference weights for numerical differentiation.
 */
private fun computeFiniteDifferenceWeights(order: Int, points: Int): DoubleArray {
    // Simple central difference weights for common cases
    return when (order) {
        1 -> when (points) {
            3 -> doubleArrayOf(-0.5, 0.0, 0.5)
            5 -> doubleArrayOf(1.0/12.0, -2.0/3.0, 0.0, 2.0/3.0, -1.0/12.0)
            else -> doubleArrayOf(-0.5, 0.0, 0.5)
        }
        2 -> when (points) {
            3 -> doubleArrayOf(1.0, -2.0, 1.0)
            5 -> doubleArrayOf(-1.0/12.0, 4.0/3.0, -5.0/2.0, 4.0/3.0, -1.0/12.0)
            else -> doubleArrayOf(1.0, -2.0, 1.0)
        }
        3 -> doubleArrayOf(-0.5, 1.0, 0.0, -1.0, 0.5)
        4 -> doubleArrayOf(1.0, -4.0, 6.0, -4.0, 1.0)
        else -> {
            // Fallback for higher orders - use simple central difference
            val weights = DoubleArray(points)
            val center = points / 2
            for (i in 0 until points) {
                weights[i] = if (i == center) -2.0 else if (kotlin.math.abs(i - center) == 1) 1.0 else 0.0
            }
            weights
        }
    }
}
/**
 * MAP count moment algorithms
 */
@Suppress("unused")
class MapCountMoment {
    companion object {
        // Class documentation marker for Dokka
    }
}