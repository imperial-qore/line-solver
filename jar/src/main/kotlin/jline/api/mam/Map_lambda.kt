package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * MAP arrival rate computation algorithms.
 * 
 * Provides methods for calculating the arrival rate (lambda) of Markovian Arrival Processes (MAP).
 * The arrival rate represents the long-run average rate at which events occur in the MAP and is
 * fundamental to many queueing theory applications.
 */

/**
 * Computes the arrival rate (lambda) of a Markovian Arrival Process (MAP).
 *
 *
 * The arrival rate is calculated using the hidden (D0) and visible (D1) transition matrices of the MAP.
 * This is done by first determining the steady-state vector (pi) of the CTMC underlying the MAP and then
 * using it to compute the rate at which visible transitions (events) occur.
 *
 *
 * @param D0 The hidden transition matrix of the MAP, representing transitions without visible events.
 * @param D1 The visible transition matrix of the MAP, representing transitions with visible events.
 * @return The arrival rate (lambda) of the MAP.
 */

fun map_lambda(D0: Matrix, D1: Matrix): Double {
    val e = Matrix.ones(D0.numRows, 1)
    var lambda = map_piq(D0, D1) // piq
    lambda = lambda.mult(D1)
    lambda = lambda.mult(e)
    return lambda.toDouble()
}

/**
 * Computes the arrival rate (lambda) of a Markovian Arrival Process (MAP) using matrices stored in a MatrixCell.
 *
 *
 * This is a convenience method that extracts the D0 and D1 matrices from a given MAP stored in a MatrixCell and computes the
 * arrival rate (lambda). The computation involves using the steady-state vector (pi) of the CTMC underlying the MAP and
 * determining the rate of visible transitions.
 *
 *
 * @param MAP The Markovian Arrival Process stored in a MatrixCell, containing the (D0, D1) matrices.
 * @return The arrival rate (lambda) of the MAP.
 */

fun map_lambda(MAP: MatrixCell): Double {
    return map_lambda(MAP[0], MAP[1])
}

/**
 * MAP arrival rate computation algorithms.
 * 
 * Provides methods for calculating the arrival rate (lambda) of Markovian Arrival Processes (MAP).
 * The arrival rate represents the long-run average rate at which events occur in the MAP and is
 * fundamental to many queueing theory applications.
 */
@Suppress("unused")
class MapLambda {
    companion object {
        // Class documentation marker for Dokka
    }
}
