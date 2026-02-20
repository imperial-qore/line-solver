/**
 * @file Absorbing Phase-type distribution assembly from parameters
 * 
 * Constructs APH(2) transition matrices from specified rates and transition probabilities.
 * Used for building phase-type distributions from fitted or prescribed parameters.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Assembles an acyclic phase-type (APH) distribution with two phases (APH(2)) using the given parameters.
 *
 *
 * This method constructs the transition matrices `D0` and `D1` for an APH(2) distribution based on the parameters:
 * `l1` and `l2` (the rates of the exponential phases), and `p1` (the probability of transitioning between the phases).
 * The resulting distribution has two phases, each represented by an exponential distribution, with transitions
 * between the phases as specified by the given parameters.
 *
 * @param l1 the rate of the first exponential phase
 * @param l2 the rate of the second exponential phase
 * @param p1 the probability of transitioning from the first phase to the second phase
 * @return a MatrixCell containing the transition matrices `D0` and `D1` of the APH(2) distribution
 */
fun aph2_assemble(l1: Double, l2: Double, p1: Double): MatrixCell {
    val APH = MatrixCell()
    val D0 = Matrix(2, 2, 4)
    D0[0, 0] = -1 / l1
    D0[0, 1] = 1 / l1 * p1
    D0[1, 0] = 0
    D0[1, 1] = -1 / l2

    val D1 = Matrix(2, 2, 4)
    D1[0, 0] = 1 / l1 * (1 - p1)
    D1[0, 1] = 0
    D1[1, 0] = 1 / l2
    D1[1, 1] = 0

    APH[0] = D0
    APH[1] = D1

    return APH
}
/**
 * APH 2 assemble algorithms
 */
@Suppress("unused")
class Aph2AssembleAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}