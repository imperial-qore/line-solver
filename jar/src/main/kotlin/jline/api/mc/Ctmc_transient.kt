/**
 * @file Continuous-time Markov chain transient analysis
 * 
 * Computes transient probabilities for CTMCs using numerical integration of the 
 * Chapman-Kolmogorov differential equation π'(t) = π(t)·Q. Employs LSODA solver 
 * for efficient adaptive time-stepping in stiff and non-stiff problems.
 * 
 * @since LINE 3.0
 */
package jline.api.mc

import jline.util.Pair
import jline.util.matrix.Matrix
import odesolver.LSODA
import org.apache.commons.math3.exception.DimensionMismatchException
import org.apache.commons.math3.exception.MaxCountExceededException
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations
import java.util.*
import java.util.stream.Collectors

/**
 * Return the transient probability distribution of the CTMC
 *
 * @param Q   infinitesimal generator of CTMC
 * @param pi0 initial state distribution vector
 * @param t0  transient analysis period start timepoint
 * @param t1  transient analysis period start boundary [t0,t1]
 * @return transient analysis results including time points t and transient probability pi
 */

fun ctmc_transient(Q: Matrix, pi0: Matrix, t0: Double, t1: Double): Pair<DoubleArray, List<DoubleArray>> {
    return ctmc_transient(Q, pi0, t0, t1, null)
}

/**
 * Return the transient probability distribution of the CTMC with optional fixed timestep
 *
 * @param Q        infinitesimal generator of CTMC
 * @param pi0      initial state distribution vector
 * @param t0       transient analysis period start timepoint
 * @param t1       transient analysis period start boundary [t0,t1]
 * @param timestep fixed timestep for equally-spaced time points (null for adaptive stepping)
 * @return transient analysis results including time points t and transient probability pi
 */
fun ctmc_transient(
    Q: Matrix,
    pi0: Matrix,
    t0: Double,
    t1: Double,
    timestep: Double?
): Pair<DoubleArray, List<DoubleArray>> {
    val lsoda = LSODA(0.0, 0.0, 1.0e-6, 1.0e-6, 12, 5)
    val ode: FirstOrderDifferentialEquations = object : FirstOrderDifferentialEquations {
        override fun getDimension(): Int {
            return pi0.length()
        }

        @Throws(MaxCountExceededException::class, DimensionMismatchException::class)
        override fun computeDerivatives(t: Double, y: DoubleArray, ydot: DoubleArray) {
            for (i in 0..<Q.numCols) {
                for (j in 0..<Q.numRows) ydot[i] += y[j] * Q[j, i]
            }
        }
    }
    val y = DoubleArray(pi0.length())

    // Use fixed timestep if specified
    if (timestep != null && timestep > 0) {
        val timePoints = mutableListOf<Double>()
        var t = t0
        while (t <= t1) {
            timePoints.add(t)
            t += timestep
        }
        if (timePoints.last() != t1) {
            timePoints.add(t1)
        }

        val piResults = mutableListOf<DoubleArray>()
        timePoints.add(0, t0) // Ensure t0 is first point
        for (i in 1 until timePoints.size) {
            lsoda.integrate(ode, timePoints[i - 1], if (i == 1) pi0.toArray1D() else y, timePoints[i], y)
            piResults.add(y.clone())
        }

        return Pair(timePoints.drop(1).toDoubleArray(), piResults)
    } else {
        // Original adaptive stepping behavior
        lsoda.integrate(ode, t0, pi0.toArray1D(), t1, y)
        val pi = lsoda.yvec.stream().map { doubleArray: Array<Double>? ->
            Arrays.stream(doubleArray).mapToDouble { obj: Double -> obj.toDouble() }.toArray()
        }.collect(Collectors.toList())
        val t = lsoda.tvec.stream().mapToDouble { obj: Double -> obj.toDouble() }.toArray()
        return Pair(t, pi)
    }
}

/**
 * Return the transient probability distribution of the CTMC
 *
 * @param Q   infinitesimal generator of CTMC
 * @param pi0 initial state distribution vector
 * @param t1  transient analysis period start boundary [0,t1]
 * @return transient analysis results including time points t and transient probability pi
 */

fun ctmc_transient(Q: Matrix, pi0: Matrix, t1: Double): Pair<DoubleArray, List<DoubleArray>> {
    return ctmc_transient(Q, pi0, 0.0, t1)
}

/**
 * Return the transient probability distribution of the CTMC
 *
 * @param Q  infinitesimal generator of CTMC
 * @param t1 initial state distribution vector
 * @return transient analysis results including time points t and transient probability pi
 */

fun ctmc_transient(Q: Matrix, t1: Double): Pair<DoubleArray, List<DoubleArray>> {
    val pi0 = DoubleArray(Q.length())
    Arrays.fill(pi0, 1.0 / Q.length())
    return ctmc_transient(Q, Matrix(pi0), 0.0, t1)
}
/**
 * CTMC transient algorithms
 */
@Suppress("unused")
class CtmcTransientAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}