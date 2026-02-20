/**
 * @file Continuous-time Markov chain Monte Carlo simulation
 * 
 * Generates sample paths for CTMCs using the standard simulation algorithm with
 * exponentially distributed sojourn times and discrete transition sampling. Essential
 * for validation of analytical results and studying rare events in complex systems.
 *
 * @since LINE 3.0
 */
package jline.api.mc

import jline.io.Ret
import jline.lang.processes.Exp
import jline.util.matrix.Matrix
import java.util.*

/**
 * Form a random infinitesimal generator of CTMC
 *
 * @param Q   infinitesimal generator of CTMC
 * @param pi0 initial state distribution vector
 * @param n   times of simulations
 * @return Infinitesimal generator of CTMC
 */

fun ctmc_simulate(Q: Matrix, pi0: DoubleArray?, n: Int): Ret.ctmcSimulation {
    return ctmc_simulate(Q, pi0, n, Random())
}

fun ctmc_simulate(Q: Matrix, pi0: DoubleArray?, n: Int, random: Random): Ret.ctmcSimulation {
    var pi0 = pi0
    val numStates = Q.length()
    if (pi0 == null || pi0.size == 0) {
        pi0 = DoubleArray(numStates)
        var sum = 0.0

        for (i in 0..<numStates) {
            pi0[i] = random.nextDouble()
            sum += pi0[i]
        }

        for (i in 0..<numStates) {
            pi0[i] /= sum
        }
    }
    var cumulative = 0.0
    var r = random.nextDouble()
    var st = 0
    for (i in pi0.indices) {
        cumulative += pi0[i]
        if (r < cumulative) {
            st = i
            break
        }
    }

    val F = Matrix(numStates, numStates)
    for (i in 0..<numStates) {
        var rowSum = 0.0
        for (j in 0..<numStates) {
            if (i != j) {
                rowSum += Q[i, j]
                F[i, j] = rowSum
            }
        }
        for (j in 0..<numStates) {
            F[i, j] = F[i, j] / rowSum
        }
    }

    val result = Ret.ctmcSimulation()
    result.states = IntArray(n)
    result.sojournTimes = DoubleArray(n)

    for (i in 0..<n) {
        result.states[i] = st
        val expDist = Exp(-Q[st, st])
        result.sojournTimes[i] = expDist.sample(1, random)[0]

        r = random.nextDouble()
        for (j in 0..<numStates) {
            if (r < F[st, j]) {
                st = j
                break
            }
        }
    }
    return result
}
/**
 * CTMC simulate algorithms
 */
@Suppress("unused")
class CtmcSimulateAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}