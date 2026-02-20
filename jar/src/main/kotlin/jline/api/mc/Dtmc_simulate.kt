/**
 * @file Discrete-time Markov chain Monte Carlo simulation
 * 
 * Generates sample trajectories for DTMCs by sampling from the transition probability
 * matrix at each step. Includes optimizations for absorbing states and provides
 * foundation for statistical analysis and validation of DTMC models.
 *
 * @since LINE 3.0
 */
package jline.api.mc

import jline.util.matrix.Matrix
import kotlin.random.Random

/**
 * Simulate a discrete-time Markov chain trajectory
 *
 * @param P Transition matrix of the DTMC
 * @param pi0 Initial state distribution
 * @param n Number of simulation steps
 * @return Array containing the state trajectory
 */
fun dtmc_simulate(P: Matrix, pi0: Matrix, n: Int): IntArray {
    val numStates = P.numRows
    val sts = IntArray(n)
    
    // Sample initial state
    val rnd0 = Random.nextDouble()
    val cpi0 = DoubleArray(numStates)
    var cumSum = 0.0
    for (i in 0 until numStates) {
        cumSum += pi0[i]
        cpi0[i] = cumSum
    }
    
    var st = 0
    for (i in 0 until numStates) {
        if (rnd0 <= cpi0[i] && pi0[i] > 0) {
            st = i
            break
        }
    }
    
    // Precompute cumulative probabilities for all states
    val F = Matrix.zeros(numStates, numStates)
    for (i in 0 until numStates) {
        var cumulative = 0.0
        for (j in 0 until numStates) {
            cumulative += P[i, j]
            F[i, j] = cumulative
        }
    }
    
    // Simulate trajectory
    for (step in 0 until n) {
        sts[step] = st
        
        // Check for absorbing state
        if (F[st, numStates - 1] == 0.0 || P[st, st] == 1.0) {
            // Fill remaining steps with same state
            for (k in step + 1 until n) {
                sts[k] = st
            }
            break
        }
        
        // Sample next state
        val rnd = Random.nextDouble()
        var nextState = 0
        for (j in 0 until numStates) {
            if (rnd <= F[st, j] && P[st, j] > 0) {
                nextState = j
                break
            }
        }
        st = nextState
    }
    
    return sts
}
/**
 * DTMC simulate algorithms
 */
@Suppress("unused")
class DtmcSimulateAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}