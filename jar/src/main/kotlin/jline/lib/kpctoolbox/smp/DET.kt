package jline.lib.kpctoolbox.smp

import jline.api.mam.map_pie
import jline.lib.kpctoolbox.basic.e
import jline.lib.kpctoolbox.mc.dtmc_solve
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import org.apache.commons.math3.linear.LUDecomposition
import org.apache.commons.math3.linear.MatrixUtils
import org.apache.commons.math3.util.FastMath
import java.util.Random

/**
 * Deterministic (Semi-Markov) Process functions.
 *
 * Ported from MATLAB: matlab/lib/kpctoolbox/smp/det/
 */

/**
 * Computes the embedded DTMC of a deterministic MAP.
 *
 * @param DET Deterministic MAP as {D0, D1}
 * @return Embedded discrete-time transition probability matrix
 */
fun det_embedded(DET: MatrixCell): Matrix {
    val D0 = DET[0]
    val D1 = DET[1]
    val n = D0.numRows

    // P = inv(-D0) * D1 is the embedded DTMC
    // Since DET is deterministic, D0 is diagonal
    val P = Matrix(n, n)

    for (i in 0 until n) {
        val rate = -D0.get(i, i)
        if (rate > 0) {
            for (j in 0 until n) {
                P.set(i, j, D1.get(i, j) / rate)
            }
        } else {
            // Absorbing state fallback
            P.set(i, i, 1.0)
        }
    }

    return P
}

/**
 * Computes the k-th moment of a deterministic process.
 *
 * @param DET Deterministic MAP as {D0, D1}
 * @param kset Array of moment orders to compute
 * @return Array of moment values
 */
fun det_moment(DET: MatrixCell, kset: IntArray): DoubleArray {
    val D0 = DET[0]
    val n = D0.numRows

    // Get embedded DTMC
    val P = det_embedded(DET)

    // Get stationary distribution
    val pi = dtmc_solve(P)

    // Compute inv(-D0)
    val invD0 = Matrix(n, n)
    for (i in 0 until n) {
        val rate = -D0.get(i, i)
        if (rate > 0) {
            invD0.set(i, i, 1.0 / rate)
        }
    }

    val moments = DoubleArray(kset.size)

    for ((idx, k) in kset.withIndex()) {
        // Compute inv(-D0)^k
        var invD0k = Matrix.eye(n)
        for (p in 0 until k) {
            val temp = Matrix(n, n)
            for (i in 0 until n) {
                for (j in 0 until n) {
                    var sum = 0.0
                    for (l in 0 until n) {
                        sum += invD0k.get(i, l) * invD0.get(l, j)
                    }
                    temp.set(i, j, sum)
                }
            }
            invD0k = temp
        }

        // M_k = pi * inv(-D0)^k * e
        var moment = 0.0
        for (i in 0 until n) {
            var sum = 0.0
            for (j in 0 until n) {
                sum += invD0k.get(i, j)
            }
            moment += pi[i] * sum
        }

        // Multiply by k! for proper moment calculation
        moments[idx] = moment * factorial(k)
    }

    return moments
}

/**
 * Helper function to compute factorial.
 */
private fun factorial(n: Int): Double {
    if (n <= 1) return 1.0
    var result = 1.0
    for (i in 2..n) {
        result *= i.toDouble()
    }
    return result
}

/**
 * Computes the squared coefficient of variation for a deterministic process.
 *
 * @param DET Deterministic MAP as {D0, D1}
 * @return SCV value
 */
fun det_scv(DET: MatrixCell): Double {
    val moments = det_moment(DET, intArrayOf(1, 2))
    val E1 = moments[0]
    val E2 = moments[1]
    return (E2 - E1 * E1) / (E1 * E1)
}

/**
 * Computes the autocorrelation function for a deterministic process.
 *
 * @param DET Deterministic MAP as {D0, D1}
 * @param kset Lags at which to compute ACF
 * @return Array of ACF values
 */
fun det_acf(DET: MatrixCell, kset: IntArray): DoubleArray {
    val D0 = DET[0]
    val n = D0.numRows

    val P = det_embedded(DET)
    val pi = dtmc_solve(P)

    // Compute holding times
    val holdTimes = DoubleArray(n)
    for (i in 0 until n) {
        holdTimes[i] = 1.0 / (-D0.get(i, i))
    }

    // Compute K matrix: K(i,j) = holdTime(i) * P(i,j)
    val K = Matrix(n, n)
    for (i in 0 until n) {
        for (j in 0 until n) {
            K.set(i, j, holdTimes[i] * P.get(i, j))
        }
    }

    val E1 = det_moment(DET, intArrayOf(1))[0]
    val E2 = det_moment(DET, intArrayOf(2))[0]
    val variance = E2 - E1 * E1

    val acf = DoubleArray(kset.size)

    for ((idx, k) in kset.withIndex()) {
        if (k <= 0) {
            acf[idx] = 1.0
            continue
        }

        // Compute P^(k-1)
        var Pk = Matrix.eye(n)
        for (p in 0 until k - 1) {
            val temp = Matrix(n, n)
            for (i in 0 until n) {
                for (j in 0 until n) {
                    var sum = 0.0
                    for (l in 0 until n) {
                        sum += Pk.get(i, l) * P.get(l, j)
                    }
                    temp.set(i, j, sum)
                }
            }
            Pk = temp
        }

        // Compute K * P^(k-1) * K
        val KPk = Matrix(n, n)
        for (i in 0 until n) {
            for (j in 0 until n) {
                var sum = 0.0
                for (l in 0 until n) {
                    sum += K.get(i, l) * Pk.get(l, j)
                }
                KPk.set(i, j, sum)
            }
        }

        val KPkK = Matrix(n, n)
        for (i in 0 until n) {
            for (j in 0 until n) {
                var sum = 0.0
                for (l in 0 until n) {
                    sum += KPk.get(i, l) * K.get(l, j)
                }
                KPkK.set(i, j, sum)
            }
        }

        // rho(k) = (pi * K * P^(k-1) * K * e - E1^2) / variance
        var joint = 0.0
        for (i in 0 until n) {
            var sum = 0.0
            for (j in 0 until n) {
                sum += KPkK.get(i, j)
            }
            joint += pi[i] * sum
        }

        if (variance > 0) {
            acf[idx] = (joint - E1 * E1) / variance
        } else {
            acf[idx] = 0.0
        }
    }

    return acf
}

/**
 * Generates samples from a deterministic MAP.
 *
 * @param DET Deterministic MAP as {D0, D1}
 * @param nSamples Number of samples to generate
 * @param initState Initial state (0-based, or null for random from stationary)
 * @return Triple of (samples, last states, first states)
 */
fun det_sample(
    DET: MatrixCell,
    nSamples: Int,
    initState: Int? = null
): Triple<DoubleArray, IntArray, IntArray> {
    val D0 = DET[0]
    val D1 = DET[1]
    val n = D0.numRows
    val random = Random()

    // Determine initial state
    val startState = if (initState != null) {
        initState
    } else {
        // Sample from interval-stationary distribution
        val pi = map_pie(DET)
        val cumPi = DoubleArray(n)
        cumPi[0] = pi[0]
        for (i in 1 until n) {
            cumPi[i] = cumPi[i - 1] + pi[i]
        }

        val r = random.nextDouble()
        var state = 0
        for (i in 0 until n) {
            if (r <= cumPi[i]) {
                state = i
                break
            }
        }
        state
    }

    val samples = DoubleArray(nSamples)
    val lastStates = IntArray(nSamples)
    val firstStates = IntArray(nSamples)

    // Build transition probabilities
    val transitionProbs = Array(n) { i ->
        val totalRate = -D0.get(i, i)
        val probs = DoubleArray(2 * n)
        if (totalRate > 0) {
            for (j in 0 until n) {
                probs[j] = D0.get(i, j) / totalRate // Hidden transitions
                probs[n + j] = D1.get(i, j) / totalRate // Arrival transitions
            }
            probs[i] = 0.0 // No self-loop without arrival
        }
        probs
    }

    // Compute cumulative probabilities
    val cumProbs = Array(n) { i ->
        val cum = DoubleArray(2 * n)
        cum[0] = transitionProbs[i][0]
        for (j in 1 until 2 * n) {
            cum[j] = cum[j - 1] + transitionProbs[i][j]
        }
        cum
    }

    // Generate samples
    var currentState = startState

    for (s in 0 until nSamples) {
        firstStates[s] = currentState
        var interarrivalTime = 0.0
        var pathStates = mutableListOf(currentState)

        while (true) {
            // Holding time in current state
            val rate = -D0.get(currentState, currentState)
            val holdTime = if (rate > 0) 1.0 / rate else 1.0

            // Sample next transition
            val rnd = random.nextDouble()
            var nextDest = currentState
            for (j in 0 until 2 * n) {
                if (rnd <= cumProbs[currentState][j]) {
                    nextDest = j
                    break
                }
            }

            if (nextDest >= n) {
                // Arrival transition
                interarrivalTime += holdTime
                currentState = nextDest - n
                break
            } else {
                // Hidden transition
                interarrivalTime += holdTime
                pathStates.add(nextDest)
                currentState = nextDest
            }
        }

        samples[s] = interarrivalTime
        lastStates[s] = currentState
    }

    return Triple(samples, lastStates, firstStates)
}

/**
 * Computes the sum of two deterministic processes.
 *
 * @param DET1 First deterministic MAP
 * @param DET2 Second deterministic MAP
 * @return Combined MAP representing sum
 */
fun det_sum(DET1: MatrixCell, DET2: MatrixCell): MatrixCell {
    val D0_1 = DET1[0]
    val D1_1 = DET1[1]
    val D0_2 = DET2[0]
    val D1_2 = DET2[1]

    val n1 = D0_1.numRows
    val n2 = D0_2.numRows
    val n = n1 * n2

    // Kronecker product construction
    // D0 = kron(D0_1, I_2) + kron(I_1, D0_2)
    // D1 = kron(D1_1, D1_2)

    val D0 = Matrix(n, n)
    val D1 = Matrix(n, n)

    for (i1 in 0 until n1) {
        for (j1 in 0 until n1) {
            for (i2 in 0 until n2) {
                for (j2 in 0 until n2) {
                    val i = i1 * n2 + i2
                    val j = j1 * n2 + j2

                    // kron(D0_1, I_2)
                    if (i2 == j2) {
                        D0.set(i, j, D0.get(i, j) + D0_1.get(i1, j1))
                    }

                    // kron(I_1, D0_2)
                    if (i1 == j1) {
                        D0.set(i, j, D0.get(i, j) + D0_2.get(i2, j2))
                    }

                    // kron(D1_1, D1_2)
                    D1.set(i, j, D1_1.get(i1, j1) * D1_2.get(i2, j2))
                }
            }
        }
    }

    val result = MatrixCell(2)
    result[0] = D0
    result[1] = D1

    return result
}
