/**
 * @file Markovian Arrival Process sample generation
 * 
 * Generates random samples from MAP distributions for simulation and empirical analysis.
 * Essential for stochastic simulation, Monte Carlo methods, and model validation studies.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import java.util.*
import kotlin.math.ln

/**
 * Generates samples of inter-arrival times from a MAP using a specified number of samples and a random generator.
 *
 *
 * This method generates `n` samples of inter-arrival times from the MAP defined by matrices D0 and D1.
 * If the MAP has a single phase, the samples are generated from an exponential distribution.
 * For multi-phase MAPs, the sampling accounts for transitions between states/phases.
 *
 * @param MAP    the MatrixCell representing the MAP, containing the D0 and D1 matrices
 * @param n      the number of samples to generate
 * @param random the random number generator to use
 * @return an array of doubles containing the generated samples
 */

fun map_sample(MAP: MatrixCell?, n: Long, random: Random?): DoubleArray {
    return map_sample(MAP, n, random)
}

/**
 * Generates samples of inter-arrival times from a MAP using a specified number of samples and a random generator.
 *
 * @param D0     the hidden transition matrix of the MAP
 * @param D1     the visible transition matrix of the MAP
 * @param n      the number of samples to generate
 * @param random the random number generator to use
 * @return an array of doubles containing the generated samples
 */

fun map_sample(D0: Matrix, D1: Matrix, n: Long, random: Random): DoubleArray {
    val samples = DoubleArray(n.toInt())
    if (D0.numElements == 1) { // exponential distribution
        val lambda = D1.value()
        for (i in 0..<n) {
            samples[i.toInt()] = -ln(random.nextDouble()) / lambda
        }
    } else {
        val nphases = D0.numCols
        val pie = map_pie(D0, D1)
        var currentState = pie.numRows
        var sum = 0.0
        var r = random.nextDouble()
        for (i in 0..<nphases) {
            sum += pie[0, i]
            if (r < sum) {
                currentState = i
                break
            }
        }

        val row = DoubleArray(2 * nphases)
        for (i in 0..<n) {
            var continue_sample = true
            samples[i.toInt()] = 0.0
            while (continue_sample) {
                val rate = -D0[currentState, currentState]
                val time = -ln(random.nextDouble()) / rate
                samples[i.toInt()] += time
                for (k in 0..<nphases) {
                    row[k] = D0[currentState, k]
                    row[nphases + k] = D1[currentState, k]
                }
                row[currentState] = 0.0
                var nextState = 2 * nphases - 1
                sum = 0.0
                r = random.nextDouble()
                for (j in 0..<2 * nphases) {
                    sum += row[j] / rate
                    if (r < sum) {
                        if (j >= nphases) {
                            nextState = j - nphases
                            continue_sample = false
                            break
                        } else {
                            nextState = j
                            continue_sample = true
                            break
                        }
                    }
                }
                currentState = nextState
            }
        }
    }
    return samples
}
/**
 * Result of BMAP sampling containing inter-arrival time and batch size.
 */
data class BmapSample(val interarrivalTime: Double, val batchSize: Int)

/**
 * Generates samples from a BMAP (Batch Markovian Arrival Process).
 *
 * This method generates `n` samples from the BMAP, returning both inter-arrival times
 * and batch sizes. The BMAP is represented as a MatrixCell containing:
 * - bmap.get(0) = D0 (transitions without arrivals)
 * - bmap.get(1) = D1_total (sum of all Dk - total arrival rate)
 * - bmap.get(2) = D1 (transitions generating 1 arrival)
 * - bmap.get(3) = D2 (transitions generating 2 arrivals)
 * - ...
 * - bmap.get(k+1) = Dk (transitions generating k arrivals)
 *
 * @param bmap   the MatrixCell representing the BMAP
 * @param n      the number of samples to generate
 * @param random the random number generator to use
 * @return an array of BmapSample containing inter-arrival times and batch sizes
 */
fun bmap_sample(bmap: MatrixCell, n: Long, random: Random): Array<BmapSample> {
    val D0 = bmap.get(0)
    val D1_total = bmap.get(1)
    val maxBatchSize = bmap.size() - 2

    // Collect individual Dk matrices
    val Dk = Array(maxBatchSize) { k -> bmap.get(k + 2) }

    return bmap_sample(D0, D1_total, Dk, n, random)
}

/**
 * Generates samples from a BMAP using D0, D1_total, and individual Dk matrices.
 *
 * @param D0        the hidden transition matrix (no arrivals)
 * @param D1_total  the sum of all Dk matrices (total arrival rate)
 * @param Dk        array of Dk matrices where Dk[k-1] generates batch size k
 * @param n         the number of samples to generate
 * @param random    the random number generator to use
 * @return an array of BmapSample containing inter-arrival times and batch sizes
 */
fun bmap_sample(D0: Matrix, D1_total: Matrix, Dk: Array<Matrix>, n: Long, random: Random): Array<BmapSample> {
    val nphases = D0.numCols
    val maxBatchSize = Dk.size

    // Handle degenerate case: single phase (exponential)
    if (D0.numElements == 1) {
        val lambda = D1_total.value()
        // For single phase, determine batch size distribution
        val batchRates = DoubleArray(maxBatchSize) { k -> Dk[k].value() }
        val totalRate = batchRates.sum()

        return Array(n.toInt()) {
            val interarrivalTime = -ln(random.nextDouble()) / lambda

            // Select batch size proportionally to rates
            var batchSize = 1
            if (totalRate > 0) {
                val r = random.nextDouble() * totalRate
                var cumSum = 0.0
                for (k in 0 until maxBatchSize) {
                    cumSum += batchRates[k]
                    if (r < cumSum) {
                        batchSize = k + 1
                        break
                    }
                }
            }
            BmapSample(interarrivalTime, batchSize)
        }
    }

    // Multi-phase BMAP sampling
    val samples = Array(n.toInt()) { BmapSample(0.0, 1) }

    // Initialize state from stationary distribution
    val pie = map_pie(D0, D1_total)
    var currentState = nphases - 1
    var sum = 0.0
    var r = random.nextDouble()
    for (i in 0 until nphases) {
        sum += pie[0, i]
        if (r < sum) {
            currentState = i
            break
        }
    }

    // Build combined transition row: [D0 transitions | D1 | D2 | ... | Dk]
    // Total columns = nphases + maxBatchSize * nphases
    val totalCols = nphases + maxBatchSize * nphases
    val row = DoubleArray(totalCols)

    for (i in 0 until n.toInt()) {
        var continuesSample = true
        var interarrivalTime = 0.0
        var batchSize = 1

        while (continuesSample) {
            val rate = -D0[currentState, currentState]
            val time = -ln(random.nextDouble()) / rate
            interarrivalTime += time

            // Build probability row
            // First nphases entries: D0 transitions (hidden, no arrival)
            for (k in 0 until nphases) {
                row[k] = D0[currentState, k]
            }
            row[currentState] = 0.0  // No self-transition

            // Next entries: Dk transitions for each batch size
            for (b in 0 until maxBatchSize) {
                for (k in 0 until nphases) {
                    row[nphases + b * nphases + k] = Dk[b][currentState, k]
                }
            }

            // Select next state and batch size
            var nextState = currentState
            sum = 0.0
            r = random.nextDouble()

            for (j in 0 until totalCols) {
                sum += row[j] / rate
                if (r < sum) {
                    if (j < nphases) {
                        // D0 transition - no arrival, continue sampling
                        nextState = j
                        continuesSample = true
                    } else {
                        // Dk transition - arrival with batch size
                        val idx = j - nphases
                        batchSize = (idx / nphases) + 1
                        nextState = idx % nphases
                        continuesSample = false
                    }
                    break
                }
            }
            currentState = nextState
        }

        samples[i] = BmapSample(interarrivalTime, batchSize)
    }

    return samples
}

/**
 * MAP sample algorithms
 */
@Suppress("unused")
class MapSampleAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}