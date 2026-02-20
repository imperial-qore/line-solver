/**
 * @file Marked Markovian Arrival Process sample generation
 * 
 * Generates random samples from MMAP distributions for each marked class.
 * Essential for simulation, Monte Carlo methods, and empirical analysis.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.io.Ret
import jline.lang.processes.DiscreteSampler
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import java.util.*
import kotlin.math.ln

/**
 * Generates samples of inter-arrival times and event types from a MMAP using a specified number of samples.
 *
 * @param MMAP the MatrixCell representing the MMAP, containing the D0, D1 matrices and additional marking matrices
 * @param n    the number of samples to generate
 * @return an mmapSampleReturn object containing the generated samples and their corresponding types
 */
fun mmap_sample(MMAP: MatrixCell, n: Long, random: Random = Random()): Ret.mamMMAPSample {
    val D0 = MMAP[0]
    val D1 = MMAP[1]
    val C = MMAP.size() - 2
    val order = D0.numRows
    val typeSampler = Array(order) { arrayOfNulls<DiscreteSampler>(order) }
    val x = Matrix(1, C)
    for (c in 0..<C) {
        x[c] = c.toDouble()
    }
    for (i in 0..<order) {
        for (j in 0..<order) {
            if (D1[i, j] > 0) {
                val pij = Matrix(1, C)
                for (c in 0..<C) {
                    pij[c] = MMAP[2 + c][i, j]
                }
                pij.scaleEq(1 / pij.elementSum())
                if (pij.elementSum() > 0) {
                    typeSampler[i][j] = DiscreteSampler(pij, x)
                }
            }
        }
    }
    val samples = DoubleArray(n.toInt())
    val types = IntArray(n.toInt())
    val states = IntArray(n.toInt())
    if (D0.numElements == 1) { // exponential distribution
        val lambda = D1.value()
        for (i in 0..<n) {
            samples[i.toInt()] = -ln(random.nextDouble()) / lambda
            types[i.toInt()] = 1
            states[i.toInt()] = 0 // Single state
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
            // Record state at the time of this event (like MATLAB STS(h) = s)
            states[i.toInt()] = currentState
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
                            types[i.toInt()] = typeSampler[currentState][nextState]!!.sample().toInt()
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
    return Ret.mamMMAPSample(samples, C, types, states)
}


/* Sample from a MMAP
            Example:
            Matrix D0 = new Matrix("[-2,1;0,-3]");
            Matrix D1 = new Matrix("[1,0;2,1]");
            Matrix D11 = new Matrix("[1,0;2,0]");
            Matrix D12 = new Matrix("[0,0;0,1]");
            MatrixCell MAP = new MatrixCell(D0, D1);
            MAP.set(0, D0);
            MAP.set(1, D1);
            MatrixCell MMAP = new MatrixCell(D0, D1);
            MMAP.set(0, D0);
            MMAP.set(1, D1);
            MMAP.set(2, D11);
            MMAP.set(3, D12);
            mmapSampleReturn samples = mmap_sample(MMAP, 100000);
         */
/**
 * MMAP sample algorithms
 */
@Suppress("unused")
class MmapSampleAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}