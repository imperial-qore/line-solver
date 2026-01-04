/**
 * @file Marked Markovian Arrival Process random generation
 * 
 * Generates random MMAP representations for testing and simulation purposes.
 * Essential for statistical validation and stochastic analysis of multiclass systems.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import java.util.*

/**
 * Generates a random MMAP (Marked Markovian Arrival Process) with a given order and number of classes.
 *
 * @param order   the number of phases (order) in the MAP
 * @param classes the number of different classes (types) of arrivals in the MMAP
 * @return a MatrixCell representing the MMAP, containing the D0, D1 matrices, and class-specific matrices
 */
fun mmap_rand(order: Int, classes: Int): MatrixCell {
    val MMAP = MatrixCell()
    for (c in 0..<2 + classes) {
        MMAP[c] = Matrix(order, order)
    }
    val MAP = map_rand(order)

    MMAP[0] = MAP[0]
    MMAP[1] = MAP[1]

    val rand = Random()
    for (i in 0..<order) {
        val p = DoubleArray(classes)
        var sum = 0.0

        for (c in 0..<classes) {
            p[c] = rand.nextDouble()
            sum += p[c]
        }

        for (c in 0..<classes) {
            p[c] /= sum
        }

        for (j in 0..<order) {
            for (c in 0..<classes) {
                MMAP[2 + c] = MMAP[1].copy()
                MMAP[2 + c].scaleEq(p[c])
            }
        }
    }

    return MMAP
}
/**
 * MMAP rand algorithms
 */
@Suppress("unused")
class MmapRandAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}