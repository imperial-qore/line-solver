/**
 * @file M3A MMAP Compression
 *
 * Compresses the representation of a Marked Markovian Arrival Process
 * based on the M3A toolbox.
 *
 * @since LINE 3.0
 */
package jline.lib.m3a

import jline.api.mam.mamap2m_fit_gamma_fb_mmap
import jline.api.mam.mmap_isfeasible
import jline.util.matrix.MatrixCell

/**
 * Compression method enumeration
 */
enum class M3aCompressMethod {
    /** Two-state acyclic MAP compression (AMAP) */
    AMAP_2STATE
}

/**
 * Options for M3A compression
 */
data class M3aCompressOptions(
    /** Compression method to use */
    val method: M3aCompressMethod = M3aCompressMethod.AMAP_2STATE,
    /** Number of states in the compressed representation */
    val numStates: Int = 2,
    /** Verbose output flag */
    val verbose: Boolean = false
)

/**
 * Compresses a Marked Markovian Arrival Process (MMAP) using M3A fitting.
 *
 * This function takes an arbitrary-order MMAP and produces a compressed
 * second-order acyclic MMAP that approximates the original process.
 *
 * The compression preserves key statistical characteristics:
 * - First three moments of inter-arrival times
 * - Autocorrelation structure (via gamma decay rate)
 * - Class probabilities
 * - Forward and backward moments
 *
 * @param mmap Input MMAP to compress, stored as MatrixCell where:
 *             - mmap[0] = D0 (transition matrix without arrivals)
 *             - mmap[1] = D1 (aggregate arrival matrix)
 *             - mmap[2+c] = D1c (class c arrival matrix) for c = 0, 1, ..., numClasses-1
 * @param options Compression options (defaults to 2-state AMAP compression)
 * @return Compressed MMAP representation
 * @throws IllegalArgumentException if the input MMAP is invalid
 *
 * Example:
 * ```kotlin
 * // Create an MMAP
 * val mmap = MatrixCell(4)
 * mmap[0] = Matrix(3, 3) // D0
 * mmap[1] = Matrix(3, 3) // D1
 * mmap[2] = Matrix(3, 3) // D11
 * mmap[3] = Matrix(3, 3) // D12
 *
 * // Compress to 2-state representation
 * val compressed = m3afit_compress(mmap)
 * ```
 *
 * References:
 * - A. Sansottera, G. Casale, P. Cremonesi. "Fitting Second-Order Acyclic
 *   Marked Markovian Arrival Processes." IEEE/IFIP DSN 2013.
 * - G. Casale, A. Sansottera, P. Cremonesi. "Compact Markov-Modulated
 *   Models for Multiclass Trace Fitting." European Journal of Operations
 *   Research, 2016.
 */
fun m3afit_compress(
    mmap: MatrixCell,
    options: M3aCompressOptions = M3aCompressOptions()
): MatrixCell {
    // Validate input
    require(mmap.size() >= 3) {
        "Input MMAP must have at least 3 matrices (D0, D1, and at least one class matrix)"
    }

    val numClasses = mmap.size() - 2

    if (options.verbose) {
        val mmapType = "${options.numStates}-state AMAP[$numClasses]"
        println("Init: M3A will search for a $mmapType")
    }

    // Perform compression based on method
    val result = when (options.method) {
        M3aCompressMethod.AMAP_2STATE -> {
            // Two-state AMAP compression using forward-backward fitting
            mamap2m_fit_gamma_fb_mmap(mmap)
        }
    }

    // Validate result
    if (options.verbose) {
        val resultType = "${result[0].numRows}-state M3PP[$numClasses]"
        if (mmap_isfeasible(result)) {
            println("Output: M3A found a valid $resultType.")
        } else {
            println("Output: M3A could *not* obtain a valid MMAP.")
        }
    }

    return result
}

/**
 * Compresses an MMAP with default options.
 * Convenience overload for Java interoperability.
 *
 * @param mmap Input MMAP to compress
 * @return Compressed MMAP representation
 */
@JvmOverloads
fun m3afit_compress(mmap: MatrixCell): MatrixCell {
    return m3afit_compress(mmap, M3aCompressOptions())
}

/**
 * M3A compress algorithms class marker for Dokka
 */
@Suppress("unused")
class M3aCompressAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}
