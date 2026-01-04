/**
 * @file M3A automatic fitting functions
 *
 * Top-level API for automatic MMAP fitting from traces using the M3A toolbox.
 * Provides m3afit_init and m3afit_auto functions for automatic model selection.
 *
 * @since LINE 3.0
 */
package jline.lib.m3a

import jline.GlobalConstants
import jline.VerboseLevel
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import jline.api.mam.*
import jline.api.mam.m3pp.*

/** Helper to check if verbose output is enabled */
private fun isVerbose(): Boolean = GlobalConstants.getVerbose() != VerboseLevel.SILENT

/**
 * Data structure for multiclass trace representation.
 *
 * @property S Inter-arrival times
 * @property C Class labels for each arrival
 * @property numClasses Number of distinct classes
 */
data class MTrace(
    val S: DoubleArray,
    val C: IntArray,
    val numClasses: Int
) {
    override fun equals(other: Any?): Boolean {
        if (this === other) return true
        if (other !is MTrace) return false
        return S.contentEquals(other.S) && C.contentEquals(other.C) && numClasses == other.numClasses
    }

    override fun hashCode(): Int {
        var result = S.contentHashCode()
        result = 31 * result + C.contentHashCode()
        result = 31 * result + numClasses
        return result
    }
}

/**
 * Options for M3A fitting algorithms.
 *
 * @property method Fitting method (0 = inter-arrival, 1 = counting process)
 * @property numStates Number of states for the fitted MMAP (mandatory)
 * @property timescale Finite time scale for counting process (auto-computed if null)
 * @property timescaleAsy Near-infinite time scale (auto-computed if null)
 */
data class M3aFitOptions(
    val method: Int = 1,
    val numStates: Int,
    val timescale: Double? = null,
    val timescaleAsy: Double? = null
)

/**
 * Prepares multiclass trace for M3A fitting.
 *
 * @param S Inter-arrival times
 * @param C Class number for each arrival
 * @return MTrace structure ready for fitting
 */
fun m3afit_init(S: DoubleArray, C: IntArray): MTrace {
    val numClasses = C.distinct().size
    return MTrace(S, C, numClasses)
}

/**
 * Prepares multiclass trace for M3A fitting from Matrix inputs.
 *
 * @param S Inter-arrival times as Matrix
 * @param C Class number for each arrival as Matrix
 * @return MTrace structure ready for fitting
 */
fun m3afit_init(S: Matrix, C: Matrix): MTrace {
    val sArray = S.toArray1D()
    val cArray = IntArray(C.numRows * C.numCols) { i -> C.toArray1D()[i].toInt() }
    return m3afit_init(sArray, cArray)
}

/**
 * Automatic fitting of trace into a Marked Markovian Arrival Process.
 *
 * Based on the M3A toolbox, this function selects the appropriate fitting algorithm
 * based on the number of classes, requested states, and fitting method.
 *
 * References:
 * [1] A. Sansottera, G. Casale, P. Cremonesi. Fitting Second-Order Acyclic
 *     Marked Markovian Arrival Processes. IEEE/IFIP DSN 2013.
 * [2] G. Casale, A. Sansottera, P. Cremonesi. Compact Markov-Modulated
 *     Models for Multiclass Trace Fitting. European Journal of Operations
 *     Research, 2016.
 *
 * @param mtrace Data structure returned by m3afit_init
 * @param options Fitting options including method and number of states
 * @return Fitted MMAP or null if fitting fails
 */
fun m3afit_auto(mtrace: MTrace, options: M3aFitOptions): MatrixCell? {
    // Compute time scales for counting process method
    val timescale = options.timescale ?: (10 * mtrace.S.average())
    val timescaleAsy = options.timescaleAsy ?: maxOf(10 * timescale, (mtrace.S.sum() - mtrace.S.first()) / 100)

    val mmapType: String
    val mmap: MatrixCell?

    when {
        // 2-class, 2-state, inter-arrival fitting
        mtrace.numClasses == 2 && options.numStates == 2 && options.method == 0 -> {
            mmapType = "${options.numStates}-state MAMAP[${mtrace.numClasses}]"
            if (isVerbose()) println("Init: M3A algorithm will search for a $mmapType")
            mmap = mamap22_fit_gamma_fs_trace(
                Matrix(mtrace.S.size, 1).also { m ->
                    for (i in mtrace.S.indices) m[i, 0] = mtrace.S[i]
                },
                Matrix(mtrace.C.size, 1).also { m ->
                    for (i in mtrace.C.indices) m[i, 0] = mtrace.C[i].toDouble()
                }
            )
        }

        // Multi-class, 2-state, inter-arrival fitting
        mtrace.numClasses > 2 && options.numStates >= 2 && options.method == 0 -> {
            mmapType = "${options.numStates}-state MAMAP[${mtrace.numClasses}]"
            if (isVerbose()) println("Init: M3A algorithm will search for a $mmapType")
            mmap = mamap2m_fit_trace(mtrace.S, mtrace.C)
        }

        // Multi-class, 2-state, counting process fitting
        mtrace.numClasses >= 2 && options.numStates == 2 && options.method == 1 -> {
            mmapType = "${options.numStates}-state M3PP[${mtrace.numClasses}]"
            if (isVerbose()) println("Init: M3A algorithm will search for a $mmapType")
            val result = m3pp2m_fitc_trace(mtrace.S, mtrace.C, "approx_ag", timescale, timescaleAsy)
            mmap = arrayToMatrixCell(result)
        }

        // Multi-class, >2-state, counting process fitting (superposition)
        mtrace.numClasses >= 2 && options.numStates > 2 && options.method == 1 -> {
            mmapType = "${options.numStates}-state M3PP[${mtrace.numClasses}]"
            if (isVerbose()) println("Init: M3A algorithm will search for a $mmapType")
            val (result, _) = m3pp_superpos_fitc_trace(mtrace.S, mtrace.C, timescale, timescaleAsy)
            mmap = arrayToMatrixCell(result)
        }

        else -> {
            if (isVerbose()) println("Output: M3A algorithm could *not* obtain a valid MMAP.")
            return null
        }
    }

    // Validate result
    val mmapResType = "${mmap[0].numRows}-state M3PP[${mmap.size() - 2}]"
    if (isVerbose()) {
        if (mmap_isfeasible(mmap)) {
            println("Output: M3A algorithm found a valid $mmapResType.")
        } else {
            println("Output: M3A algorithm could *not* obtain a valid MMAP.")
        }
    }

    return mmap
}

/**
 * Automatic fitting with simple parameters.
 *
 * @param S Inter-arrival times
 * @param C Class labels
 * @param numStates Number of states for the fitted MMAP
 * @param method Fitting method (0 = inter-arrival, 1 = counting process)
 * @return Fitted MMAP or null if fitting fails
 */
@JvmOverloads
fun m3afit_auto(S: DoubleArray, C: IntArray, numStates: Int, method: Int = 1): MatrixCell? {
    val mtrace = m3afit_init(S, C)
    val options = M3aFitOptions(method = method, numStates = numStates)
    return m3afit_auto(mtrace, options)
}

/**
 * Converts Array<Matrix> to MatrixCell.
 */
private fun arrayToMatrixCell(array: Array<Matrix>): MatrixCell {
    val result = MatrixCell(array.size)
    for (i in array.indices) {
        result[i] = array[i]
    }
    return result
}

/**
 * Wrapper for mamap2m_fit_trace using arrays.
 */
private fun mamap2m_fit_trace(S: DoubleArray, C: IntArray): MatrixCell {
    return mamap2m_fit_gamma_fb_trace(S, C)
}

