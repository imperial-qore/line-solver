/**
 * @file M3PP superposition trace-based counting process fitting
 *
 * Superposes k M3PP processes to fit a multi-class trace with m classes.
 * Each class is fitted independently and then superposed.
 *
 * @since LINE 3.0
 */
package jline.api.mam.m3pp

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import jline.api.mam.Mmpp2_fitc

/**
 * Superposes k M3PP processes to fit a multi-class trace with m classes.
 *
 * @param T Inter-arrival times as array
 * @param A Class labels as array
 * @param t Optional finite time scale
 * @param tinf Optional near-infinite time scale
 * @return Pair of (fitted superposed M3PP, individual M3PP processes)
 */
@JvmOverloads
fun m3pp_superpos_fitc_trace(
    T: DoubleArray,
    A: IntArray,
    t: Double? = null,
    tinf: Double? = null
): Pair<Array<Matrix>, List<Array<Matrix>>> {
    // Default time scales
    val tUsed = t ?: (10 * T.average())
    val tinfUsed = tinf ?: maxOf(10 * tUsed, (T.sum() - T.first()) / 100)

    // Total rate
    val a = 1.0 / T.average()

    // Per-class rates
    val L = A.distinct().sorted()
    val m = L.size
    val pv = DoubleArray(m) { i ->
        A.count { it == L[i] } / A.size.toDouble()
    }
    val av = DoubleArray(m) { pv[it] * a }

    // Compute counting process at resolution t and tinf
    val Nt = mtrace_iat2counts(T, A, tUsed)
    val Ninf = mtrace_iat2counts(T, A, tinfUsed)

    // Compute IDC(t) and IDC(inf) for each class
    val btv = DoubleArray(m)
    val binfv = DoubleArray(m)
    val m3tv = DoubleArray(m)

    for (i in 0 until m) {
        val classCol = getColumn(Nt, i)
        val classColInf = getColumn(Ninf, i)

        btv[i] = variance(classCol) / (av[i] * tUsed)
        binfv[i] = variance(classColInf) / (av[i] * tinfUsed)

        // Third centered moment
        val mt = doubleArrayOf(
            classCol.average(),
            classCol.map { it * it }.average(),
            classCol.map { it * it * it }.average()
        )
        m3tv[i] = mt[2] - 3 * mt[1] * mt[0] + 2 * mt[0] * mt[0] * mt[0]
    }

    // Fit superposition
    return m3pp_superpos_fitc(av, btv, binfv, m3tv, tUsed, tinfUsed)
}

/**
 * Superposes k M3PP processes to fit a multi-class trace from Matrix inputs.
 *
 * @param T Inter-arrival times as Matrix
 * @param A Class labels as Matrix
 * @param t Optional finite time scale
 * @param tinf Optional near-infinite time scale
 * @return Pair of (fitted superposed M3PP, individual M3PP processes)
 */
@JvmOverloads
fun m3pp_superpos_fitc_trace(
    T: Matrix,
    A: Matrix,
    t: Double? = null,
    tinf: Double? = null
): Pair<Array<Matrix>, List<Array<Matrix>>> {
    val tArray = T.toArray1D()
    val aArray = IntArray(A.numRows * A.numCols) { i -> A.toArray1D()[i].toInt() }
    return m3pp_superpos_fitc_trace(tArray, aArray, t, tinf)
}

/**
 * Superposes k individual M3PP processes.
 *
 * @param av Per-class arrival rates
 * @param btv Per-class IDC at time t
 * @param binfv Per-class IDC at time infinity
 * @param m3tv Per-class third centered moments
 * @param t Finite time scale
 * @param tinf Near-infinite time scale
 * @return Pair of (fitted superposed M3PP, individual M3PP processes)
 */
fun m3pp_superpos_fitc(
    av: DoubleArray,
    btv: DoubleArray,
    binfv: DoubleArray,
    m3tv: DoubleArray,
    t: Double,
    tinf: Double
): Pair<Array<Matrix>, List<Array<Matrix>>> {
    val m = av.size
    val m3pps = mutableListOf<Array<Matrix>>()

    // Fit individual M3PP for each class
    for (i in 0 until m) {
        try {
            // For single-class fitting, use MMPP2 fitting
            val mmpp = Mmpp2_fitc.mmpp2_fitc(av[i], btv[i], btv[i], binfv[i], m3tv[i], t, tinf)
            m3pps.add(mmpp)
        } catch (e: Exception) {
            // Fallback to exponential process if fitting fails
            val D0 = Matrix(1, 1)
            D0[0, 0] = -av[i]
            val D1 = Matrix(1, 1)
            D1[0, 0] = av[i]
            m3pps.add(arrayOf(D0, D1))
        }
    }

    // Superpose individual M3PPs
    val fit = superposeMaps(m3pps)

    return Pair(fit, m3pps)
}

/**
 * Superposes multiple MAPs into one.
 */
private fun superposeMaps(maps: List<Array<Matrix>>): Array<Matrix> {
    if (maps.isEmpty()) {
        throw IllegalArgumentException("Cannot superpose empty list of MAPs")
    }

    if (maps.size == 1) {
        return maps[0]
    }

    // Start with the first MAP
    var result = maps[0]

    // Superpose each additional MAP
    for (i in 1 until maps.size) {
        result = superposeTwoMaps(result, maps[i])
    }

    return result
}

/**
 * Superposes two MAPs using Kronecker sum/product.
 */
private fun superposeTwoMaps(map1: Array<Matrix>, map2: Array<Matrix>): Array<Matrix> {
    val D0a = map1[0]
    val D1a = map1[1]
    val D0b = map2[0]
    val D1b = map2[1]

    val na = D0a.numRows
    val nb = D0b.numRows
    val n = na * nb

    // Kronecker sum for D0: D0a ⊕ D0b = D0a \otimes I_b + I_a \otimes D0b
    val Ia = Matrix.eye(na)
    val Ib = Matrix.eye(nb)
    val D0 = D0a.kron(Ib).add(Ia.kron(D0b))

    // Kronecker product for D1: D1a \otimes I_b + I_a \otimes D1b
    val D1 = D1a.kron(Ib).add(Ia.kron(D1b))

    return arrayOf(D0, D1)
}

// Helper functions

private fun getColumn(matrix: Array<DoubleArray>, col: Int): DoubleArray {
    return DoubleArray(matrix.size) { row -> matrix[row][col] }
}

private fun variance(arr: DoubleArray): Double {
    if (arr.isEmpty()) return 0.0
    val mean = arr.average()
    return arr.map { (it - mean) * (it - mean) }.average()
}

/**
 * M3PP superposition fitc trace algorithms
 */
@Suppress("unused")
class M3ppSuperposFitcTraceAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}
