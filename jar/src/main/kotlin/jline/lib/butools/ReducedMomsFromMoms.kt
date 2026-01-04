package jline.lib.butools

import jline.util.matrix.Matrix

/**
 * Returns the reduced moments given the raw moments.
 *
 * The raw moments are: m_i = E(X^i)
 * The reduced moments are: r_i = m_i / i!
 *
 * @param m The list of raw moments (starting with the first moment)
 * @return The list of reduced moments
 */
fun reducedMomsFromMoms(m: Matrix): Matrix {
    val rm = Matrix(m.numRows, m.numCols, m.numRows * m.numCols)
    var invFactorial = 1.0

    for (i in 0 until m.length()) {
        invFactorial /= (i + 1) // Calculate 1/(i+1)!
        rm[i] = m[i] * invFactorial
    }

    return rm
}

/**
 * Returns the reduced moments given the raw moments (DoubleArray version).
 *
 * The raw moments are: m_i = E(X^i)
 * The reduced moments are: r_i = m_i / i!
 *
 * @param m The list of raw moments (starting with the first moment)
 * @return The list of reduced moments
 */
fun ReducedMomsFromMoms(m: DoubleArray): DoubleArray {
    val rm = DoubleArray(m.size)
    var invFactorial = 1.0

    for (i in m.indices) {
        invFactorial /= (i + 1) // Calculate 1/(i+1)!
        rm[i] = m[i] * invFactorial
    }

    return rm
}