package jline.lib.butools

import jline.util.matrix.Matrix

fun MomsFromFactorialMoms(fm: Matrix): Matrix {
    val n = fm.length()
    val m = Matrix(1, n, n)
    m[0] = fm[0]

    for (i in 1..<n) {
        // Compute polynomial with roots [0, 1, 2, ..., i-1]
        val rootsMatrix = Matrix(1, i, i)
        for (k in 0 until i) {
            rootsMatrix[k] = k.toDouble()
        }

        // Get polynomial coefficients and negate
        val ehMatrix = poly(rootsMatrix)
        val eh = DoubleArray(ehMatrix.numCols) { k -> -ehMatrix.get(0, k) }

        // Extract coefficients: eh(end-1:-1:2) in MATLAB
        // means from index size-1 down to 1 in 0-indexed, which gives size-1 elements
        // But we want from size-2 down to 1, which gives size-2 elements = i elements
        // Reverse and skip first: take indices i down to 1
        val eh_coeff = DoubleArray(i) { k -> eh[i - k] }

        // Compute dot product: eh_coeff[0..i-1] dot m[0..i-1]
        var sum = 0.0
        for (k in 0..<i) {
            sum += eh_coeff[k] * m[k]
        }

        m[i] = fm[i] + sum
    }

    if (fm.numRows > fm.numCols) {
        return m.transpose()
    }
    return m
}