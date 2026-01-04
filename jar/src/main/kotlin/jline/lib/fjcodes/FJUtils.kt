package jline.lib.fjcodes

import jline.util.matrix.Matrix

/**
 * FJ_codes utility functions
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 *
 * Based on FJ_codes MATLAB toolkit:
 * Z. Qiu, J.F. Pérez, and P. Harrison, "Beyond the Mean in Fork-Join Queues:
 * Efficient Approximation for Response-Time Tails", IFIP Performance 2015.
 */

/**
 * Compute Kronecker sum of two matrices
 *
 * Returns: A ⊕ B = A \otimes I + I \otimes B
 *
 * @param A First matrix
 * @param B Second matrix
 * @return Kronecker sum A ⊕ B
 */
fun kronsum(A: Matrix, B: Matrix): Matrix {
    val IA = Matrix.eye(A.numRows)
    val IB = Matrix.eye(B.numRows)
    return A.kron(IB).add(1.0, IA.kron(B))
}

/**
 * Find the row index in a matrix that matches a given row vector
 *
 * Assumes there is exactly one matching row in the matrix.
 *
 * @param row Row vector to search for
 * @param matrix Matrix to search in
 * @return 1-based index of the matching row, or -1 if not found
 */
fun vectmatch(row: DoubleArray, matrix: Matrix): Int {
    val m = matrix.numRows
    val n = matrix.numCols

    for (outer in 0 until m) {
        var matches = true
        for (inner in 0 until n) {
            if (Math.abs(matrix.get(outer, inner) - row[inner]) > 1e-10) {
                matches = false
                break
            }
        }
        if (matches) {
            return outer + 1  // Return 1-based index (MATLAB convention)
        }
    }

    return -1  // Not found
}

/**
 * Build combinatorial index patterns
 *
 * Generates all ways to distribute 'cr' items into 'm' bins.
 * For example, build_index(2, 1) with m=2, cr=1 gives:
 * [1, 0]
 * [0, 1]
 *
 * @param m Number of bins/phases
 * @param cr Number of items to distribute
 * @return Matrix where each row is a distribution pattern
 */
fun build_index(m: Int, cr: Int): Matrix {
    // Compute total number of combinations: C(cr+m-1, cr)
    val totalDim = binomialCoeff(cr + m - 1, cr)
    val indexes = Matrix(totalDim, m)

    // First row: all items in first bin
    indexes.set(0, 0, cr.toDouble())

    // Generate remaining rows
    for (row in 1 until totalDim) {
        // Find first non-zero element in previous row
        var k = -1
        for (col in 0 until m) {
            if (indexes.get(row - 1, col) > 0.0) {
                k = col
                break
            }
        }

        if (k >= 0 && k < m - 1) {
            // Copy previous row
            for (col in 0 until m) {
                indexes.set(row, col, indexes.get(row - 1, col))
            }

            // Shift one item from position k to k+1
            val currentVal = indexes.get(row, k + 1)
            indexes.set(row, k + 1, currentVal + 1.0)
            val kVal = indexes.get(row, k)
            indexes.set(row, 0, kVal - 1.0)

            // Zero out positions 1 to k
            for (col in 1..k) {
                indexes.set(row, col, 0.0)
            }
        }
    }

    return indexes
}

/**
 * Extract a row from a Matrix as a DoubleArray
 *
 * @param matrix Source matrix
 * @param row Row index (0-based)
 * @return Row as DoubleArray
 */
fun getRowAsArray(matrix: Matrix, row: Int): DoubleArray {
    val rowMatrix = matrix.getRow(row)
    val result = DoubleArray(rowMatrix.length())
    for (i in 0 until rowMatrix.length()) {
        result[i] = rowMatrix.get(i)
    }
    return result
}

/**
 * Copy a submatrix into a target matrix at specified position
 *
 * @param target Target matrix to copy into
 * @param startRow Starting row in target (0-based)
 * @param startCol Starting column in target (0-based)
 * @param source Source matrix to copy from
 */
fun setSubMatrix(target: Matrix, startRow: Int, startCol: Int, source: Matrix) {
    for (i in 0 until source.numRows) {
        for (j in 0 until source.numCols) {
            target.set(startRow + i, startCol + j, source.get(i, j))
        }
    }
}

/**
 * Compute binomial coefficient C(n, k) = n! / (k! * (n-k)!)
 */
private fun binomialCoeff(n: Int, k: Int): Int {
    if (k > n) return 0
    if (k == 0 || k == n) return 1

    var result: Long = 1
    val kMin = Math.min(k, n - k)

    for (i in 0 until kMin) {
        result = result * (n - i).toLong() / (i + 1).toLong()
    }

    return result.toInt()
}
