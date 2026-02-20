package jline.lib.butools

import jline.util.matrix.Matrix

/**
 * Constructs a Hankel matrix from the given first column and last row.
 * A Hankel matrix has constant anti-diagonals.
 */
private fun hankel(col: DoubleArray, row: DoubleArray): Matrix {
    val n = col.size
    val m = row.size
    val H = Matrix.zeros(n, m)
    for (i in 0 until n) {
        for (j in 0 until m) {
            if (i + j < n) {
                H[i, j] = col[i + j]
            } else {
                H[i, j] = row[i + j - n + 1]
            }
        }
    }
    return H
}

/**
 * Returns the Hankel moments given the raw moments.
 *
 * The raw moments are: m_i = E(X^i)
 *
 * The ith Hankel moment is the determinant of matrix Delta_{i/2}, if i is even,
 * and it is the determinant of Delta^(1)_{(i+1)/2}, if i is odd.
 * For the definition of matrices Delta and Delta^(1) see the Stieltjes moment problem.
 *
 * @param m The list of raw moments (starting with the first moment)
 * @return The list of Hankel moments
 *
 * Reference: http://en.wikipedia.org/wiki/Stieltjes_moment_problem
 */
fun hankelMomsFromMoms(m: Matrix): Matrix {
    val n = m.length()
    val hm = Matrix(m.numRows, m.numCols, n)

    // Extract raw moments as array for convenience
    val mArr = DoubleArray(n) { k -> m[k] }

    for (i in 0 until n) {
        if (i % 2 == 0) {
            // Even: N = i/2 + 1
            val N = i / 2 + 1
            // H = hankel(m(1:N), m(N:2*N-1))  [MATLAB 1-based]
            // 0-based: col = m[0..N-1], row = m[N-1..2*N-2]
            val col = DoubleArray(N) { k -> mArr[k] }
            val row = DoubleArray(N) { k -> mArr[N - 1 + k] }
            val H = hankel(col, row)
            hm[i] = H.det()
        } else {
            // Odd: N = (i+1)/2 + 1
            val N = (i + 1) / 2 + 1
            // H = hankel([1, m(1:N-1)], m(N-1:2*N-2))  [MATLAB 1-based]
            // 0-based: col = [1, m[0], m[1], ..., m[N-2]], row = m[N-2..2*N-3]
            val col = DoubleArray(N)
            col[0] = 1.0
            for (k in 1 until N) {
                col[k] = mArr[k - 1]
            }
            val row = DoubleArray(N) { k -> mArr[N - 2 + k] }
            val H = hankel(col, row)
            hm[i] = H.det()
        }
    }

    return hm
}
