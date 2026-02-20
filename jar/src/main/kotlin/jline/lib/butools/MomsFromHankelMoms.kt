package jline.lib.butools

import jline.util.matrix.Matrix

/**
 * Constructs a Hankel matrix from the given first column and last row.
 */
private fun hankelMH(col: DoubleArray, row: DoubleArray): Matrix {
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
 * Returns the raw moments given the Hankel moments.
 *
 * @param hm The list of Hankel moments (starting with the first moment)
 * @return The list of raw moments
 *
 * Reference: http://en.wikipedia.org/wiki/Stieltjes_moment_problem
 */
fun momsFromHankelMoms(hm: Matrix): Matrix {
    val n = hm.length()
    val mList = mutableListOf(hm[0])

    for (i in 1 until n) {
        val mArr = mList.toDoubleArray()
        val mSize = mArr.size
        val N: Int
        val H: Matrix

        if (i % 2 == 0) {
            N = i / 2 + 1
            // H = hankel(m(1:N), [m(N:2*N-2), 0])  [MATLAB 1-based]
            val col = DoubleArray(N) { k -> mArr[k] }
            val row = DoubleArray(N)
            for (k in 0 until N - 1) {
                row[k] = mArr[N - 1 + k]
            }
            row[N - 1] = 0.0
            H = hankelMH(col, row)
        } else {
            N = (i + 1) / 2 + 1
            // H = hankel([1, m(1:N-1)], [m(N-1:2*N-3), 0])  [MATLAB 1-based]
            val col = DoubleArray(N)
            col[0] = 1.0
            for (k in 1 until N) {
                col[k] = mArr[k - 1]
            }
            val row = DoubleArray(N)
            for (k in 0 until N - 1) {
                row[k] = mArr[N - 2 + k]
            }
            row[N - 1] = 0.0
            H = hankelMH(col, row)
        }

        // Solve for the unknown moment using cofactor expansion along the last row
        var h = hm[i]
        val rH = Matrix.zeros(N - 1, N) // first N-1 rows of H
        for (r in 0 until N - 1) {
            for (c in 0 until N) {
                rH[r, c] = H[r, c]
            }
        }

        var lastCofactor = 0.0
        for (j in 0 until N) {
            // Compute cofactor: (-1)^(N+j-1) * det(rH with column j removed)
            // Build submatrix by removing column j
            val sub = Matrix.zeros(N - 1, N - 1)
            for (r in 0 until N - 1) {
                var cc = 0
                for (c in 0 until N) {
                    if (c != j) {
                        sub[r, cc] = rH[r, c]
                        cc++
                    }
                }
            }
            val sign = if ((N + j - 1) % 2 == 0) 1.0 else -1.0
            val cofactor = sign * sub.det()

            if (j < N - 1) {
                h -= cofactor * H[N - 1, j]
            } else {
                lastCofactor = cofactor
            }
        }

        mList.add(h / lastCofactor)
    }

    val m = Matrix(hm.numRows, hm.numCols, n)
    for (i in 0 until n) {
        m[i] = mList[i]
    }
    return m
}
