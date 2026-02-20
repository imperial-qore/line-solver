package jline.lib.butools

import jline.util.matrix.Matrix

fun poly(x: Matrix): Matrix {
    val n = x.length()
    val c = Matrix.concatColumns(Matrix.singleton(1.0), Matrix(1, n, n), null)
    for (j in 0..<n) {
        for (i in j + 1 downTo 1) {
            c[i] = c[i] - x[j] * c[i - 1]
        }
    }
    return c
}