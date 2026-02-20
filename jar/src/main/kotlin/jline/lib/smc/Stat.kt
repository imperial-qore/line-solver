package jline.lib.smc

import jline.util.matrix.Matrix

fun stat(A: Matrix): Matrix {
    val S = A.numRows
    val e = Matrix.ones(S, 1)
    val B = Matrix.concatColumns(A.add(-1.0, Matrix.eye(S)), e, null)
    val y = Matrix.concatColumns(Matrix(1, S), Matrix.singleton(1.0), null)
    // MATLAB: theta = y / B solves X * B = y for X (row vector)
    return y.rightMatrixDivide(B)
}