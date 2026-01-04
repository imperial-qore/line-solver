package jline.solvers.mam.handlers

import jline.util.matrix.Matrix

fun qna_superpos(lambda: Matrix, a2: Matrix): Double {
    val lambda_finite_idx: MutableList<Int?> = ArrayList<Int?>()
    for (i in 0..<lambda.length()) {
        if (java.lang.Double.isFinite(lambda.get(i))) {
            lambda_finite_idx.add(i)
        }
    }
    val a2_new = Matrix(1, lambda_finite_idx.size, lambda_finite_idx.size)
    for (i in lambda_finite_idx.indices) {
        a2_new.set(i, a2.get(lambda_finite_idx.get(i)!!))
    }
    val lambda_new = Matrix(1, lambda_finite_idx.size, lambda_finite_idx.size)
    for (i in lambda_finite_idx.indices) {
        lambda_new.set(i, lambda.get(lambda_finite_idx.get(i)!!))
    }
    return a2_new.mult(lambda_new.transpose()).get(0) / lambda_new.elementSum()
}