package jline.lib.butools

import jline.util.matrix.Matrix
import java.util.stream.IntStream

fun SimilarityMatrixForVectors(vecA: Matrix, vecB: Matrix): Matrix {
    val m = vecA.length()
    val neg_vecA = vecA.copy()
    neg_vecA.scaleEq(-1.0)
    val neg_vecA_array = neg_vecA.toArray1D()
    val ix = IntStream.range(0, neg_vecA_array.size).boxed()
        .sorted(Comparator.comparingDouble { i: Int -> neg_vecA_array[i] }).mapToInt { it }.toArray()
    val P = Matrix(m, m, m)
    for (i in 0 until m) {
        P[i, ix[i]] = 1.0
    }
    val cp = P.mult(vecA)

    val B = Matrix(m, m, m * m)
    for (i in 0 until m) {
        var cp_sum = 0.0
        for (j in 0..i) {
            cp_sum += cp[j, 0]
        }
        for (j in 0..i) {
            B[i, j] = vecB[i] / cp_sum
        }
    }
    return B.mult(P)
}