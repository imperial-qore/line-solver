/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 *
 * Reference:
 * A. van de Liefvoort. The moment problem for continuous distributions.
 * Technical report, University of Missouri, WP-CM-1990-02, Kansas City, 1990.
 */
package jline.lib.butools.dph

import jline.lib.butools.reducedMomsFromMoms
import jline.lib.butools.factorialMomsFromMoms
import jline.lib.butools.ph.meFromMoments
import jline.util.matrix.Matrix

/**
 * Result class for MGFromMoments containing both alpha and A.
 */
data class MGRepresentation(val alpha: Matrix, val A: Matrix)

/**
 * Creates a matrix-geometric distribution that has the same moments as given.
 *
 * @param moms The list of moments. The order of the resulting matrix-geometric distribution
 *        is determined based on the number of moments given. To obtain a matrix-geometric
 *        distribution of order M, 2*M-1 moments are required.
 * @return The MGRepresentation containing alpha (initial vector) and A (matrix parameter).
 */
fun mgFromMoments(moms: DoubleArray): MGRepresentation {
    // Convert raw moments to factorial moments, then to reduced moments
    val fmoms = factorialMomsFromMoms(Matrix(moms))
    val rfmoms = reducedMomsFromMoms(fmoms)

    // Prepend 1 to rfmoms
    val rfmomsArr = DoubleArray(moms.size + 1)
    rfmomsArr[0] = 1.0
    for (i in 0 until moms.size) {
        rfmomsArr[i + 1] = rfmoms[i]
    }

    // Compute vlist
    val vlist = DoubleArray(moms.size)
    val tmpVec = DoubleArray(moms.size + 1)
    var k = 1.0
    tmpVec[0] = rfmomsArr[0]

    for (i in 1..moms.size) {
        val sign = if (i % 2 == 0) 1.0 else -1.0
        tmpVec[i] = sign * rfmomsArr[i]
        k *= i
        var sum = 0.0
        for (j in 0..i) {
            sum += tmpVec[j]
        }
        vlist[i - 1] = k * sum
    }

    // Use MEFromMoments to get alpha and C
    val meRep = meFromMoments(vlist)
    val alpha = meRep.alpha
    val C = meRep.A

    // A = inv(C) * inv(inv(C) + I)
    val iC = C.inv()
    val I = Matrix.eye(C.numRows)
    val A = iC.mult(iC.add(I).inv())

    return MGRepresentation(alpha, A)
}
