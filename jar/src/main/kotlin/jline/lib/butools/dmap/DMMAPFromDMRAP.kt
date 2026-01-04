/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dmap

import jline.lib.butools.reptrans.findMarkovianRepresentation
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Obtains a Markovian representation of a discrete rational arrival process
 * of the same size, if possible.
 *
 * @param H The H0...HN matrices of the DMRAP to transform (as MatrixCell)
 * @param prec A representation is considered to be Markovian if it is closer than this precision
 * @return The D0...DN matrices of the DMMAP (if found)
 */
fun dmmapFromDMRAP(H: MatrixCell, prec: Double = 1e-14): MatrixCell {
    if (!checkDMRAPRepresentation(H, prec)) {
        throw IllegalArgumentException("DMMAPFromDMRAP: Input isn't a valid DMRAP representation!")
    }

    // Convert MatrixCell to List<Matrix>
    val HList = (0 until H.size()).map { H[it] }

    // Transformation function: nH = inv(B) * oH * B for each matrix
    val transfun = { oH: List<Matrix>, B: Matrix ->
        val Binv = B.inv()
        oH.map { Binv.mult(it).mult(B) }
    }

    // Evaluation function: measures distance from Markovian representation
    val evalfun = { oH: List<Matrix>, k: Int ->
        val ones = Matrix.ones(oH[0].numRows, oH[0].numCols)
        if (k % 2 == 0) {
            var dist = Double.POSITIVE_INFINITY
            for (i in oH.indices) {
                val minOH = oH[i].elementMin()
                val minOnesMinusOH = ones.sub(oH[i]).elementMin()
                dist = minOf(dist, minOf(minOH, minOnesMinusOH))
            }
            -dist
        } else {
            var dist = 0.0
            for (i in oH.indices) {
                var sumNegOH = 0.0
                var sumNegOnesMinusOH = 0.0
                for (r in 0 until oH[i].numRows) {
                    for (c in 0 until oH[i].numCols) {
                        val v = oH[i][r, c]
                        if (v < 0) sumNegOH += v
                        val vOnes = 1.0 - v
                        if (vOnes < 0) sumNegOnesMinusOH += vOnes
                    }
                }
                dist += minOf(sumNegOH, sumNegOnesMinusOH)
            }
            -dist
        }
    }

    val result = findMarkovianRepresentation(HList, transfun, evalfun, prec)

    // Convert back to MatrixCell
    val resultCell = MatrixCell(result.size)
    for (i in result.indices) {
        resultCell[i] = result[i]
    }
    return resultCell
}

/**
 * Overload for Array<Matrix>.
 */
fun dmmapFromDMRAP(H: Array<Matrix>, prec: Double = 1e-14): MatrixCell {
    val cell = MatrixCell(H.size)
    for (i in H.indices) {
        cell[i] = H[i]
    }
    return dmmapFromDMRAP(cell, prec)
}
