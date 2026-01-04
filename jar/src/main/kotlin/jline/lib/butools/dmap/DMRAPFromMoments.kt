/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dmap

import jline.lib.butools.factorialMomsFromMoms
import jline.lib.butools.jFactorialMomsFromJMoms
import jline.lib.butools.dph.mgFromMoments
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Creates a discrete marked rational arrival process that
 * has the same marginal and lag-1 joint moments as given.
 *
 * @param moms The list of marginal moments. To obtain a discrete
 *             marked rational process of order M, 2*M-1 marginal
 *             moments are required.
 * @param Nm The list of lag-1 joint moment matrices. The length
 *           determines K, the number of arrival types.
 * @return The H0, H1, ..., HK matrices of the discrete marked rational process
 */
fun dmrapFromMoments(moms: DoubleArray, Nm: MatrixCell): MatrixCell {
    val mgResult = mgFromMoments(moms)
    val v = mgResult.alpha
    val H0 = mgResult.A

    val N = H0.numRows
    val I = Matrix.eye(N)

    val H0i = I.sub(H0).inv()
    val Ge = Matrix(N, N)
    val G1 = Matrix(N, N)

    var H0ip = Matrix.eye(N)
    for (i in 0 until N) {
        // Ge(i,:) = v * H0ip
        val row = v.mult(H0ip)
        for (j in 0 until N) {
            Ge[i, j] = row[0, j]
        }

        // G1(:,i) = sum(H0ip, 2)
        for (j in 0 until N) {
            var sum = 0.0
            for (k in 0 until N) {
                sum += H0ip[j, k]
            }
            G1[j, i] = sum
        }

        // H0ip = H0ip * (i+1) * H0i
        H0ip = H0ip.scale((i + 1).toDouble()).mult(H0i)
        if (i > 0) {
            H0ip = H0ip.mult(H0)
        }
    }

    val Gei = Ge.inv()
    val G1i = G1.inv()

    val numTypes = Nm.size()
    val H = MatrixCell(numTypes + 1)
    H[0] = H0

    for (i in 0 until numTypes) {
        val Nmi = Nm[i]

        // Extract submatrices and convert moments
        // row1 = FactorialMomsFromMoms(Nmi(1,2:end))
        val row1Input = Matrix(1, N - 1)
        for (j in 0 until N - 1) {
            row1Input[0, j] = Nmi[0, j + 1]
        }
        val row1 = factorialMomsFromMoms(row1Input)

        // col1 = FactorialMomsFromMoms(Nmi(2:end,1))
        val col1Input = Matrix(N - 1, 1)
        for (j in 0 until N - 1) {
            col1Input[j, 0] = Nmi[j + 1, 0]
        }
        val col1 = factorialMomsFromMoms(col1Input)

        // mid = JFactorialMomsFromJMoms(Nmi(2:end,2:end))
        val midInput = Matrix(N - 1, N - 1)
        for (r in 0 until N - 1) {
            for (c in 0 until N - 1) {
                midInput[r, c] = Nmi[r + 1, c + 1]
            }
        }
        val mid = jFactorialMomsFromJMoms(midInput)

        // Construct the transformed Nmi
        val NmiTransformed = Matrix(N, N)
        NmiTransformed[0, 0] = Nmi[0, 0]
        for (j in 0 until N - 1) {
            NmiTransformed[0, j + 1] = row1[0, j]
            NmiTransformed[j + 1, 0] = col1[j, 0]
        }
        for (r in 0 until N - 1) {
            for (c in 0 until N - 1) {
                NmiTransformed[r + 1, c + 1] = mid[r, c]
            }
        }

        // H{i} = (eye(N)-H0)*Gei*Nmi*G1i
        H[i + 1] = I.sub(H0).mult(Gei).mult(NmiTransformed).mult(G1i)
    }

    return H
}

/**
 * Overload for Array<Matrix>.
 */
fun dmrapFromMoments(moms: DoubleArray, Nm: Array<Matrix>): MatrixCell {
    val cell = MatrixCell(Nm.size)
    for (i in Nm.indices) {
        cell[i] = Nm[i]
    }
    return dmrapFromMoments(moms, cell)
}
