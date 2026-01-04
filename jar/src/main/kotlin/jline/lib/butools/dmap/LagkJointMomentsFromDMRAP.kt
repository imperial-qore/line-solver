/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dmap

import jline.lib.butools.MomsFromFactorialMoms
import jline.lib.butools.jMomsFromJFactorialMoms
import jline.lib.butools.mc.drpSolve
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Returns the lag-L joint moments of a discrete marked rational arrival process.
 *
 * @param H The H0...HN matrices of the DMRAP (as MatrixCell)
 * @param K The dimension of the matrix of joint moments to compute.
 *          If K=0, the MxM joint moments will be computed.
 * @param L The lag at which the joint moments are computed. Default is 1.
 * @param prec Numerical precision to check if the input is valid.
 * @return List of matrices containing the lag-L joint moments
 */
fun lagkJointMomentsFromDMRAP(
    H: MatrixCell,
    K: Int = 0,
    L: Int = 1,
    prec: Double = 1e-14
): MatrixCell {
    if (!checkDMRAPRepresentation(H, prec)) {
        throw IllegalArgumentException("LagkJointMomentsFromDMRAP: Input isn't a valid DMRAP representation!")
    }

    val M = H.size() - 1
    val N = H[0].numRows
    val actualK = if (K == 0) N - 1 else K

    // sumH = sum of H[1] ... H[M]
    var sumH = Matrix(N, N)
    for (i in 1..M) {
        sumH = sumH.add(H[i])
    }

    val I = Matrix.eye(N)
    val iH0 = I.sub(H[0]).inv()

    // pi = DRPSolve(iH0*sumH)
    val pi = drpSolve(iH0.mult(sumH))

    // Compute H0 powers
    val H0p = arrayOfNulls<Matrix>(actualK + 1)
    var Pw = Matrix.eye(N)
    H0p[0] = Pw

    Pw = Pw.mult(iH0)
    H0p[1] = Pw

    for (i in 2..actualK) {
        Pw = Pw.scale(i.toDouble()).mult(iH0).mult(H[0])
        H0p[i] = Pw
    }

    // Pl = (iH0*sumH)^(L-1)
    var Pl = Matrix.eye(N)
    val transitionMatrix = iH0.mult(sumH)
    for (i in 0 until L - 1) {
        Pl = Pl.mult(transitionMatrix)
    }

    // Compute joint moments for each type
    val Nm = MatrixCell(M)
    for (m in 0 until M) {
        val Nmm = Matrix(actualK + 1, actualK + 1)
        for (i in 0..actualK) {
            for (j in 0..actualK) {
                // Nmm(i,j) = sum(pi*H0p{i}*iH0*H{m+1}*Pl*H0p{j})
                val temp = pi.mult(H0p[i]!!).mult(iH0).mult(H[m + 1]).mult(Pl).mult(H0p[j]!!)
                Nmm[i, j] = temp.elementSum()
            }
        }

        // Convert factorial moments to moments
        // row1 = MomsFromFactorialMoms(Nmm(1,2:end))
        val row1Input = Matrix(1, actualK)
        for (j in 0 until actualK) {
            row1Input[0, j] = Nmm[0, j + 1]
        }
        val row1 = MomsFromFactorialMoms(row1Input)

        // col1 = MomsFromFactorialMoms(Nmm(2:end,1))
        val col1Input = Matrix(actualK, 1)
        for (j in 0 until actualK) {
            col1Input[j, 0] = Nmm[j + 1, 0]
        }
        val col1 = MomsFromFactorialMoms(col1Input)

        // mid = JMomsFromJFactorialMoms(Nmm(2:end,2:end))
        val midInput = Matrix(actualK, actualK)
        for (r in 0 until actualK) {
            for (c in 0 until actualK) {
                midInput[r, c] = Nmm[r + 1, c + 1]
            }
        }
        val mid = jMomsFromJFactorialMoms(midInput)

        // Construct result matrix
        val result = Matrix(actualK + 1, actualK + 1)
        result[0, 0] = Nmm[0, 0]
        for (j in 0 until actualK) {
            result[0, j + 1] = row1[0, j]
            result[j + 1, 0] = col1[j, 0]
        }
        for (r in 0 until actualK) {
            for (c in 0 until actualK) {
                result[r + 1, c + 1] = mid[r, c]
            }
        }

        Nm[m] = result
    }

    return Nm
}

/**
 * Overload for Array<Matrix>.
 */
fun lagkJointMomentsFromDMRAP(
    H: Array<Matrix>,
    K: Int = 0,
    L: Int = 1,
    prec: Double = 1e-14
): MatrixCell {
    val cell = MatrixCell(H.size)
    for (i in H.indices) {
        cell[i] = H[i]
    }
    return lagkJointMomentsFromDMRAP(cell, K, L, prec)
}
