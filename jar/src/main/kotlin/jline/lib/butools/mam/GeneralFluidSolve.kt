/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.mam

import jline.lib.butools.FluidFundamentalMatrices
import jline.util.matrix.Matrix
import kotlin.math.abs

/**
 * Result class for GeneralFluidSolve containing the parameters of the
 * matrix-exponentially distributed stationary distribution.
 */
data class GeneralFluidSolution(
    val mass0: Matrix,  // Stationary probability vector of zero level
    val ini: Matrix,    // Initial vector of stationary density
    val K: Matrix,      // Matrix parameter of stationary density
    val clo: Matrix     // Closing matrix of stationary density
)

/**
 * Returns the parameters of the matrix-exponentially distributed stationary
 * distribution of a general Markovian fluid model, where the fluid rates
 * associated with the states of the background process can be arbitrary
 * (zero is allowed as well).
 *
 * Using the returned 4 parameters the stationary solution can be obtained as follows:
 * - The probability that the fluid level is zero while being in different states
 *   of the background process is given by vector mass0.
 * - The density that the fluid level is x while being in different states of
 *   the background process is: pi(x) = ini * e^(K*x) * clo
 *
 * @param Q The generator of the background Markov chain (N x N)
 * @param R The diagonal matrix of fluid rates associated with different states (N x N)
 * @param Q0 The generator of the background Markov chain at level 0 (optional, defaults to Q)
 * @param prec Numerical precision for computing the fundamental matrix (default 1e-14)
 * @return GeneralFluidSolution containing mass0, ini, K, and clo
 */
fun generalFluidSolve(
    Q: Matrix,
    R: Matrix,
    Q0: Matrix? = null,
    prec: Double = 1e-14
): GeneralFluidSolution {
    val N = Q.numRows

    // Partition the state space according to zero, positive and negative fluid rates
    val ixz = mutableListOf<Int>()  // zero rate states
    val ixp = mutableListOf<Int>()  // positive rate states
    val ixn = mutableListOf<Int>()  // negative rate states

    for (i in 0 until N) {
        val rate = R[i, i]
        when {
            abs(rate) <= prec -> ixz.add(i)
            rate > prec -> ixp.add(i)
            rate < -prec -> ixn.add(i)
        }
    }

    val Nz = ixz.size
    val Np = ixp.size
    val Nn = ixn.size

    // Build permutation matrix P that converts between original and partitioned state ordering
    // Ordering: zero states, positive states, negative states
    val P = Matrix.zeros(N, N)
    for (i in 0 until Nz) {
        P[i, ixz[i]] = 1.0
    }
    for (i in 0 until Np) {
        P[Nz + i, ixp[i]] = 1.0
    }
    for (i in 0 until Nn) {
        P[Nz + Np + i, ixn[i]] = 1.0
    }
    val iP = P.inv()

    // Reorder states with permutation matrix P: 0, +, -
    val Qv = P.mult(Q).mult(iP)
    val Rv = P.mult(R).mult(iP)

    // Extract submatrices for censoring zero-rate states
    // Qv is partitioned as: [Qv_00, Qv_0pm; Qv_pm0, Qv_pmpm]
    val Qv00 = if (Nz > 0) Matrix.getSubMatrix(Qv, 0, Nz, 0, Nz) else Matrix.zeros(0, 0)
    val Qv0pm = if (Nz > 0) Matrix.getSubMatrix(Qv, 0, Nz, Nz, N) else Matrix.zeros(0, Np + Nn)
    val Qvpm0 = if (Nz > 0) Matrix.getSubMatrix(Qv, Nz, N, 0, Nz) else Matrix.zeros(Np + Nn, 0)
    val Qvpmpm = Matrix.getSubMatrix(Qv, Nz, N, Nz, N)

    // New fluid process censored to states + and -
    // Qbar = Qv[+,-;+,-] + Qv[+,-;0] * pinv(-Qv[0;0]) * Qv[0;+,-]
    val Qbar = if (Nz > 0) {
        val negQv00 = Qv00.neg()
        val iQv00 = negQv00.pinv()
        Qvpmpm.add(1.0, Qvpm0.mult(iQv00).mult(Qv0pm))
    } else {
        Qvpmpm.copy()
    }

    // Build diagonal matrix of absolute inverse rates for +/- states
    val absRi = Matrix.zeros(Np + Nn, Np + Nn)
    for (i in 0 until Np + Nn) {
        absRi[i, i] = abs(1.0 / Rv[Nz + i, Nz + i])
    }

    // Normalize by rates
    val Qz = absRi.mult(Qbar)

    // Calculate fundamental matrices on partitioned blocks
    // Qz is partitioned as: [Qz_pp, Qz_pn; Qz_np, Qz_nn]
    val Qzpp = if (Np > 0) Matrix.getSubMatrix(Qz, 0, Np, 0, Np) else Matrix.zeros(0, 0)
    val Qzpn = if (Np > 0 && Nn > 0) Matrix.getSubMatrix(Qz, 0, Np, Np, Np + Nn) else Matrix.zeros(Np, Nn)
    val Qznp = if (Np > 0 && Nn > 0) Matrix.getSubMatrix(Qz, Np, Np + Nn, 0, Np) else Matrix.zeros(Nn, Np)
    val Qznn = if (Nn > 0) Matrix.getSubMatrix(Qz, Np, Np + Nn, Np, Np + Nn) else Matrix.zeros(0, 0)

    val result = FluidFundamentalMatrices(Qzpp, Qzpn, Qznp, Qznn, prec, null, null)
    val Psi = result["P"]!!
    val K = result["K"]!!
    val U = result["U"]!!

    // Closing matrix Pm = [I_Np, Psi]
    val Pm = Matrix.zeros(Np, Np + Nn)
    for (i in 0 until Np) {
        Pm[i, i] = 1.0
        for (j in 0 until Nn) {
            Pm[i, Np + j] = Psi[i, j]
        }
    }

    // absRi partitioned
    val iCn = if (Nn > 0) Matrix.getSubMatrix(absRi, Np, Np + Nn, Np, Np + Nn) else Matrix.zeros(0, 0)
    val iCp = if (Np > 0) Matrix.getSubMatrix(absRi, 0, Np, 0, Np) else Matrix.zeros(0, 0)

    // Submatrices from Qv for closing matrix
    val QvPz = if (Nz > 0) Matrix.getSubMatrix(Qv, Nz, Nz + Np, 0, Nz) else Matrix.zeros(Np, 0)
    val QvNz = if (Nz > 0) Matrix.getSubMatrix(Qv, Nz + Np, N, 0, Nz) else Matrix.zeros(Nn, 0)

    // clo = [(iCp*Qv[+,0] + Psi*iCn*Qv[-,0])*pinv(-Qv[0,0]), Pm*absRi]
    var clo: Matrix
    if (Nz > 0) {
        val negQv00 = Qv00.neg()
        val iQv00 = negQv00.pinv()
        val leftPart = iCp.mult(QvPz).add(1.0, Psi.mult(iCn).mult(QvNz)).mult(iQv00)
        val rightPart = Pm.mult(absRi)
        clo = Matrix.zeros(Np, N)
        for (i in 0 until Np) {
            for (j in 0 until Nz) {
                clo[i, j] = leftPart[i, j]
            }
            for (j in 0 until Np + Nn) {
                clo[i, Nz + j] = rightPart[i, j]
            }
        }
    } else {
        clo = Pm.mult(absRi)
    }

    val mass0: Matrix
    val ini: Matrix

    if (Q0 == null) {
        // Regular boundary behavior
        // Go back to the original state ordering
        clo = clo.mult(P)

        // Calculate boundary vector
        // Ua = iCn*Qv[-,0]*iQv00*ones(Nz,1) + iCn*ones(Nn,1) + Qz[-,+]*inv(-K)*clo*ones(N,1)
        val onesN = Matrix.ones(N, 1)
        val onesNn = Matrix.ones(Nn, 1)
        val onesNz = if (Nz > 0) Matrix.ones(Nz, 1) else Matrix.zeros(0, 1)

        val negK = K.neg()
        val invNegK = negK.inv()

        var Ua: Matrix
        if (Nz > 0) {
            val negQv00 = Qv00.neg()
            val iQv00 = negQv00.pinv()
            val term1 = iCn.mult(QvNz).mult(iQv00).mult(onesNz)
            val term2 = iCn.mult(onesNn)
            val term3 = Qznp.mult(invNegK).mult(clo).mult(onesN)
            Ua = term1.add(1.0, term2).add(1.0, term3)
        } else {
            val term2 = iCn.mult(onesNn)
            val term3 = Qznp.mult(invNegK).mult(clo).mult(onesN)
            Ua = term2.add(1.0, term3)
        }

        // Solve [U, Ua]' * pm' = [0, ..., 0, 1]'
        // pm = linsolve([U, Ua]', [zeros(1,Nn), 1]')'
        val UAugmented = Matrix.zeros(Nn, Nn + 1)
        for (i in 0 until Nn) {
            for (j in 0 until Nn) {
                UAugmented[i, j] = U[i, j]
            }
            UAugmented[i, Nn] = Ua[i, 0]
        }
        val rhs = Matrix.zeros(Nn + 1, 1)
        rhs[Nn, 0] = 1.0

        // Solve UAugmented' * x = rhs, then pm = x'
        val UAugT = UAugmented.transpose()
        val pmT = UAugT.leftMatrixDivide(rhs)
        val pm = pmT.transpose()

        // mass0 = [pm*iCn*Qv[-,0]*iQv00, zeros(1,Np), pm*iCn] * P
        val mass0Unsorted = Matrix.zeros(1, N)
        if (Nz > 0) {
            val negQv00 = Qv00.neg()
            val iQv00 = negQv00.pinv()
            val zeroTerm = pm.mult(iCn).mult(QvNz).mult(iQv00)
            for (j in 0 until Nz) {
                mass0Unsorted.set(0, j, zeroTerm[0, j])
            }
        }
        // zeros(1, Np) is already set
        val negTerm = pm.mult(iCn)
        for (j in 0 until Nn) {
            mass0Unsorted.set(0, Nz + Np + j, negTerm[0, j])
        }
        mass0 = mass0Unsorted.mult(P)

        // ini = pm * Qz[-,+]
        ini = pm.mult(Qznp)
    } else {
        // Custom boundary behavior with Q0
        val Q0v = P.mult(Q0).mult(iP)

        // Build M matrix
        // M = [-clo*Rv; Q0v[-,:]; Q0v[0,:]]
        val M = Matrix.zeros(N, N)

        // First Np rows: -clo*Rv
        val negCloRv = clo.mult(Rv).neg()
        for (i in 0 until Np) {
            for (j in 0 until N) {
                M[i, j] = negCloRv[i, j]
            }
        }

        // Next Nn rows: Q0v[-, :]
        for (i in 0 until Nn) {
            for (j in 0 until N) {
                M[Np + i, j] = Q0v[Nz + Np + i, j]
            }
        }

        // Last Nz rows: Q0v[0, :]
        for (i in 0 until Nz) {
            for (j in 0 until N) {
                M[Np + Nn + i, j] = Q0v[i, j]
            }
        }

        // Ma = [sum(inv(-K)*clo, 2); ones(Nz+Nn, 1)]
        val negK = K.neg()
        val invNegK = negK.inv()
        val invKClo = invNegK.mult(clo)
        val Ma = Matrix.zeros(N, 1)
        for (i in 0 until Np) {
            var rowSum = 0.0
            for (j in 0 until N) {
                rowSum += invKClo[i, j]
            }
            Ma[i, 0] = rowSum
        }
        for (i in Np until N) {
            Ma[i, 0] = 1.0
        }

        // Solve [M, Ma]' * sol' = [0, ..., 0, 1]'
        val MAugmented = Matrix.zeros(N, N + 1)
        for (i in 0 until N) {
            for (j in 0 until N) {
                MAugmented[i, j] = M[i, j]
            }
            MAugmented[i, N] = Ma[i, 0]
        }
        val rhs = Matrix.zeros(N + 1, 1)
        rhs[N, 0] = 1.0

        val MAugT = MAugmented.transpose()
        val solT = MAugT.leftMatrixDivide(rhs)
        val sol = solT.transpose()

        // ini = sol(1:Np)
        ini = Matrix.zeros(1, Np)
        for (j in 0 until Np) {
            ini.set(0, j, sol[0, j])
        }

        // clo = clo * P (go back to original ordering)
        clo = clo.mult(P)

        // mass0 = [sol(Np+Nn+1:end), zeros(1,Np), sol(Np+1:Np+Nn)] * P
        val mass0Unsorted = Matrix.zeros(1, N)
        // Zero states: sol(Np+Nn+1:end) = sol(Np+Nn:N-1) (0-indexed)
        for (j in 0 until Nz) {
            mass0Unsorted.set(0, j, sol[0, Np + Nn + j])
        }
        // Positive states: zeros
        // Negative states: sol(Np+1:Np+Nn) = sol(Np:Np+Nn-1) (0-indexed)
        for (j in 0 until Nn) {
            mass0Unsorted.set(0, Nz + Np + j, sol[0, Np + j])
        }
        mass0 = mass0Unsorted.mult(P)
    }

    return GeneralFluidSolution(mass0, ini, K, clo)
}
