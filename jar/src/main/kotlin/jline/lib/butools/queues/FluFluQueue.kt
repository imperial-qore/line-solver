/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.queues

import jline.lib.butools.mam.GeneralFluidSolution
import jline.lib.butools.mam.generalFluidSolve
import jline.lib.butools.mc.ctmcSolve
import jline.util.matrix.Matrix

/**
 * Result class for FluFluQueue containing ME distribution parameters
 * and computed performance measures.
 */
data class FluFluResult(
    // Fluid level distribution parameters (mass0, ini, K, clo from GeneralFluidSolve)
    val fluidSolution: GeneralFluidSolution?,
    // Sojourn time distribution parameters (from GeneralFluidSolve on transposed problem)
    val sojournSolution: GeneralFluidSolution?,
    // Fluid level moments (if requested)
    val fluidMoments: DoubleArray?,
    // Sojourn time moments (if requested)
    val sojournMoments: DoubleArray?,
    // Mean input rate (lambda)
    val lambda: Double,
    // Mean output rate (mu)
    val mu: Double
)

/**
 * Returns various performance measures of a fluid queue with independent
 * fluid arrival and service processes.
 *
 * Two types of boundary behavior are available:
 * - If srv0stop=false, the output process evolves continuously even if the queue is empty.
 * - If srv0stop=true, the output process slows down if there is less fluid in the queue
 *   than it can serve. If the queue is empty and the fluid input rate is zero, the output
 *   process freezes until fluid arrives.
 *
 * @param Qin Generator of the background Markov chain for the input process (N x N)
 * @param Rin Diagonal matrix of fluid input rates (N x N)
 * @param Qout Generator of the background Markov chain for the output process (M x M)
 * @param Rout Diagonal matrix of fluid output rates (M x M)
 * @param srv0stop If true, service slows down when queue is nearly empty
 * @param numFluidMoments Number of fluid level moments to compute (0 to skip)
 * @param numSojournMoments Number of sojourn time moments to compute (0 to skip)
 * @param prec Numerical precision (default 1e-14)
 * @return FluFluResult containing requested performance measures
 *
 * References:
 * Horvath G, Telek M, "Sojourn times in fluid queues with independent and
 * dependent input and output processes", PERFORMANCE EVALUATION 79: pp. 160-181, 2014.
 */
fun fluFluQueue(
    Qin: Matrix,
    Rin: Matrix,
    Qout: Matrix,
    Rout: Matrix,
    srv0stop: Boolean,
    numFluidMoments: Int = 0,
    numSojournMoments: Int = 0,
    prec: Double = 1e-14
): FluFluResult {
    val Nin = Qin.numRows
    val Nout = Qout.numRows

    val Iin = Matrix.eye(Nin)
    val Iout = Matrix.eye(Nout)

    // Compute mean input and output rates for normalization
    val piIn = ctmcSolve(Qin, prec)
    val piOut = ctmcSolve(Qout, prec)
    val lambda = piIn.mult(Rin).elementSum()
    val mu = piOut.mult(Rout).elementSum()

    var fluidSolution: GeneralFluidSolution? = null
    var sojournSolution: GeneralFluidSolution? = null
    var fluidMoments: DoubleArray? = null
    var sojournMoments: DoubleArray? = null

    // Compute fluid level distribution if needed
    if (numFluidMoments > 0) {
        // Q = kron(Qin, Iout) + kron(Iin, Qout)
        val Q = Qin.kron(Iout).add(1.0, Iin.kron(Qout))

        // R = kron(Rin, Iout) - kron(Iin, Rout)
        val R = Rin.kron(Iout).add(-1.0, Iin.kron(Rout))

        // Q0 depends on srv0stop
        val Q0: Matrix? = if (srv0stop) {
            // Q0 = kron(Qin, Iout) + kron(Rin, pinv(Rout)*Qout)
            val RoutPinv = Rout.pinv()
            Qin.kron(Iout).add(1.0, Rin.kron(RoutPinv.mult(Qout)))
        } else {
            null  // Q0 = Q (regular boundary)
        }

        fluidSolution = generalFluidSolve(Q, R, Q0, prec)

        // Compute fluid level moments: E[X^m] = m! * ini * (-K)^{-(m+1)} * clo * ones
        val ini = fluidSolution.ini
        val K = fluidSolution.K
        val clo = fluidSolution.clo

        val negK = K.neg()
        val invNegK = negK.inv()
        val onesN = Matrix.ones(clo.numCols, 1)

        fluidMoments = DoubleArray(numFluidMoments)
        var invKPower = invNegK.copy()  // (-K)^{-1}

        for (m in 1..numFluidMoments) {
            // (-K)^{-(m+1)} = (-K)^{-m} * (-K)^{-1}
            invKPower = invKPower.mult(invNegK)
            val moment = factorial(m) * ini.mult(invKPower).mult(clo).mult(onesN).elementSum()
            fluidMoments[m - 1] = moment
        }
    }

    // Compute sojourn time distribution if needed
    if (numSojournMoments > 0) {
        // For sojourn time, use different formulation:
        // Rh = kron(Rin, Iout) - kron(Iin, Rout)
        // Qh = kron(Qin, Rout) + kron(Rin, Qout)
        val Rh = Rin.kron(Iout).add(-1.0, Iin.kron(Rout))
        val Qh = Qin.kron(Rout).add(1.0, Rin.kron(Qout))

        sojournSolution = generalFluidSolve(Qh, Rh, null, prec)

        val inih = sojournSolution.ini
        val Kh = sojournSolution.K
        val cloh = sojournSolution.clo

        // kclo depends on srv0stop
        val kclo: Matrix = if (srv0stop) {
            // kclo = cloh * kron(Rin, Rout) / (lambda * mu)
            cloh.mult(Rin.kron(Rout)).scale(1.0 / (lambda * mu))
        } else {
            // kclo = cloh * kron(Rin, Iout) / lambda
            cloh.mult(Rin.kron(Iout)).scale(1.0 / lambda)
        }

        val negKh = Kh.neg()
        val invNegKh = negKh.inv()

        sojournMoments = DoubleArray(numSojournMoments)
        var invKhPower = invNegKh.copy()  // (-Kh)^{-1}

        for (m in 1..numSojournMoments) {
            // (-Kh)^{-(m+1)} = (-Kh)^{-m} * (-Kh)^{-1}
            invKhPower = invKhPower.mult(invNegKh)
            val moment = factorial(m) * inih.mult(invKhPower).mult(kclo).elementSum()
            sojournMoments[m - 1] = moment
        }
    }

    return FluFluResult(
        fluidSolution = fluidSolution,
        sojournSolution = sojournSolution,
        fluidMoments = fluidMoments,
        sojournMoments = sojournMoments,
        lambda = lambda,
        mu = mu
    )
}

/**
 * Computes factorial of n.
 */
private fun factorial(n: Int): Double {
    var result = 1.0
    for (i in 2..n) {
        result *= i
    }
    return result
}
