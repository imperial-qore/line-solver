/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.mam

import jline.lib.butools.FluidFundamentalMatrices
import jline.lib.butools.mc.ctmcSolve
import jline.util.matrix.Matrix

/**
 * Result class for FluidSolve containing the parameters of the
 * matrix-exponentially distributed stationary distribution.
 */
data class FluidSolution(
    val mass0: Matrix,  // Stationary probability vector of zero level
    val ini: Matrix,    // Initial vector of stationary density
    val K: Matrix,      // Matrix parameter of stationary density
    val clo: Matrix     // Closing matrix of stationary density
)

/**
 * Returns the parameters of the matrix-exponentially distributed stationary
 * distribution of a canonical Markovian fluid model.
 *
 * The canonical Markov fluid model is defined by the matrix blocks of the
 * generator of the background Markov chain partitioned according to the
 * sign of the associated fluid rates (i.e., there are "+" and "-" states).
 *
 * Using the returned 4 parameters the stationary solution can be obtained as follows:
 * - The probability that the fluid level is zero while being in different states
 *   of the background process is given by vector mass0.
 * - The density that the fluid level is x while being in different states of
 *   the background process is: pi(x) = ini * e^(K*x) * clo
 *
 * @param Fpp The matrix of transition rates between states having positive fluid rates
 * @param Fpm The matrix of transition rates from positive to negative fluid rate states
 * @param Fmp The matrix of transition rates from negative to positive fluid rate states
 * @param Fmm The matrix of transition rates between states having negative fluid rates
 * @param prec Numerical precision for computing the fundamental matrix (default 1e-14)
 * @return FluidSolution containing mass0, ini, K, and clo
 */
fun fluidSolve(
    Fpp: Matrix,
    Fpm: Matrix,
    Fmp: Matrix,
    Fmm: Matrix,
    prec: Double = 1e-14
): FluidSolution {
    val Np = Fpp.numRows
    val Nm = Fmm.numRows

    // Compute fundamental matrices
    val result = FluidFundamentalMatrices(Fpp, Fpm, Fmp, Fmm, prec, null, null)
    val Psi = result["P"]!!  // Note: FluidFundamentalMatrices uses "P" for Psi
    val K = result["K"]!!
    val U = result["U"]!!

    // Compute stationary distribution of U
    var mass0Minus = ctmcSolve(U)

    // Normalize
    val Ki = K.neg().inv()
    val nr = mass0Minus.elementSum() + 2 * mass0Minus.mult(Fmp).mult(Ki).elementSum()
    mass0Minus = mass0Minus.scale(1.0 / nr)

    // Compute ini
    val ini = mass0Minus.mult(Fmp)

    // Compute clo = [I, Psi]
    val clo = Matrix.zeros(Np, Np + Nm)
    for (i in 0 until Np) {
        clo[i, i] = 1.0
        for (j in 0 until Nm) {
            clo[i, Np + j] = Psi[i, j]
        }
    }

    // mass0 = [zeros(1,Np), mass0Minus]
    val mass0 = Matrix.zeros(1, Np + Nm)
    for (i in 0 until Nm) {
        mass0[0, Np + i] = mass0Minus[0, i]
    }

    return FluidSolution(mass0, ini, K, clo)
}
