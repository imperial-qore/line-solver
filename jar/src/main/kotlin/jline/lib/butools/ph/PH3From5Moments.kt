/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 *
 * Reference:
 * G. Horvath and M. Telek, "On the canonical representation of phase type
 * distributions," Performance Evaluation, vol. 66, no. 8, pp. 396-409, 2009.
 */
package jline.lib.butools.ph

import jline.lib.butools.ReducedMomsFromMoms
import jline.util.Maths
import jline.util.matrix.Matrix
import kotlin.math.abs
import kotlin.math.max
import kotlin.math.pow
import kotlin.math.sqrt

/**
 * Result class for PH3From5Moments.
 */
data class PH3Representation(val alpha: Matrix, val A: Matrix)

/**
 * Returns a PH(3) which has the same 5 moments as given.
 *
 * @param moms The moments to match (length 5).
 * @param prec Numerical precision, default value is 1e-10.
 * @return The PH3Representation containing alpha (initial vector) and A (generator).
 * @throws IllegalArgumentException if the moments are not feasible with a PH(3).
 */
fun ph3From5Moments(moms: DoubleArray, prec: Double = 1e-10): PH3Representation {
    // Convert the moments to reduced moments
    val rmoms = ReducedMomsFromMoms(moms)
    for (i in 0 until 5) {
        rmoms[i] = rmoms[i] / moms[0].pow(i + 1)
    }

    // Solve linear system of equations for a0 a1 a2
    // M = [rmoms(3) -rmoms(2) rmoms(1); rmoms(4) -rmoms(3) rmoms(2); rmoms(5) -rmoms(4) rmoms(3)]
    val M = Matrix(3, 3)
    M[0, 0] = rmoms[2]; M[0, 1] = -rmoms[1]; M[0, 2] = rmoms[0]
    M[1, 0] = rmoms[3]; M[1, 1] = -rmoms[2]; M[1, 2] = rmoms[1]
    M[2, 0] = rmoms[4]; M[2, 1] = -rmoms[3]; M[2, 2] = rmoms[2]

    val b = Matrix(3, 1)
    b[0, 0] = 1.0
    b[1, 0] = rmoms[0]
    b[2, 0] = rmoms[1]

    val a = M.inv().mult(b)
    val a0 = a[0, 0]
    val a1 = a[1, 0]
    val a2 = a[2, 0]

    val discr = a2 * a2 - 3 * a1
    if (discr < 0) {
        throw IllegalArgumentException("Invalid characteristic polynomial!")
    }

    val gu = (a2 + 2 * sqrt(discr)) / 3
    val g0 = (a2 + sqrt(discr)) / 3

    // Find roots of characteristic polynomial: s^3 + a2*s^2 + a1*s + a0
    val coeffs = doubleArrayOf(a0, a1, a2, 1.0)
    val roots = Maths.roots(coeffs)

    // Sort by real part ascending
    val sortedRoots = roots.sortedBy { it.real }
    val lambda = sortedRoots.map { it.negate() }

    val d1 = a1 - a2 - a0 * rmoms[1]
    val d2 = a0 - a1 - a2 * d1
    val d3 = -a0 - a1 * d1 - a2 * d2

    if (d1 > prec || (abs(d1) < prec && d2 > 0)) {
        throw IllegalArgumentException("Negative density around 0!")
    }

    if (lambda[2].real < 0) {
        throw IllegalArgumentException("Invalid eigenvalues!")
    }

    val gl = if (abs(lambda[0].imaginary) < prec) {
        lambda[0].real
    } else {
        g0
    }

    if (gl > gu + prec) {
        throw IllegalArgumentException("Invalid eigenvalues (gl>gu detected)!")
    }
    val glAdj = if (gl > gu) gu else gl

    val g2 = if (abs(d1) < prec) {
        0.0
    } else {
        -d2 / d1
    }

    if (g2 > gu + prec) {
        throw IllegalArgumentException("alpha_2 is negative!")
    }
    val g2Adj = if (g2 > gu) gu else g2

    val x1 = max(g2Adj, glAdj)

    val x13 = if (lambda[0].real == lambda[0].real && abs(lambda[0].imaginary) < prec && g2Adj < glAdj) {
        0.0
    } else {
        x1 - a0 / (x1 * x1 - a2 * x1 + a1)
    }

    var bels = (a2 - x1).pow(2) - 4 * (x1 * x1 - a2 * x1 + a1)
    if (bels < 0 && bels > -prec) {
        bels = 0.0
    }

    val x2 = (a2 - x1 + sqrt(bels)) / 2
    val x3 = (a2 - x1 - sqrt(bels)) / 2
    val p1 = d1 / (x13 - x1)
    val p2 = (x1 * d1 + d2) / (x13 - x1) / x2
    val p3 = (x1 * x2 * d1 + x2 * d2 + x1 * d2 + d3) / (x13 - x1) / x2 / x3

    val A = Matrix(3, 3)
    A[0, 0] = -x1 / moms[0]; A[0, 1] = 0.0; A[0, 2] = x13 / moms[0]
    A[1, 0] = x2 / moms[0]; A[1, 1] = -x2 / moms[0]; A[1, 2] = 0.0
    A[2, 0] = 0.0; A[2, 1] = x3 / moms[0]; A[2, 2] = -x3 / moms[0]

    val alpha = Matrix(1, 3)
    alpha[0, 0] = p1
    alpha[0, 1] = p2
    alpha[0, 2] = p3

    if (x13 < -prec || x13 > x1) {
        throw IllegalArgumentException("Invalid generator!")
    }

    if (alpha.toArray1D().minOrNull()!! < -prec) {
        throw IllegalArgumentException("Initial vector has negative entries!")
    }

    if (alpha.toArray1D().maxOrNull()!! > 1 + prec) {
        throw IllegalArgumentException("Initial vector has entries that are greater than 1!")
    }

    return PH3Representation(alpha, A)
}
