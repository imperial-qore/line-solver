/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dph

import jline.lib.butools.ph.canonicalFromPH3
import jline.util.matrix.Matrix
import org.apache.commons.math3.complex.Complex
import kotlin.math.abs

/**
 * Result class for CanonicalFromDPH3 containing both beta and B.
 */
data class DPH3Representation(val beta: Matrix, val B: Matrix)

/**
 * Returns the canonical form of an order-3 discrete phase-type distribution.
 *
 * @param alpha Initial vector of the discrete phase-type distribution
 * @param A Transition probability matrix of the discrete phase-type distribution
 * @param prec Numerical precision for checking the input, default value is 1e-14
 * @return The DPH3Representation containing beta (canonical initial vector) and B (canonical transition matrix)
 */
fun canonicalFromDPH3(alpha: Matrix, A: Matrix, prec: Double = 1e-14): DPH3Representation {
    if (A.numRows != 3 || A.numCols != 3) {
        throw IllegalArgumentException("CanonicalFromDPH3: Dimension must be 3!")
    }

    if (!checkMGRepresentation(alpha, A, prec)) {
        throw IllegalArgumentException("CanonicalFromDPH3: Input isn't a valid MG distribution!")
    }

    // Get eigenvalues sorted by absolute value (descending)
    val eigenvalues: List<Complex> = A.eig()
    val lambda = eigenvalues.sortedByDescending { ev: Complex -> abs(ev.real * ev.real + ev.imaginary * ev.imaginary) }

    // Compute characteristic polynomial coefficients
    val a0 = -lambda[0].real * lambda[1].real * lambda[2].real
    val a1 = lambda[0].real * lambda[1].real + lambda[0].real * lambda[2].real + lambda[1].real * lambda[2].real
    val a2 = -lambda[0].real - lambda[1].real - lambda[2].real

    val N = A.numRows
    val e = Matrix(N, 1)
    for (i in 0 until N) {
        e[i, 0] = 1.0
    }

    var alphaOut: Matrix
    var Aout: Matrix

    if (lambda[0].real > 0 && lambda[1].real >= 0 && lambda[2].real >= 0) {
        // PPP case: Use continuous PH3 canonical form
        val I = Matrix.eye(3)
        val ph3Rep = canonicalFromPH3(alpha, A.sub(I), prec)
        alphaOut = ph3Rep.alpha
        Aout = ph3Rep.A.add(I)
    } else if (lambda[0].real > 0 && lambda[1].real >= 0 && lambda[2].real < 0) {
        // PPN case
        val x1 = lambda[0].real
        val x2 = lambda[1].real + lambda[2].real
        val x3 = lambda[1].real * lambda[2].real / (lambda[1].real + lambda[2].real - 1)

        Aout = Matrix(3, 3)
        Aout[0, 0] = x1
        Aout[0, 1] = 1 - x1
        Aout[0, 2] = 0.0
        Aout[1, 0] = 0.0
        Aout[1, 1] = x2
        Aout[1, 2] = 1 - x2
        Aout[2, 0] = 0.0
        Aout[2, 1] = x3
        Aout[2, 2] = 0.0

        // b3 = (e - A*e) / (1 - x3)
        val Ae = A.mult(e)
        val eMinusAe = e.sub(Ae)
        val b3 = eMinusAe.scale(1.0 / (1 - x3))

        // b2 = A * b3 / (1 - x2)
        val b2 = A.mult(b3).scale(1.0 / (1 - x2))

        // b1 = e - b2 - b3
        val b1 = e.sub(b2).sub(b3)

        // B = [b1, b2, b3]
        val B = Matrix(3, 3)
        for (i in 0 until 3) {
            B[i, 0] = b1[i, 0]
            B[i, 1] = b2[i, 0]
            B[i, 2] = b3[i, 0]
        }

        alphaOut = alpha.mult(B)
    } else if (lambda[0].real > 0 && lambda[1].real < 0 && lambda[2].real >= 0) {
        // PNP case
        val x1 = -a2
        val x2 = (a0 - a1 * a2) / (a2 * (1 + a2))
        val x3 = a0 * (1 + a2) / (a0 - a2 - a1 * a2 - a2 * a2)

        Aout = Matrix(3, 3)
        Aout[0, 0] = x1
        Aout[0, 1] = 1 - x1
        Aout[0, 2] = 0.0
        Aout[1, 0] = x2
        Aout[1, 1] = 0.0
        Aout[1, 2] = 1 - x2
        Aout[2, 0] = 0.0
        Aout[2, 1] = x3
        Aout[2, 2] = 0.0

        val Ae = A.mult(e)
        val eMinusAe = e.sub(Ae)
        val b3 = eMinusAe.scale(1.0 / (1 - x3))
        val b2 = A.mult(b3).scale(1.0 / (1 - x2))
        val b1 = e.sub(b2).sub(b3)

        // Check if first element of alpha*b1 is non-negative
        val alphab1 = alpha.mult(b1)[0, 0]

        if (alphab1 >= 0) {
            val B = Matrix(3, 3)
            for (i in 0 until 3) {
                B[i, 0] = b1[i, 0]
                B[i, 1] = b2[i, 0]
                B[i, 2] = b3[i, 0]
            }
            alphaOut = alpha.mult(B)
        } else {
            // Try to find valid parameters using firstInitElem
            var foundValid = false
            var x33 = 0.0
            var validX1 = 0.0
            var validX2 = 0.0
            var validX3 = 0.0
            var validB = Matrix(3, 3)
            var validA1 = 0.0

            while (x33 <= 1 && !foundValid) {
                val result = firstInitElem(x33, lambda.map { ev: Complex -> ev.real }.toDoubleArray(), alpha, A)
                validA1 = result.a1
                validX1 = result.m1
                validX2 = result.m2
                validX3 = result.m3
                validB = result.B

                if (validA1 >= 0 && validX1 >= 0 && validX2 >= 0 && validX3 >= 0 && validX3 + x33 < 1) {
                    foundValid = true
                    break
                }
                x33 += 0.01
            }

            if (foundValid) {
                Aout[0, 0] = validX1
                Aout[0, 1] = 1 - validX1
                Aout[0, 2] = 0.0
                Aout[1, 0] = validX2
                Aout[1, 1] = 0.0
                Aout[1, 2] = 1 - validX2
                Aout[2, 0] = 0.0
                Aout[2, 1] = validX3
                Aout[2, 2] = x33
                alphaOut = alpha.mult(validB)
            } else {
                // PNP+ case
                val px1 = lambda[2].real
                val px2 = lambda[0].real + lambda[1].real
                val px3 = lambda[0].real * lambda[1].real / (lambda[0].real + lambda[1].real - 1)

                Aout[0, 0] = px1
                Aout[0, 1] = 0.0
                Aout[0, 2] = 0.0
                Aout[1, 0] = 0.0
                Aout[1, 1] = px2
                Aout[1, 2] = 1 - px2
                Aout[2, 0] = 0.0
                Aout[2, 1] = px3
                Aout[2, 2] = 0.0

                val Ae2 = A.mult(e)
                val eMinusAe2 = e.sub(Ae2)
                val p1 = alpha.mult(eMinusAe2)[0, 0]
                val p2 = alpha.mult(A).mult(eMinusAe2)[0, 0]

                val l1 = lambda[0].real
                val l2 = lambda[1].real
                val l3 = lambda[2].real

                val d1 = (1 - l1) * ((1 - l2) * (1 - l3) + (-1 + l2 + l3) * p1 - p2) / ((l1 - l2) * (l1 - l3))
                val d2 = (l2 - 1) * ((1 - l1) * (1 - l3) + (-1 + l1 + l3) * p1 - p2) / ((l1 - l2) * (l2 - l3))
                val d3 = (l3 - 1) * ((1 - l1) * (1 - l2) + (-1 + l1 + l2) * p1 - p2) / ((l2 - l3) * (l3 - l1))

                alphaOut = Matrix(1, 3)
                alphaOut[0, 0] = d3 / (1 - l3)
                alphaOut[0, 1] = (d1 * l1 + d2 * l2) / ((1 - l1) * (1 - l2))
                alphaOut[0, 2] = (d1 + d2) * (1 - l1 - l2) / ((1 - l1) * (1 - l2))

                if (alphaOut.elementMin() < 0 || Aout.elementMin() < 0) {
                    throw IllegalArgumentException("CanonicalFromDPH3: Unhandled PNP case!")
                }
            }
        }
    } else if (lambda[0].real > 0 && lambda[1].real < 0 && lambda[2].real < 0) {
        // PNN case
        val absLambda2 = abs(lambda[1].real * lambda[1].real + lambda[1].imaginary * lambda[1].imaginary)
        if (lambda.all { ev: Complex -> abs(ev.imaginary) < prec } || absLambda2 <= 2 * lambda[0].real * (-lambda[1].real)) {
            val x1 = -a2
            val x2 = -a1 / (1 + a2)
            val x3 = -a0 / (1 + a1 + a2)

            Aout = Matrix(3, 3)
            Aout[0, 0] = x1
            Aout[0, 1] = 1 - x1
            Aout[0, 2] = 0.0
            Aout[1, 0] = x2
            Aout[1, 1] = 0.0
            Aout[1, 2] = 1 - x2
            Aout[2, 0] = x3
            Aout[2, 1] = 0.0
            Aout[2, 2] = 0.0

            val Ae = A.mult(e)
            val eMinusAe = e.sub(Ae)
            val b3 = eMinusAe.scale(1.0 / (1 - x3))
            val b2 = A.mult(b3).scale(1.0 / (1 - x2))
            val b1 = e.sub(b2).sub(b3)

            val B = Matrix(3, 3)
            for (i in 0 until 3) {
                B[i, 0] = b1[i, 0]
                B[i, 1] = b2[i, 0]
                B[i, 2] = b3[i, 0]
            }
            alphaOut = alpha.mult(B)
        } else {
            val I = Matrix.eye(3)
            val ph3Rep = canonicalFromPH3(alpha, A.sub(I), prec)
            alphaOut = ph3Rep.alpha
            Aout = ph3Rep.A.add(I)
        }
    } else {
        throw IllegalArgumentException("CanonicalFromDPH3: Unhandled eigenvalue configuration!")
    }

    return DPH3Representation(alphaOut, Aout)
}

/**
 * Helper function for PNP case.
 */
private data class FirstInitElemResult(val a1: Double, val m1: Double, val m2: Double, val m3: Double, val B: Matrix)

private fun firstInitElem(m33: Double, sortEigs: DoubleArray, alpha: Matrix, A: Matrix): FirstInitElemResult {
    val l1 = sortEigs[0]
    val l2 = sortEigs[1]
    val l3 = sortEigs[2]

    val m1 = -m33 + l1 + l2 + l3

    val m2Num = -((l2 - l3) * (l1 * l1 - l1 * l2 - l1 * l3 + l2 * l3) *
        (m33 * m33 * m33 - 2 * m33 * m33 * l1 + m33 * l1 * l1 - 2 * m33 * m33 * l2 + 3 * m33 * l1 * l2 - l1 * l1 * l2 +
            m33 * l2 * l2 - l1 * l2 * l2 - 2 * m33 * m33 * l3 + 3 * m33 * l1 * l3 - l1 * l1 * l3 + 3 * m33 * l2 * l3 -
            2 * l1 * l2 * l3 - l2 * l2 * l3 + m33 * l3 * l3 - l1 * l3 * l3 - l2 * l3 * l3))

    val m2Den = 2 * m33 * l1 * l1 * l2 + 2 * m33 * m33 * l1 * l1 * l2 - l1 * l1 * l1 * l2 - 3 * m33 * l1 * l1 * l1 * l2 +
        l1 * l1 * l1 * l1 * l2 - 2 * m33 * l1 * l2 * l2 - 2 * m33 * m33 * l1 * l2 * l2 + l1 * l1 * l1 * l2 * l2 +
        l1 * l2 * l2 * l2 + 3 * m33 * l1 * l2 * l2 * l2 - l1 * l1 * l2 * l2 * l2 - l1 * l2 * l2 * l2 * l2 -
        2 * m33 * l1 * l1 * l3 - 2 * m33 * m33 * l1 * l1 * l3 + l1 * l1 * l1 * l3 + 3 * m33 * l1 * l1 * l1 * l3 -
        l1 * l1 * l1 * l1 * l3 + 2 * m33 * l2 * l2 * l3 + 2 * m33 * m33 * l2 * l2 * l3 - l2 * l2 * l2 * l3 -
        3 * m33 * l2 * l2 * l2 * l3 + l2 * l2 * l2 * l2 * l3 + 2 * m33 * l1 * l3 * l3 + 2 * m33 * m33 * l1 * l3 * l3 -
        l1 * l1 * l1 * l3 * l3 - 2 * m33 * l2 * l3 * l3 - 2 * m33 * m33 * l2 * l3 * l3 + l2 * l2 * l2 * l3 * l3 -
        l1 * l3 * l3 * l3 - 3 * m33 * l1 * l3 * l3 * l3 + l1 * l1 * l3 * l3 * l3 + l2 * l3 * l3 * l3 +
        3 * m33 * l2 * l3 * l3 * l3 - l2 * l2 * l3 * l3 * l3 + l1 * l3 * l3 * l3 * l3 - l2 * l3 * l3 * l3 * l3

    val m2 = if (abs(m2Den) > 1e-14) m2Num / m2Den else 0.0

    val m3Num = (l2 - l3) * (l1 * l1 - l1 * l2 - l1 * l3 + l2 * l3) *
        (m33 * m33 * m33 + m33 * m33 * m33 * m33 - m33 * m33 * l1 - 2 * m33 * m33 * m33 * l1 + m33 * m33 * l1 * l1 -
            m33 * m33 * l2 - 2 * m33 * m33 * m33 * l2 + m33 * l1 * l2 + 3 * m33 * m33 * l1 * l2 - m33 * l1 * l1 * l2 +
            m33 * m33 * l2 * l2 - m33 * l1 * l2 * l2 - m33 * m33 * l3 - 2 * m33 * m33 * m33 * l3 + m33 * l1 * l3 +
            3 * m33 * m33 * l1 * l3 - m33 * l1 * l1 * l3 + m33 * l2 * l3 + 3 * m33 * m33 * l2 * l3 - l1 * l2 * l3 -
            4 * m33 * l1 * l2 * l3 + l1 * l1 * l2 * l3 - m33 * l2 * l2 * l3 + l1 * l2 * l2 * l3 + m33 * m33 * l3 * l3 -
            m33 * l1 * l3 * l3 - m33 * l2 * l3 * l3 + l1 * l2 * l3 * l3)

    val m3Den = -2 * m33 * l1 * l1 * l2 - 2 * m33 * m33 * l1 * l1 * l2 - m33 * m33 * m33 * l1 * l1 * l2 +
        l1 * l1 * l1 * l2 + 3 * m33 * l1 * l1 * l1 * l2 + 2 * m33 * m33 * l1 * l1 * l1 * l2 - l1 * l1 * l1 * l1 * l2 -
        m33 * l1 * l1 * l1 * l1 * l2 + 2 * m33 * l1 * l2 * l2 + 2 * m33 * m33 * l1 * l2 * l2 + m33 * m33 * m33 * l1 * l2 * l2 -
        l1 * l1 * l1 * l2 * l2 - 2 * m33 * l1 * l1 * l1 * l2 * l2 + l1 * l1 * l1 * l1 * l2 * l2 - l1 * l2 * l2 * l2 -
        3 * m33 * l1 * l2 * l2 * l2 - 2 * m33 * m33 * l1 * l2 * l2 * l2 + l1 * l1 * l2 * l2 * l2 +
        2 * m33 * l1 * l1 * l2 * l2 * l2 + l1 * l2 * l2 * l2 * l2 + m33 * l1 * l2 * l2 * l2 * l2 - l1 * l1 * l2 * l2 * l2 * l2 +
        2 * m33 * l1 * l1 * l3 + 2 * m33 * m33 * l1 * l1 * l3 + m33 * m33 * m33 * l1 * l1 * l3 - l1 * l1 * l1 * l3 -
        3 * m33 * l1 * l1 * l1 * l3 - 2 * m33 * m33 * l1 * l1 * l1 * l3 + l1 * l1 * l1 * l1 * l3 + m33 * l1 * l1 * l1 * l1 * l3 -
        2 * m33 * l2 * l2 * l3 - 2 * m33 * m33 * l2 * l2 * l3 - m33 * m33 * m33 * l2 * l2 * l3 + l2 * l2 * l2 * l3 +
        3 * m33 * l2 * l2 * l2 * l3 + 2 * m33 * m33 * l2 * l2 * l2 * l3 - l2 * l2 * l2 * l2 * l3 - m33 * l2 * l2 * l2 * l2 * l3 -
        2 * m33 * l1 * l3 * l3 - 2 * m33 * m33 * l1 * l3 * l3 - m33 * m33 * m33 * l1 * l3 * l3 + l1 * l1 * l1 * l3 * l3 +
        2 * m33 * l1 * l1 * l1 * l3 * l3 - l1 * l1 * l1 * l1 * l3 * l3 + 2 * m33 * l2 * l3 * l3 + 2 * m33 * m33 * l2 * l3 * l3 +
        m33 * m33 * m33 * l2 * l3 * l3 - l2 * l2 * l2 * l3 * l3 - 2 * m33 * l2 * l2 * l2 * l3 * l3 + l2 * l2 * l2 * l2 * l3 * l3 +
        l1 * l3 * l3 * l3 + 3 * m33 * l1 * l3 * l3 * l3 + 2 * m33 * m33 * l1 * l3 * l3 * l3 - l1 * l1 * l3 * l3 * l3 -
        2 * m33 * l1 * l1 * l3 * l3 * l3 - l2 * l3 * l3 * l3 - 3 * m33 * l2 * l3 * l3 * l3 - 2 * m33 * m33 * l2 * l3 * l3 * l3 +
        l2 * l2 * l3 * l3 * l3 + 2 * m33 * l2 * l2 * l3 * l3 * l3 - l1 * l3 * l3 * l3 * l3 - m33 * l1 * l3 * l3 * l3 * l3 +
        l1 * l1 * l3 * l3 * l3 * l3 + l2 * l3 * l3 * l3 * l3 + m33 * l2 * l3 * l3 * l3 * l3 - l2 * l2 * l3 * l3 * l3 * l3

    val m3 = if (abs(m3Den) > 1e-14) m3Num / m3Den else 0.0

    // Compute B matrix
    val N = A.numRows
    val e = Matrix(N, 1)
    for (i in 0 until N) {
        e[i, 0] = 1.0
    }
    val I = Matrix.eye(N)

    val eMinusAe = I.sub(A).mult(e)
    val b3 = eMinusAe.scale(1.0 / (1 - m3 - m33))
    val m33I = I.scale(m33)
    val b2 = Matrix.negative(A.sub(m33I)).mult(b3).scale(1.0 / (1 - m2))
    val b1 = A.mult(b2).sub(b3.scale(m3)).scale(1.0 / (1 - m1))

    val B = Matrix(3, 3)
    for (i in 0 until 3) {
        B[i, 0] = b1[i, 0]
        B[i, 1] = b2[i, 0]
        B[i, 2] = b3[i, 0]
    }

    val a1 = alpha.mult(b1)[0, 0]

    return FirstInitElemResult(a1, m1, m2, m3, B)
}

/**
 * Overload for DoubleArray alpha.
 */
fun canonicalFromDPH3(alpha: DoubleArray, A: Matrix, prec: Double = 1e-14): DPH3Representation {
    return canonicalFromDPH3(Matrix(alpha), A, prec)
}
