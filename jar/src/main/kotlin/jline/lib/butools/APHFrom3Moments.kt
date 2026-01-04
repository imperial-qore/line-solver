/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.lib.butools

import jline.GlobalConstants
import jline.lang.processes.APH
import jline.util.Maths
import jline.util.matrix.Matrix
import org.apache.commons.math3.complex.Complex
import kotlin.math.abs
import kotlin.math.pow


fun APHFrom3Moments(moms: DoubleArray, maxSize: Int = 100, prec: Double = 1e-14): APH {
    var maxSize = maxSize
    var prec = prec
    if (maxSize <= 0) {
        maxSize = 100
    }
    if (prec <= 0) {
        prec = 1e-14
    }

    val m1 = moms[0]
    var m2 = moms[1]
    var m3 = moms[2]

    // Detect number of phases needed
    var n = 2
    while (n < maxSize && (APH2ndMomentLowerBound(m1, n) > m2 || APH3rdMomentLowerBound(m1,
            m2,
            n) >= m3 || APH3rdMomentUpperBound(m1, m2, n) <= m3)) {
        n++
    }

    // If PH is too large, adjust moment to bounds
    if (APH2ndMomentLowerBound(m1, n) > m2) {
        m2 = APH2ndMomentLowerBound(m1, n)
    }

    if (APH3rdMomentLowerBound(m1, m2, n) > m3) {
        m3 = APH3rdMomentLowerBound(m1, m2, n)
    }

    if (APH3rdMomentUpperBound(m1, m2, n) < m3) {
        m3 = APH3rdMomentUpperBound(m1, m2, n)
    }

    // Compute normalized moments
    val nmoms = NormMomsFromMoms(doubleArrayOf(m1, m2, m3))
    val n1 = nmoms[0]
    val n2 = nmoms[1]
    val n3 = nmoms[2]


    if (n2 > 2 || n3 < 2 * n2 - 1) {
        val nComplex = Complex(n.toDouble(), 0.0)
        val n2Complex = Complex(n2, 0.0)
        val n3Complex = Complex(n3, 0.0)
        val four = Complex(4.0, 0.0)
        val two = Complex(2.0, 0.0)
        val three = Complex(3.0, 0.0)
        val twelve = Complex(12.0, 0.0)
        val sixteen = Complex(16.0, 0.0)
        val fifteen = Complex(15.0, 0.0)
        val eight = Complex(8.0, 0.0)
        Complex(-1.0, 0.0)
        val one = Complex.ONE

        val numeratorB = two.multiply(four.subtract(nComplex.multiply(three.multiply(n2Complex).subtract(four))))
        val denominatorB = n2Complex.multiply(four.add(nComplex).subtract(nComplex.multiply(n3Complex)))
            .add((nComplex.multiply(n2Complex)).sqrt()
                .multiply((twelve.multiply(n2Complex.multiply(n2Complex)).multiply(nComplex.add(one))
                    .add(sixteen.multiply(n3Complex).multiply(nComplex.add(one)))
                    .add(n2Complex.multiply(nComplex.multiply(n3Complex.subtract(fifteen)).multiply(n3Complex.add(one))
                        .subtract(eight.multiply(n3Complex.add(Complex(3.0, 0.0)))))).sqrt())))

        val b = numeratorB.divide(denominatorB)

        val numeratorA = (b.multiply(n2Complex).subtract(two)).multiply(nComplex.subtract(one)).multiply(b)
        val denominatorA = b.subtract(one).multiply(nComplex)
        val a = numeratorA.divide(denominatorA)

        val p = (b.subtract(one)).divide(a)

        val lambda = (p.multiply(a).add(one)).divide(n1)

        val mu = (nComplex.subtract(one)).multiply(lambda).divide(a)

        // Construct representation
        val alpha = DoubleArray(n)
        alpha[0] = p.real
        alpha[n - 1] = 1.0 - p.real

        val A = Matrix(n, n)
        A[n - 1, n - 1] = -lambda.real
        for (i in 0..<n - 1) {
            A[i, i] = -mu.real
            A[i, i + 1] = mu.real
        }
        return APH(Matrix(alpha).transpose(), A)
    } else {
        val c4 = n2 * (3.0 * n2 - 2.0 * n3) * (n - 1).toDouble().pow(2.0)
        val c3 = 2.0 * n2 * (n3 - 3.0) * (n - 1).toDouble().pow(2.0)
        val c2 = 6.0 * (n - 1.0) * (n - n2)
        val c1 = 4.0 * n * (2.0 - n)
        val c0 = n * (n - 2.0)

        val coefficients = doubleArrayOf(c0, c1, c2, c3, c4)
        val fs = Maths.roots(coefficients)

        for (f in fs) {
            if (f.isNaN || f.isInfinite) {
                continue
            }

            val temp1 = f.multiply(f).multiply(n2).subtract(f.multiply(2.0)).add(2.0)
            val temp2 = temp1.multiply(n - 1)
            val temp3 = temp2.subtract(n.toDouble())

            if (temp3.abs() < prec) {
                continue
            }

            val a = f.subtract(1.0).multiply(2.0).multiply(n - 1).divide(temp3)
            val p = f.subtract(1.0).multiply(a)
            val lambda = a.add(p).divide(n1)
            val mu = Complex((n - 1).toDouble()).divide(n1 - p.divide(lambda).real)

            if (abs(p.imaginary) < GlobalConstants.Zero) if (abs(lambda.imaginary) < GlobalConstants.Zero) {
                if (abs(mu.imaginary) < GlobalConstants.Zero) {
                    if (p.real >= 0 && p.real <= 1 && lambda.real > 0 && mu.real > 0) {
                        val alpha = DoubleArray(n)
                        alpha[0] = p.real
                        alpha[1] = 1.0 - p.real

                        val A = Matrix(n, n)
                        A[0, 0] = -lambda.real
                        A[0, 1] = lambda.real
                        for (j in 1..<n) {
                            A[j, j] = -mu.real
                            if (j < n - 1) {
                                A[j, j + 1] = mu.real
                            }
                        }
                        // Return or further process alpha and A as needed
                        return APH(Matrix(alpha).transpose(), A)
                    }
                }
            }
        }
    }

    throw IllegalArgumentException("No APH found for the given 3 moments!")
}

