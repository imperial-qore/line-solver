/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 *
 * Reference:
 * Mocanu, S., Commault, C.: "Sparse representations of phase-type distributions,"
 * Stoch. Models 15, 759-778 (1999)
 */
package jline.lib.butools.reptrans

import jline.util.matrix.Matrix
import org.apache.commons.math3.complex.Complex
import kotlin.math.*

/**
 * Data class representing a Feedback Erlang Block (FEB).
 */
private data class FEB(
    val lambda: Double,
    val z: Double,
    val n: Int,
    val multip: Int,
    val evals: List<Complex>,
    val emuls: List<Int>
)

/**
 * Transforms an arbitrary matrix to a Markovian monocyclic matrix.
 *
 * @param A Matrix parameter of the initial representation (shape N,N)
 * @param maxSize The maximal order of the resulting Markovian representation (default 100)
 * @param precision Matrix entries smaller than the precision are considered to be zeros (default 1e-14)
 * @return Transient generator matrix of the Markovian monocyclic representation.
 *         Note that M>N if there are complex eigenvalues.
 * @throws IllegalArgumentException if no Markovian monocyclic generator has been found up to the given size.
 */
fun transformToMonocyclic(A: Matrix, maxSize: Int = 100, precision: Double = 1e-14): Matrix {
    // Determine eigenvalues and their multiplicities
    val (evalues, repeats) = eigvalc(A, precision)

    // Build monocyclic representation
    val febs = generateFEBs(evalues, repeats, maxSize, precision)

    // Calculate total size
    var totalSize = 0
    for (feb in febs) {
        totalSize += feb.n * feb.multip
    }

    // Assemble generator matrix of the hyper-feb
    val B = Matrix.zeros(totalSize, totalSize)
    var pos = 0
    for (i in febs.indices) {
        val Ni = febs[i].n * febs[i].multip
        val febGen = febGenerator(febs[i].lambda, febs[i].z, febs[i].n, febs[i].multip)

        for (r in 0 until Ni) {
            for (c in 0 until Ni) {
                B[pos + r, pos + c] = febGen[r, c]
            }
        }

        if (i < febs.size - 1) {
            var rowSum = 0.0
            for (c in 0 until totalSize) {
                rowSum += B[pos + Ni - 1, c]
            }
            B[pos + Ni - 1, pos + Ni] = -rowSum
        }
        pos += Ni
    }

    return B
}

/**
 * Computes eigenvalues and their algebraic multiplicity.
 */
private fun eigvalc(A: Matrix, precision: Double): Pair<List<Complex>, List<Int>> {
    val tol = sqrt(Math.ulp(1.0))  // sqrt(machine epsilon), matches MATLAB's sqrt(eps)

    val eigenvalues = A.eig().toMutableList()
    eigenvalues.sortBy { it.real }

    // Round eigenvalues
    for (i in eigenvalues.indices) {
        val rounded = Complex(
            (eigenvalues[i].real / tol).roundToLong() * tol,
            (eigenvalues[i].imaginary / tol).roundToLong() * tol
        )
        eigenvalues[i] = rounded
    }

    // Get unique eigenvalues and their multiplicities
    val uniqueEvals = mutableListOf<Complex>()
    val repeats = mutableListOf<Int>()

    for (ev in eigenvalues) {
        var found = false
        for (j in uniqueEvals.indices) {
            if (ev.subtract(uniqueEvals[j]).abs() <= tol) {
                repeats[j] = repeats[j] + 1
                found = true
                break
            }
        }
        if (!found) {
            uniqueEvals.add(ev)
            repeats.add(1)
        }
    }

    return Pair(uniqueEvals, repeats)
}

/**
 * Generates Feedback Erlang Blocks from eigenvalues.
 */
private fun generateFEBs(evalues: List<Complex>, repeats: List<Int>, maxSize: Int, precision: Double): List<FEB> {
    val febs = mutableListOf<FEB>()
    var i = 0
    var size = 0

    while (i < evalues.size) {
        val multip = repeats[i]
        val evalimag = -abs(evalues[i].imaginary)

        if (-evalimag < precision) {
            // Real eigenvalue
            val n = 1
            val sigma = -evalues[i].real
            val z = 0.0
            val ev = listOf(Complex(evalues[i].real, 0.0))
            val em = listOf(multip)
            size += 1

            febs.add(FEB(sigma, z, n, multip, ev, em))
            i++
        } else {
            // Complex eigenvalue
            var n = 3
            size += 3
            while (evalimag / evalues[i].real >= 1.0 / tan(PI / n)) {
                n++
                size++
                if (size > maxSize) {
                    throw IllegalArgumentException("The representation is too large (>maxSize). No result returned.")
                }
            }

            val sigma = -(2 * evalues[i].real + evalimag * (1.0 / tan(PI / n) - tan(PI / n))) / 2
            val z = (-evalimag * (1.0 / tan(PI / n) + tan(PI / n)) / (2 * sigma)).pow(n)

            val ev = mutableListOf<Complex>()
            val em = mutableListOf<Int>()
            for (k in 0 until n) {
                val zRoot = z.pow(1.0 / n)
                val angle = 2 * (k) * PI / n
                ev.add(Complex(
                    -(1 - zRoot * cos(angle)) * sigma,
                    zRoot * sin(angle) * sigma
                ))
                em.add(multip)
            }

            febs.add(FEB(sigma, z, n, multip, ev, em))
            i += 2 // Skip conjugate pair
        }
    }

    // Order according to the dominant eigenvalue (descending)
    // MATLAB max() on complex arrays selects by magnitude, then < compares real parts
    return febs.sortedByDescending { feb -> feb.evals.maxByOrNull { it.abs() }?.real ?: Double.MIN_VALUE }
}

/**
 * Generates the generator matrix for a single FEB.
 */
private fun febGenerator(lambda: Double, z: Double, n: Int, multip: Int): Matrix {
    val size = multip * n
    val A = Matrix.zeros(size, size)

    for (ii in 0 until multip) {
        // Fill the n x n block
        for (j in 0 until n) {
            A[ii * n + j, ii * n + j] = -lambda
            if (j < n - 1) {
                A[ii * n + j, ii * n + j + 1] = lambda
            }
        }
        // Feedback connection
        A[ii * n + n - 1, ii * n] = A[ii * n + n - 1, ii * n] + z * lambda

        // Connection to next block
        if (ii < multip - 1) {
            A[ii * n + n - 1, (ii + 1) * n] = (1 - z) * lambda
        }
    }

    return A
}
