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

/**
 * Result class for ExtendToMarkovian.
 */
data class MarkovianRepresentation(val beta: Matrix, val B: Matrix)

/**
 * Extends a non-Markovian initial vector to a Markovian one by appending an Erlang tail.
 *
 * Assume we have an existing monocyclic (or acyclic) representation of a matrix-exponential
 * distribution described by matrix A and vector alpha such that A is Markovian but alpha is not.
 * This procedure appends an appropriate Erlang tail to the representation that makes
 * the result Markovian (both the generator matrix and the initial vector parameter),
 * while keeping the distribution the same.
 *
 * @param alpha The (non-Markovian) initial vector
 * @param A The (Markovian) transient generator
 * @param maxSize The procedure stops if more than maxSize new phases are required (default 100)
 * @param precision The initial vector is considered valid if the smallest entry is greater than -precision (default 1e-14)
 * @return The MarkovianRepresentation containing beta (Markovian initial vector) and B (Markovian generator)
 * @throws IllegalArgumentException if no positive representation is found up to the given size
 */
fun extendToMarkovian(alpha: Matrix, A: Matrix, maxSize: Int = 100, precision: Double = 1e-14): MarkovianRepresentation {
    val N = A.numRows

    // Initial value of t0upper via bisection
    var t0lower = 0.0
    var t0upper = 1.0
    var beta = alpha.mult(A.scale(t0upper).expm())

    while (beta.toArray1D().minOrNull()!! < -precision) {
        t0upper *= 2.0
        beta = alpha.mult(A.scale(t0upper).expm())
    }

    // Interval bisection to find t0
    while ((t0upper - t0lower) / (t0upper + t0lower) > precision) {
        val t0 = (t0upper + t0lower) / 2.0
        beta = alpha.mult(A.scale(t0).expm())
        if (beta.toArray1D().minOrNull()!! < -precision) {
            t0lower = t0
        } else {
            t0upper = t0
        }
    }
    var t0 = t0upper

    // Find optimal length and rate parameters of the Erlang tail
    val increment = 1.1
    var bestT0 = -1.0
    var bestLupper = -1

    for (iter in 0 until 100) {
        // Initial value of Lupper
        var Llower = 1
        var Lupper = 1
        beta = inivecWithTail(alpha, A, Lupper, Lupper.toDouble() / t0)

        while (beta.toArray1D().minOrNull()!! < -precision && Lupper < maxSize) {
            Lupper *= 2
            beta = inivecWithTail(alpha, A, Lupper, Lupper.toDouble() / t0)
        }

        val success = beta.toArray1D().minOrNull()!! >= -precision

        if (success) {
            // Interval bisection for L
            while (Lupper - Llower > 1) {
                val L = (Lupper + Llower) / 2
                beta = inivecWithTail(alpha, A, L, L.toDouble() / t0)
                if (beta.toArray1D().minOrNull()!! < -precision) {
                    Llower = L
                } else {
                    Lupper = L
                }
            }
        }

        if (success) {
            if (bestLupper >= 0 && Lupper > bestLupper) {
                // Previous attempt was better, stop
                break
            } else {
                bestLupper = Lupper
                bestT0 = t0
                t0 *= increment
            }
        } else {
            if (bestLupper >= 0) {
                break
            } else {
                t0 *= increment
            }
        }
    }

    if (bestLupper < 0) {
        throw IllegalArgumentException("No positive representation found up to the given size!")
    }

    t0 = bestT0
    val Lupper = bestLupper

    // Final result
    beta = inivecWithTail(alpha, A, Lupper, Lupper.toDouble() / t0)
    val B = addErlangTail(A, Lupper, Lupper.toDouble() / t0)

    return MarkovianRepresentation(beta, B)
}

/**
 * Computes the initial vector with Erlang tail.
 */
private fun inivecWithTail(gamma: Matrix, G: Matrix, tailLength: Int, mu: Double): Matrix {
    val vlen = G.numRows + tailLength
    val beta = Matrix.zeros(1, vlen)

    val WG = Matrix.eye(G.numRows).add(G.scale(1.0 / mu))
    var opv = gamma.copy()

    // Compute closing vector: clv = -sum(G/mu, 2)
    val clv = Matrix(G.numRows, 1)
    for (i in 0 until G.numRows) {
        var sum = 0.0
        for (j in 0 until G.numCols) {
            sum += G[i, j] / mu
        }
        clv[i, 0] = -sum
    }

    for (k in vlen - 1 downTo G.numRows) {
        beta[0, k] = opv.mult(clv)[0, 0]
        opv = opv.mult(WG)
    }

    for (i in 0 until G.numRows) {
        beta[0, i] = opv[0, i]
    }

    return beta
}

/**
 * Adds an Erlang tail to a generator matrix.
 */
private fun addErlangTail(D: Matrix, len: Int, mu: Double): Matrix {
    val DN = D.numRows
    val E = Matrix.zeros(DN + len, DN + len)

    // Copy D to upper-left
    for (i in 0 until DN) {
        for (j in 0 until DN) {
            E[i, j] = D[i, j]
        }
    }

    // Connection from D to Erlang tail
    var rowSum = 0.0
    for (j in 0 until DN) {
        rowSum += D[DN - 1, j]
    }
    E[DN - 1, DN] = -rowSum

    // Erlang tail
    for (ei in 0 until len) {
        E[DN + ei, DN + ei] = -mu
        if (ei < len - 1) {
            E[DN + ei, DN + ei + 1] = mu
        }
    }

    return E
}
