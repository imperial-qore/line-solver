/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 *
 * Reference:
 * G Horváth, M Telek, "A minimal representation of Markov arrival processes
 * and a moments matching method," Performance Evaluation 64:(9-12) pp. 1153-1168. (2007)
 */
package jline.lib.butools.reptrans

import jline.util.matrix.Matrix

/**
 * Obtains a Markovian representation from a non-Markovian one while keeping
 * the size the same, by applying a series of elementary transformations.
 *
 * @param rep The initial non-Markovian representation (list of matrices)
 * @param transfun A function that transforms the representation using the given similarity transformation matrix
 * @param evalfun A function that returns how far the representation is from the Markovian one
 * @param precision A representation is considered Markovian if it is closer than the precision (default 1e-7)
 * @return The Markovian representation, if found. If not found, the closest one is returned.
 *
 * Note: This function should not be called directly. It is used by 'phFromME', 'mapFromRAP', etc.
 */
fun findMarkovianRepresentation(
    rep: List<Matrix>,
    transfun: (List<Matrix>, Matrix) -> List<Matrix>,
    evalfun: (List<Matrix>, Int) -> Double,
    precision: Double = 1e-7
): List<Matrix> {
    if (evalfun(rep, 0) < precision) {
        return rep
    }

    var nrep = rep.map { it.copy() }
    val M = nrep[0].numCols
    var b = 0.5
    var odist = Double.MAX_VALUE

    while (b > precision / 2) {
        for (m in 0 until M * M) {
            for (k in 0 until 4) {
                val (newRep, ddist) = minimize(nrep, M * M, b, k, evalfun, transfun)
                nrep = newRep
                if (ddist < precision) {
                    return nrep
                }
            }
            if (odist <= evalfun(nrep, 0)) {
                break
            }
            odist = evalfun(nrep, 0)
        }
        b /= 2.0
    }

    return nrep
}

/**
 * Minimization helper function.
 */
private fun minimize(
    orep: List<Matrix>,
    iters: Int,
    b: Double,
    k: Int,
    evalfun: (List<Matrix>, Int) -> Double,
    transfun: (List<Matrix>, Matrix) -> List<Matrix>
): Pair<List<Matrix>, Double> {
    var lastdist = evalfun(orep, k)
    var bestrep = orep.map { it.copy() }
    var currentRep = orep.map { it.copy() }

    for (i in 0 until iters) {
        val (newRep, dist) = elementary(currentRep, b, k, evalfun, transfun)
        if (dist >= lastdist) {
            break
        } else {
            lastdist = dist
            bestrep = newRep.map { it.copy() }
            currentRep = newRep
        }
    }

    return Pair(bestrep, lastdist)
}

/**
 * Elementary transformation helper function.
 */
private fun elementary(
    erep: List<Matrix>,
    b: Double,
    k: Int,
    evalfun: (List<Matrix>, Int) -> Double,
    transfun: (List<Matrix>, Matrix) -> List<Matrix>
): Pair<List<Matrix>, Double> {
    var bestdist = evalfun(erep, k)
    var bestrep = erep.map { it.copy() }
    val repSize = erep[0].numCols

    for (i in 0 until repSize) {
        for (j in 0 until repSize) {
            if (i != j) {
                // Create elementary transformation matrix with +b
                val Bpos = Matrix.eye(repSize)
                Bpos[i, j] = b
                Bpos[i, i] = 1.0 - b

                // Apply similarity transform
                val newrepPos = transfun(erep, Bpos)
                val newdistPos = evalfun(newrepPos, k)

                // Store result if better
                if (newdistPos < bestdist) {
                    bestrep = newrepPos.map { it.copy() }
                    bestdist = newdistPos
                }

                // Create elementary transformation matrix with -b
                val Bneg = Matrix.eye(repSize)
                Bneg[i, j] = -b
                Bneg[i, i] = 1.0 + b

                // Apply similarity transform
                val newrepNeg = transfun(erep, Bneg)
                val newdistNeg = evalfun(newrepNeg, k)

                // Store result if better
                if (newdistNeg < bestdist) {
                    bestrep = newrepNeg.map { it.copy() }
                    bestdist = newdistNeg
                }
            }
        }
    }

    return Pair(bestrep, bestdist)
}
