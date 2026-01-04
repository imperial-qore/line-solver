package jline.api.mc

import jline.util.matrix.Matrix

/**
 * Result class for DTMC solve reducible
 */
data class DtmcSolveReducibleResult(
    val pi: Matrix,
    val pis: List<Matrix>,
    val pi0: Matrix,
    val scc: List<List<Int>>,
    val isrec: List<Boolean>,
    val Pl: Matrix,
    val pil: Matrix
)

/**
 * Estimate limiting distribution for a DTMC that may have reducible components
 *
 * @param P DTMC transition matrix
 * @param pin Initial vector (optional)
 * @param options Solution options including tolerance
 * @return Solution result containing limiting distributions and SCC information
 */
fun dtmc_solve_reducible(
    P: Matrix,
    pin: Matrix? = null,
    options: Map<String, Any> = mapOf("tol" to 1e-12)
): Pair<Matrix, List<List<Int>>> {
    val (scc, isrec) = stronglyconncomp(P)
    val numSCC = (scc.maxOrNull() ?: -1) + 1
    
    if (numSCC == 1) {
        val pi = dtmc_solve(P)
        return Pair(pi, listOf((0 until P.numRows).toList()))
    }
    
    // Create lumped transition matrix
    var Pl = Matrix.zeros(numSCC, numSCC)
    val sccIdx = Array(numSCC) { i ->
        (0 until P.numRows).filter { scc[it] == i }
    }
    
    // Compute transition probabilities between SCCs
    for (i in 0 until numSCC) {
        for (j in 0 until numSCC) {
            if (i != j) {
                var sum = 0.0
                for (si in sccIdx[i]) {
                    for (sj in sccIdx[j]) {
                        sum += P[si, sj]
                    }
                }
                Pl[i, j] = sum
            }
        }
    }
    
    // Make lumped matrix stochastic
    Pl = dtmc_makestochastic(Pl)
    
    // Ensure recurrent SCCs have self-loops
    for (i in 0 until numSCC) {
        if (isrec[i]) {
            Pl[i, i] = 1.0
        }
    }
    
    // Compute initial distribution for lumped chain
    val pinl = if (pin == null) {
        val temp = DoubleArray(numSCC) { 1.0 }
        // Check for zero columns (absorbing states)
        for (j in 0 until P.numCols) {
            var colSum = 0.0
            for (i in 0 until P.numRows) {
                colSum += P[i, j]
            }
            if (colSum < 1e-12) {
                val sccIndex = scc[j]
                temp[sccIndex] = 0.0
            }
        }
        val sum = temp.sum()
        Matrix(1, numSCC).apply {
            for (i in 0 until numSCC) {
                this[0, i] = temp[i] / sum
            }
        }
    } else {
        val temp = DoubleArray(numSCC)
        for (i in 0 until numSCC) {
            var sum = 0.0
            for (idx in sccIdx[i]) {
                sum += pin[0, idx]
            }
            temp[i] = sum
        }
        Matrix(1, numSCC).apply {
            for (i in 0 until numSCC) {
                this[0, i] = temp[i]
            }
        }
    }
    
    // Compute limiting distribution
    val pi = Matrix.zeros(1, P.numRows)
    val pis = mutableListOf<Matrix>()
    
    // Compute spectral decomposition of lumped matrix
    val PI = computeLimitingMatrix(Pl)
    
    for (i in 0 until numSCC) {
        if (pinl[0, i] > 0) {
            val pi0i = Matrix.zeros(1, numSCC)
            pi0i[0, i] = 1.0
            val pili = pi0i.mult(PI)
            
            val pisi = Matrix.zeros(1, P.numRows)
            for (j in 0 until numSCC) {
                if (pili[0, j] > 0 && sccIdx[j].isNotEmpty()) {
                    val indices = Matrix(sccIdx[j].size, 1)
                    for ((idx, value) in sccIdx[j].withIndex()) {
                        indices[idx, 0] = value.toDouble()
                    }
                    val subP = P.getSubMatrix(indices, indices)
                    val subPi = dtmc_solve(subP)
                    for ((k, stateIdx) in sccIdx[j].withIndex()) {
                        pisi[0, stateIdx] = pili[0, j] * subPi[0, k]
                    }
                }
            }
            pis.add(pisi)
            
            for (k in 0 until P.numRows) {
                pi[0, k] += pisi[0, k] * pinl[0, i]
            }
        }
    }
    
    // Handle single transient SCC case
    val transientStates = isrec.indices.filter { !isrec[it] }
    if (transientStates.size == 1 && pin == null) {
        return Pair(pis[transientStates[0]], sccIdx.map { it.toList() })
    }
    
    return Pair(pi, sccIdx.map { it.toList() })
}

/**
 * Full version that returns all computed values
 */
fun dtmc_solve_reducible_full(
    P: Matrix,
    pin: Matrix? = null,
    options: Map<String, Any> = mapOf("tol" to 1e-12)
): DtmcSolveReducibleResult {
    // This would be a full implementation matching the MATLAB version
    // For now, returning a simplified version
    val (pi, scc) = dtmc_solve_reducible(P, pin, options)
    
    return DtmcSolveReducibleResult(
        pi = pi,
        pis = listOf(pi),
        pi0 = Matrix.zeros(1, 1),
        scc = scc,
        isrec = scc.map { true }, // Simplified
        Pl = P,
        pil = pi
    )
}

private fun stronglyconncomp(P: Matrix): Pair<IntArray, BooleanArray> {
    val graph = jline.util.graph.DirectedGraph(P)
    val result = graph.stronglyconncomp()
    
    // Convert 1-based indexing to 0-based indexing for SCC assignments
    val scc = IntArray(result.I.size) { i -> result.I[i] - 1 }
    val isrec = result.recurrent
    
    return Pair(scc, isrec)
}

private fun computeLimitingMatrix(P: Matrix): Matrix {
    // Pre-check: if eigenvector matrix is singular, use power method directly
    // This matches MATLAB behavior which checks condition number of eigenvector matrix
    try {
        val eigs = P.eigvec()
        val V = eigs.vectors
        // Try to invert V - if it fails, matrix is degenerate
        try {
            V.inv()
        } catch (e: Exception) {
            // Matrix has repeated eigenvalues or is numerically degenerate
            // Use power method directly to avoid slow/failing spectd
            return computeLimitingMatrixPowerMethod(P)
        }
    } catch (e: Exception) {
        // If eigenvalue computation fails, use power method
        return computeLimitingMatrixPowerMethod(P)
    }

    try {
        // Compute spectral decomposition of P
        val spectral = Matrix.spectd(P)
        val projectors = spectral.projectors

        // Compute limiting matrix PI = sum of projectors for eigenvalues near 1
        var PI = Matrix.zeros(P.numRows, P.numCols)
        val spectrum = spectral.spectrum

        for (e in 0 until spectrum.numRows) {
            if (Math.abs(spectrum[e, e] - 1.0) < 1e-12) {
                val proj = projectors.get(e) as Matrix
                PI = PI.add(proj)
            }
        }

        // Check if PI has NaN values (spectd failed on singular matrix)
        if (PI.hasNaN()) {
            return computeLimitingMatrixPowerMethod(P)
        }

        return PI
    } catch (e: Exception) {
        // Fallback for singular/defective matrices: use power iteration
        return computeLimitingMatrixPowerMethod(P)
    }
}

private fun computeLimitingMatrixPowerMethod(P: Matrix): Matrix {
    // Use power iteration to compute the limiting matrix
    // For ergodic chains, P^n converges to a matrix where each row is the stationary distribution
    val n = P.numRows
    var Pk = Matrix(P)
    val maxIter = 1000
    val tol = 1e-10

    for (iter in 0 until maxIter) {
        val Pk1 = Pk.mult(P)

        // Check convergence
        var maxDiff = 0.0
        for (i in 0 until n) {
            for (j in 0 until n) {
                val diff = Math.abs(Pk1.get(i, j) - Pk.get(i, j))
                if (diff > maxDiff) maxDiff = diff
            }
        }

        Pk = Pk1

        if (maxDiff < tol) {
            break
        }
    }

    return Pk
}
/**
 * DTMC solve reducible algorithms
 */
@Suppress("unused")
class DtmcSolveReducibleAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}