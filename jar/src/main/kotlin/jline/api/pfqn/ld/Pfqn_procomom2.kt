/**
 * @file Probabilistic class-oriented method of moments for two-station load-dependent models
 * 
 * Implements the probabilistic class-oriented method of moments for analyzing
 * load-dependent queueing models consisting of a single queueing station and delay station.
 * Provides exact marginal state probability computation for specialized load-dependent topologies.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.ld

import jline.io.Ret
import jline.util.Maths
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.CombinatoricsUtils
import org.apache.commons.math3.util.FastMath

/**
 * Compute marginal state probabilities for the queue in a model consisting of a
 * queueing station and a delay station only.
 *
 * @param L  - demands at all stations
 * @param N  - number of jobs for each class
 * @param Z  - think times for each class
 * @param mu - load-dependent scalings (optional)
 * @param m  - number of servers (optional, default 1)
 * @return marginal state probabilities and related matrices
 */
fun pfqn_procomom2(L: Matrix, N: Matrix, Z: Matrix, mu: Matrix?, m: Int?): Ret.pfqnProcomom2 {
    val M = L.numRows
    val R = L.numCols
    val Nsum = N.elementSum().toInt()
    val mServers = m ?: 1
    
    // Initialize mu if not provided or reshape if provided
    val muMatrix: Matrix
    if (mu == null || mu.isEmpty) {
        muMatrix = Matrix.ones(mServers, Nsum + 1)
    } else {
        // Reshape mu as [1, mu(:)']
        val muFlat = Matrix(1, mu.length() + 1)
        muFlat[0, 0] = 1.0
        for (i in 0 until mu.length()) {
            muFlat[0, i + 1] = mu[i]
        }
        muMatrix = muFlat
    }
    
    // Compute solution for [1,0,0,...,0]
    val p0 = Matrix(Nsum + 1, 1)
    p0[Nsum, 0] = 1.0
    
    // Generate T matrices for each class
    val T = ArrayList<Matrix>()
    for (r in 0 until R) {
        val Tr = Matrix(1 + Nsum, 1 + Nsum)
        for (n in Nsum downTo 1) {
            val row = Nsum - n
            Tr.set(row, row, Z[r])
            Tr.set(row, row + 1, (n + mServers - 1) * L[r] / muMatrix[0, n])
        }
        Tr.set(Nsum, Nsum, Z[r])
        T.add(Tr)
    }
    
    // Initialize F and B as identity matrices
    var F = Matrix.eye(Nsum + 1)
    var B = Matrix.eye(Nsum + 1)
    
    // Compute F and B matrices
    for (r in 0 until R) {
        val Nr = N[r].toInt()
        // F = F * T{r}^N(r) / factorial(N(r))
        val Tpower = Matrix.pow(T[r], Nr)
        F = F.mult(Tpower).scale(1.0 / CombinatoricsUtils.factorial(Nr))
        // B = B * T{r}
        B = B.mult(T[r])
    }
    
    // Compute pk = (F*p0)'
    var pk = F.mult(p0).transpose()
    
    // Compute G
    val G = pk.elementSum()
    
    // Compute lG
    val lG: Double
    if (!pk[0, 0].isFinite() || !G.isFinite()) {
        // Use logsumexp for numerical stability
        val logValues = Matrix(1, pk.numCols)
        for (i in 0 until pk.numCols) {
            logValues[0, i] = FastMath.log(pk[0, i])
        }
        lG = Matrix.logsumexp(logValues)
    } else {
        lG = FastMath.log(G)
    }
    
    // Normalize pk
    pk = pk.scale(1.0 / G)
    
    // Reverse the order of pk
    val pkReversed = Matrix(pk.numRows, pk.numCols)
    for (i in 0 until pk.numCols) {
        pkReversed[0, i] = pk[0, pk.numCols - 1 - i]
    }
    
    return Ret.pfqnProcomom2(pkReversed, lG, G, T, F, B)
}

/**
 * Overloaded version with default values
 */
fun pfqn_procomom2(L: Matrix, N: Matrix, Z: Matrix, mu: Matrix?): Ret.pfqnProcomom2 {
    return pfqn_procomom2(L, N, Z, mu, null)
}

/**
 * Overloaded version with default values
 */
fun pfqn_procomom2(L: Matrix, N: Matrix, Z: Matrix): Ret.pfqnProcomom2 {
    return pfqn_procomom2(L, N, Z, null, null)
}
/**
 * PFQN procomom2 algorithms
 */
@Suppress("unused")
class PfqnProcomom2Algo {
    companion object {
        // Class documentation marker for Dokka
    }
}