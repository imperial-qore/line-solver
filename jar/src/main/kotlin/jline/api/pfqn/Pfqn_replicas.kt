/**
 * Station replica consolidation utilities for Product-Form Queueing Networks
 *
 * Provides functions to identify replicated stations (stations with identical demand
 * rows and service rates) and consolidate them into unique stations with multiplicity.
 * This optimization reduces computational cost for networks with symmetric station replicas.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn

import jline.GlobalConstants
import jline.io.Ret
import jline.util.matrix.Matrix

/**
 * Result class for pfqn_unique containing all output matrices and mapping information.
 *
 * @property L_unique Reduced demand matrix (M' x R) with M' <= M unique stations
 * @property mu_unique Reduced load-dependent rates (M' x Ntot), null if mu was empty
 * @property gamma_unique Reduced class-dependent rates (M' x R), null if gamma was empty
 * @property mi Multiplicity vector (1 x M'), mi[j] = count of stations mapping to unique station j
 * @property mapping Mapping vector (1 x M), mapping[i] = unique station index for original station i
 */
data class PfqnUniqueResult(
    val L_unique: Matrix,
    val mu_unique: Matrix?,
    val gamma_unique: Matrix?,
    val mi: Matrix,
    val mapping: IntArray
)

/**
 * Consolidate replicated stations into unique stations with multiplicity.
 *
 * Identifies stations with identical demand rows L(i,:) and (if present)
 * identical load-dependent rates mu(i,:) or class-dependent rates gamma(i,:).
 * Returns reduced matrices with only unique stations plus a multiplicity vector.
 *
 * @param L Service demand matrix (M x R)
 * @param mu Load-dependent rate matrix (M x Ntot), optional - pass null if not used
 * @param gamma Class-dependent service rate matrix (M x R), optional - pass null if not used
 * @return PfqnUniqueResult containing reduced matrices and mapping information
 */
fun pfqn_unique(L: Matrix, mu: Matrix? = null, gamma: Matrix? = null): PfqnUniqueResult {
    val M = L.numRows
    val R = L.numCols
    val tol = GlobalConstants.FineTol

    // Build fingerprint matrix for comparison
    // Concatenate L, mu, and gamma horizontally if they exist
    val fingerprintCols = R + (mu?.numCols ?: 0) + (gamma?.numCols ?: 0)
    val fingerprint = Matrix(M, fingerprintCols)

    for (i in 0 until M) {
        var col = 0
        // Copy L
        for (j in 0 until R) {
            fingerprint[i, col++] = L[i, j]
        }
        // Copy mu if present
        if (mu != null) {
            for (j in 0 until mu.numCols) {
                fingerprint[i, col++] = mu[i, j]
            }
        }
        // Copy gamma if present
        if (gamma != null) {
            for (j in 0 until gamma.numCols) {
                fingerprint[i, col++] = gamma[i, j]
            }
        }
    }

    // Find unique rows using tolerance-based comparison
    val mapping = IntArray(M) { -1 }
    val uniqueIdx = ArrayList<Int>()
    val miList = ArrayList<Int>()

    for (i in 0 until M) {
        if (mapping[i] == -1) {  // not yet assigned
            // This is a new unique station
            uniqueIdx.add(i)
            val groupIdx = uniqueIdx.size - 1
            mapping[i] = groupIdx
            var count = 1

            // Find all stations identical to this one
            for (j in (i + 1) until M) {
                if (mapping[j] == -1) {  // not yet assigned
                    // Compute infinity norm of difference
                    var maxDiff = 0.0
                    for (k in 0 until fingerprintCols) {
                        val diff = kotlin.math.abs(fingerprint[i, k] - fingerprint[j, k])
                        if (diff > maxDiff) maxDiff = diff
                    }
                    if (maxDiff < tol) {
                        mapping[j] = groupIdx
                        count++
                    }
                }
            }
            miList.add(count)
        }
    }

    // Extract unique rows
    val M_unique = uniqueIdx.size
    val L_unique = Matrix(M_unique, R)
    for (i in 0 until M_unique) {
        for (j in 0 until R) {
            L_unique[i, j] = L[uniqueIdx[i], j]
        }
    }

    val mu_unique = if (mu != null) {
        val result = Matrix(M_unique, mu.numCols)
        for (i in 0 until M_unique) {
            for (j in 0 until mu.numCols) {
                result[i, j] = mu[uniqueIdx[i], j]
            }
        }
        result
    } else null

    val gamma_unique = if (gamma != null) {
        val result = Matrix(M_unique, gamma.numCols)
        for (i in 0 until M_unique) {
            for (j in 0 until gamma.numCols) {
                result[i, j] = gamma[uniqueIdx[i], j]
            }
        }
        result
    } else null

    // Create multiplicity vector
    val mi = Matrix(1, M_unique)
    for (i in 0 until M_unique) {
        mi[0, i] = miList[i].toDouble()
    }

    return PfqnUniqueResult(L_unique, mu_unique, gamma_unique, mi, mapping)
}

/**
 * Expand per-station metrics from reduced model to original dimensions.
 *
 * Expands performance metrics computed on a reduced model (with unique stations)
 * back to the original model dimensions by replicating values according to mapping.
 *
 * @param QN Queue lengths from reduced model (M' x R)
 * @param UN Utilizations from reduced model (M' x R)
 * @param CN Cycle times from reduced model (M' x R)
 * @param mapping Mapping vector from pfqn_unique (length M), mapping[i] = unique station index
 * @return Triple of (QN_full, UN_full, CN_full) in original dimensions (M x R)
 */
fun pfqn_expand(QN: Matrix, UN: Matrix, CN: Matrix, mapping: IntArray): Triple<Matrix, Matrix, Matrix> {
    val M_original = mapping.size
    val R = QN.numCols

    val QN_full = Matrix(M_original, R)
    val UN_full = Matrix(M_original, R)
    val CN_full = Matrix(M_original, R)

    for (i in 0 until M_original) {
        val uniqueIdx = mapping[i]
        for (r in 0 until R) {
            QN_full[i, r] = QN[uniqueIdx, r]
            UN_full[i, r] = UN[uniqueIdx, r]
            CN_full[i, r] = CN[uniqueIdx, r]
        }
    }

    return Triple(QN_full, UN_full, CN_full)
}

/**
 * Combine user-provided multiplicity vector with detected replica multiplicity.
 *
 * For each unique station j, sums the mi values of all original stations mapping to it.
 *
 * @param mi User-provided multiplicity vector (1 x M_original)
 * @param mapping Mapping vector from pfqn_unique (length M_original)
 * @param M_unique Number of unique stations
 * @return Combined multiplicity vector (1 x M_unique)
 */
fun pfqn_combine_mi(mi: Matrix, mapping: IntArray, M_unique: Int): Matrix {
    val mi_combined = Matrix(1, M_unique)
    for (i in 0 until mapping.size) {
        mi_combined[0, mapping[i]] = mi_combined[0, mapping[i]] + mi[0, i]
    }
    return mi_combined
}

/**
 * PFQN replicas algorithms
 */
@Suppress("unused")
class PfqnReplicasAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}
