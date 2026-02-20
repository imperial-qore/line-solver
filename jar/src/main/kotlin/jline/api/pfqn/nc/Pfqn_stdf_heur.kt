/**
 * @file Heuristic Sojourn Time Distribution Function for multi-server FCFS stations
 * 
 * Implements a heuristic variant of McKenna's 1987 method for computing sojourn time
 * distributions at multi-server FCFS stations. Uses class-dependent MAP approximations
 * and heuristic rate adjustments for enhanced computational efficiency.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.nc

import jline.api.mam.map_cdf
import jline.api.mam.map_exponential
import jline.api.mam.map_sumind
import jline.GlobalConstants
import jline.io.line_warning
import jline.api.pfqn.ld.pfqn_comomrm_ld
import jline.api.pfqn.ld.pfqn_mu_ms
import jline.api.pfqn.ld.pfqn_mushift
import jline.api.pfqn.ld.pfqn_mvald
import jline.solvers.SolverOptions
import jline.util.Maths
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import kotlin.math.*

/**
 * Heuristic sojourn time distribution analysis at multiserver FCFS nodes
 * based on a variant of the method in J. McKenna 1987 JACM
 *
 * @param L - routing matrix
 * @param N - population vector
 * @param Z - think times
 * @param S - number of servers per station
 * @param fcfsNodes - indices of FCFS nodes
 * @param rates - service rates
 * @param tset - time points for evaluation
 * @return response time distributions
 */
fun pfqn_stdf_heur(L: Matrix,
                   N: Matrix,
                   Z: Matrix,
                   S: Matrix,
                   fcfsNodes: Matrix,
                   rates: Matrix,
                   tset: Matrix): Array<Array<Matrix?>> {
    var stabilityWarnIssued = false
    val M = L.numRows
    val R = L.numCols
    val T = tset.length()

    // Initialize service rate matrix
    val mu = Matrix(M, N.elementSum().toInt())
    for (k in 0..<M) {
        for (n in 0..<N.elementSum().toInt()) {
            mu[k, n] = min(S[k], (n + 1).toDouble())
        }
    }

    // Ensure t > 0 for stability of hkc function
    val tsetMod = tset.copy()
    for (i in 0..<tsetMod.length()) {
        if (tsetMod[i] == 0.0) {
            tsetMod[i] = GlobalConstants.FineTol
        }
    }

    val RD = Array<Array<Matrix?>>(M) { Array(R) { null } }
    val hkc = Array<Matrix?>(R) { null }

    for (k in 0..<fcfsNodes.length()) {
        val kIdx = fcfsNodes[k].toInt()
        for (r in 0..<R) {
            if (L[kIdx, r] > GlobalConstants.FineTol) {
                val Nr = Matrix.oner(N, arrayListOf(r))

                val lGr: Double
                val Q1: Matrix
                val isNumStable: Boolean

                if (M == 1) {
                    val options = SolverOptions()
                    val result1 = pfqn_comomrm_ld(L, Nr, Z, mu, options)
                    val m = 2
                    val muMs = pfqn_mu_ms(N.elementSum().toInt(), m, S[kIdx].toInt())
                    val result2 = pfqn_comomrm_ld(L, Nr, Z, muMs, options)

                    val Q1Local = Matrix(M, R)
                    for (s in 0..<R) {
                        Q1Local[kIdx, s] = L[kIdx, s] * exp(result2.lG - result1.lG)
                    }
                    lGr = result1.lG
                    Q1 = Q1Local
                    isNumStable = true
                } else {
                    val result = pfqn_mvald(L, Nr, Z, mu)
                    lGr = result.lG[result.lG.size - 1]
                    Q1 = result.Q
                    isNumStable = result.isNumStable
                }

                // Build heuristic MAPs for each population level
                val hMAPr = Array<MatrixCell?>(N.elementSum().toInt() + 1) { null }
                for (n in 0..N.elementSum().toInt()) {
                    hMAPr[n] = if (n < S[kIdx]) {
                        // Use average rates for n < S(k)
                        map_exponential(1.0 / rates[kIdx, r])
                    } else {
                        // Use rates per class for n >= S(k)
                        val mapList = mutableListOf<MatrixCell>()
                        mapList.add(map_exponential(1.0 / rates[kIdx, r]))

                        for (s in 0..<R) {
                            val qSum = Q1.sumRows(kIdx)
                            if (qSum > GlobalConstants.FineTol) {
                                val rate = Q1[kIdx, s] * (n - S[kIdx] + 1) / qSum / rates[kIdx, s]
                                mapList.add(map_exponential(rate))
                            }
                        }
                        map_sumind(mapList.toTypedArray())
                    }
                }

                // Compute CDFs
                hkc[r] = Matrix(T, N.elementSum().toInt() + 1)
                for (n in 0..N.elementSum().toInt()) {
                    val cdfValues = map_cdf(hMAPr[n]!!, tsetMod)
                    for (t in 0..<T) {
                        hkc[r]!![t, n] = cdfValues[t]
                    }
                }

                if (!isNumStable && !stabilityWarnIssued) {
                    stabilityWarnIssued = true
                    line_warning("pfqn_stdf_heur", "The computation of the sojourn time distribution is numerically unstable.")
                }

                RD[kIdx][r] = Matrix(tset.length(), 2)
                for (t in 0..<tset.length()) {
                    RD[kIdx][r]!![t, 1] = tset[t]
                }

                // Compute response time distribution using recursive form
                val Hkrt = DoubleArray(T)
                val LReduced = Matrix(L.numRows - 1, L.numCols)
                val muReduced = Matrix(mu.numRows - 1, mu.numCols)
                var rowIdx = 0
                for (i in 0..<M) {
                    if (i != kIdx) {
                        for (j in 0..<R) {
                            LReduced[rowIdx, j] = L[i, j]
                        }
                        for (j in 0..<mu.numCols) {
                            muReduced[rowIdx, j] = mu[i, j]
                        }
                        rowIdx++
                    }
                }

                val lGk = if (M == 1) {
                    0.0
                } else {
                    val result = pfqn_mvald(LReduced, Nr, Z, muReduced)
                    result.lG[result.lG.size - 1]
                }

                for (t in 0..<T) {
                    val gammat = mu.copy()
                    for (m in 0..<Nr.elementSum().toInt()) {
                        if (m + 1 < hkc[r]!!.numCols) {
                            val denominator = hkc[r]!![t, m + 1]
                            if (abs(denominator) > GlobalConstants.FineTol) {
                                gammat[kIdx, m] = mu[kIdx, m] * hkc[r]!![t, m] / denominator
                            }
                        }
                    }

                    val gammak = pfqn_mushift(gammat, kIdx)
                    Hkrt[t] = hkc[r]!![t, 0] * exp(lGk)

                    for (s in 0..<R) {
                        if (Nr[s] > 0) {
                            val lYks_t = pfqn_rd(L, Matrix.oner(Nr, arrayListOf(s)), Z, gammak, null).lG
                            if (gammat[kIdx, 0] > GlobalConstants.FineTol) {
                                Hkrt[t] += (L[kIdx, s] * hkc[r]!![t, 0] / gammat[kIdx, 0]) * exp(lYks_t!!)
                            }
                        }
                    }
                }

                // Handle NaN values
                for (t in 0..<T) {
                    if (Hkrt[t].isNaN()) {
                        Hkrt[t] = GlobalConstants.FineTol
                    }
                }

                val lHkrt = Hkrt.map { ln(it) }.toDoubleArray()

                for (t in 0..<T) {
                    RD[kIdx][r]!![t, 0] = exp(lHkrt[t] - lGr)
                }
            }
        }
    }

    return RD
}

/**
 * Helper function Fm(m,x) - not used in current implementation but kept for reference
 */
private fun Fm(m: Int, x: Double): Double {
    return if (m == 1) {
        1.0 - exp(-x)
    } else {
        var A = 0.0
        for (j in 0..<m) {
            A += x.pow(j) / Maths.fact(j)
        }
        1.0 - exp(-x) * A
    }
}
/**
 * PFQN stdf heur algorithms
 */
@Suppress("unused")
class PfqnStdfHeurAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}