/**
 * @file Sojourn Time Distribution Function (STDF) for multi-server FCFS stations
 * 
 * Implements McKenna's 1987 method for computing sojourn time distributions at
 * multi-server FCFS stations in closed queueing networks. Uses Markovian Arrival Process
 * (MAP) representations and recursive computation for exact distribution analysis.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.nc

import jline.api.mam.map_cdf
import jline.api.mam.map_erlang
import jline.api.mam.map_exponential
import jline.api.mam.map_sumind
import jline.GlobalConstants
import jline.io.line_warning
import jline.api.pfqn.ld.pfqn_comomrm_ld
import jline.api.pfqn.ld.pfqn_mushift
import jline.api.pfqn.ld.pfqn_mvald
import jline.solvers.SolverOptions
import jline.util.Maths
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import kotlin.math.*

/**
 * Sojourn time distribution function at multiserver FCFS nodes
 * J. McKenna 1987 JACM
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
fun pfqn_stdf(L: Matrix,
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

    for (k in 0..<fcfsNodes.length()) {
        val kIdx = fcfsNodes[k].toInt()
        // Check that FCFS station has uniform service rates
        val rateRow = Matrix(1, R)
        Matrix.extract(rates, kIdx, kIdx + 1, 0, R, rateRow, 0, 0)
        val rateRange = rateRow.elementMax() - rateRow.elementMin()
        if (rateRange > GlobalConstants.FineTol) {
            throw IllegalArgumentException("The FCFS stations has distinct service rates, the model is invalid.")
        }

        // Build MAPs for each population level
        val hMAP = Array<MatrixCell?>(N.elementSum().toInt() + 1) { null }
        for (n in 0..N.elementSum().toInt()) {
            hMAP[n] = if (n < S[kIdx]) {
                map_exponential(1.0 / rates[kIdx, 0])
            } else {
                val mapList = arrayOf(map_exponential(1.0 / rates[kIdx, 0]),
                    map_erlang((n - S[kIdx] + 1) / (S[kIdx] * rates[kIdx, 0]), (n - S[kIdx] + 1).toInt()))
                map_sumind(mapList)
            }
        }

        // Compute CDFs
        val hkc = Matrix(T, N.elementSum().toInt() + 1)
        for (n in 0..N.elementSum().toInt()) {
            val cdfValues = map_cdf(hMAP[n]!!, tsetMod)
            for (t in 0..<T) {
                hkc[t, n] = cdfValues[t]
            }
        }

        for (r in 0..<R) {
            if (L[kIdx, r] > GlobalConstants.FineTol) {
                val Nr = Matrix.oner(N, arrayListOf(r))

                val lGr: Double
                val isNumStable: Boolean

                if (L.numRows == 1) {
                    val options = SolverOptions()
                    val result = pfqn_comomrm_ld(L, Nr, Z, mu, options)
                    lGr = result.lG
                    isNumStable = true
                } else {
                    val result = pfqn_mvald(L, Nr, Z, mu)
                    lGr = result.lG[result.lG.size - 1]
                    isNumStable = result.isNumStable
                }

                if (!isNumStable && !stabilityWarnIssued) {
                    stabilityWarnIssued = true
                    line_warning("pfqn_stdf", "The computation of the sojourn time distribution is numerically unstable.")
                }

                RD[kIdx][r] = Matrix(tset.length(), 2)
                for (t in 0..<tset.length()) {
                    RD[kIdx][r]!![t, 1] = tset[t]
                }

                // Compute response time distribution using recursive form
                val Hkrt = DoubleArray(T)

                // Create reduced matrices excluding station k
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

                val lGk = if (L.numRows == 1) {
                    0.0
                } else {
                    val result = pfqn_mvald(LReduced, Nr, Z, muReduced)
                    result.lG[result.lG.size - 1]
                }

                for (t in 0..<T) {
                    val gammat = mu.copy()
                    for (m in 0..<Nr.elementSum().toInt()) {
                        if (m + 1 < hkc.numCols) {
                            val denominator = hkc[t, m + 1]
                            if (abs(denominator) > GlobalConstants.FineTol) {
                                gammat[kIdx, m] = mu[kIdx, m] * hkc[t, m] / denominator
                            }
                        }
                    }

                    val gammak = pfqn_mushift(gammat, kIdx)
                    // Extract reduced gammak matrix (excluding one column)
                    val gammakCols = min(gammak.numCols, Nr.elementSum().toInt() - 1)
                    val gammakReduced = Matrix(gammak.numRows, gammakCols)
                    if (gammakCols > 0) {
                        Matrix.extract(gammak, 0, gammak.numRows, 0, gammakCols, gammakReduced, 0, 0)
                    }

                    Hkrt[t] = hkc[t, 0] * exp(lGk)

                    for (s in 0..<R) {
                        if (Nr[s] > 0) {
                            val lYks_t = if (L.numRows == 1) {
                                val options = SolverOptions()
                                val result =
                                    pfqn_comomrm_ld(L, Matrix.oner(Nr, arrayListOf(s)), Z, gammakReduced, options)
                                result.lG
                            } else {
                                val result = pfqn_mvald(L, Matrix.oner(Nr, arrayListOf(s)), Z, gammakReduced)
                                result.lG[result.lG.size - 1]
                            }

                            if (gammat[kIdx, 0] > GlobalConstants.FineTol) {
                                Hkrt[t] += (L[kIdx, s] * hkc[t, 0] / gammat[kIdx, 0]) * exp(lYks_t)
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
                    RD[kIdx][r]!![t, 0] = min(1.0, exp(lHkrt[t] - lGr))
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
 * PFQN stdf algorithms
 */
@Suppress("unused")
class PfqnStdfAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}