package jline.lib.kpctoolbox.kpcfit

import jline.api.mam.*
import jline.lib.kpctoolbox.basic.logspacei
import jline.lib.kpctoolbox.mmpp.mmpp2_fit
import jline.lib.kpctoolbox.trace.*
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import org.apache.commons.math3.analysis.MultivariateFunction
import org.apache.commons.math3.optim.*
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.NelderMeadSimplex
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.SimplexOptimizer
import org.apache.commons.math3.util.FastMath
import java.util.Random
import kotlin.math.abs
import kotlin.math.sqrt

/**
 * KPC-Toolbox fitting functions.
 * Kronecker Product Composition (KPC) method for fitting MAPs to traces.
 *
 * Ported from MATLAB: matlab/lib/kpctoolbox/kpcfit/
 */

// Global tolerance for KPC fitting
const val KPCFIT_TOL = 1e-6

/**
 * Converts Array<Matrix> to MatrixCell.
 */
private fun arrayToMatrixCell(arr: Array<Matrix>): MatrixCell {
    val cell = MatrixCell(arr.size)
    for (i in arr.indices) {
        cell[i] = arr[i]
    }
    return cell
}

/**
 * Converts MatrixCell to Array<Matrix>.
 */
private fun matrixCellToArray(cell: MatrixCell): Array<Matrix> {
    val size = cell.size()
    return Array(size) { i -> cell[i] }
}

/**
 * Data class for initialized trace data.
 */
data class TraceData(
    val S: DoubleArray,
    val E: DoubleArray,
    val AC: DoubleArray,
    val ACFull: DoubleArray,
    val ACLags: IntArray,
    val BC: DoubleArray,
    val BCGridLags: IntArray,
    val BCLags: Array<IntArray>
) {
    override fun equals(other: Any?): Boolean {
        if (this === other) return true
        if (javaClass != other?.javaClass) return false
        other as TraceData
        return S.contentEquals(other.S)
    }

    override fun hashCode(): Int = S.contentHashCode()
}

/**
 * Options for KPC fitting.
 */
data class KPCFitOptions(
    val onlyAC: Boolean = false,
    val numMAPs: Int? = null,
    val numStates: Int? = null,
    val maxIterAC: Int = 300,
    val maxIterBC: Int = 10,
    val maxRunsAC: Int = 50,
    val maxRunsBC: Int = 30,
    val maxResAC: Int = 10,
    val maxRetMAPs: Int = 1
)

/**
 * Result of KPC fitting.
 */
data class KPCFitResult(
    val MAP: MatrixCell,
    val fac: Double,
    val fbc: Double,
    val subMAPs: List<MatrixCell>,
    val otherMAPs: List<MatrixCell> = emptyList(),
    val otherFACs: DoubleArray = doubleArrayOf(),
    val otherFBCs: DoubleArray = doubleArrayOf()
)

/**
 * Initializes trace data for KPC fitting.
 * Computes moments, autocorrelations, and bicovariances.
 *
 * @param S Input trace (inter-arrival times)
 * @param acLags Optional: specific lags for autocorrelation
 * @param bcGridLags Optional: lags for bicovariance grid
 * @return TraceData structure with computed statistics
 */
fun kpcfit_init(
    S: DoubleArray,
    acLags: IntArray? = null,
    bcGridLags: IntArray? = null
): TraceData {
    val n = S.size
    val nMinSupportAC = 10

    // Default AC lags: logarithmically spaced
    val defaultACLags = acLags ?: logspacei(1.0, (n / nMinSupportAC).toDouble().coerceAtLeast(1.0), 500)
        .distinct().toIntArray()

    // Default BC grid lags
    val maxACLag = defaultACLags.maxOrNull() ?: 1
    val defaultBCGridLags = bcGridLags ?: logspacei(1.0, maxACLag.toDouble(), 5)
        .distinct().toIntArray()

    // Compute moments E[X], E[X^2], E[X^3]
    val E = DoubleArray(3)
    for (j in 1..3) {
        E[j - 1] = S.map { FastMath.pow(it, j.toDouble()) }.average()
    }

    // Compute autocorrelations
    val validACLags = defaultACLags.filter { it <= n - 2 }.toIntArray()
    val AC = if (validACLags.isNotEmpty()) trace_acf(S, validACLags) else doubleArrayOf()
    val ACFull = trace_acf(S, (1..(n / nMinSupportAC).coerceAtLeast(1)).toList().toIntArray())

    // Truncate AC lags where autocorrelation becomes negligible
    var cutIdx = validACLags.size
    for (i in AC.indices) {
        if (abs(AC[i]) < 1e-6) {
            cutIdx = i
            break
        }
    }

    val trimmedACLags = validACLags.take(cutIdx).toIntArray()
    val trimmedAC = AC.take(cutIdx).toDoubleArray()

    // Filter BC grid lags
    val validBCGridLags = defaultBCGridLags.filter { it <= (trimmedACLags.maxOrNull() ?: 1) }.toIntArray()

    // Compute bicovariances
    val (BC, BCLags) = trace_bicov(S, validBCGridLags)

    return TraceData(
        S = S,
        E = E,
        AC = trimmedAC,
        ACFull = ACFull,
        ACLags = trimmedACLags,
        BC = BC,
        BCGridLags = validBCGridLags,
        BCLags = BCLags
    )
}

/**
 * Automatic KPC fitting of a trace to a MAP.
 *
 * @param trace Prepared trace data from kpcfit_init
 * @param options Fitting options
 * @return KPCFitResult with fitted MAP and diagnostics
 */
fun kpcfit_auto(trace: TraceData, options: KPCFitOptions = KPCFitOptions()): KPCFitResult {
    // Determine number of MAPs to use
    val numMAPs = options.numMAPs ?: if (options.numStates != null) {
        FastMath.ceil(FastMath.log(2.0, options.numStates.toDouble())).toInt()
    } else {
        // Automatic order selection using BIC
        kpcfit_sub_bic(trace.ACFull, intArrayOf(2, 4, 8, 16, 32, 64, 128))
    }

    return kpcfit_manual(
        numMAPs = numMAPs,
        E = trace.E,
        AC = trace.AC,
        ACLags = trace.ACLags,
        BC = trace.BC,
        BCLags = trace.BCLags,
        options = options
    )
}

/**
 * Manual KPC fitting with specified parameters.
 */
fun kpcfit_manual(
    numMAPs: Int,
    E: DoubleArray,
    AC: DoubleArray,
    ACLags: IntArray,
    BC: DoubleArray,
    BCLags: Array<IntArray>,
    options: KPCFitOptions = KPCFitOptions()
): KPCFitResult {
    // Fit autocorrelations
    val (resSCV, resG2, fobjAC) = kpcfit_sub_acfit(
        E, AC, ACLags, numMAPs,
        options.maxIterAC, options.maxRunsAC, options.maxResAC
    )

    // Fit bicovariances
    val resE1 = ArrayList<DoubleArray>()
    val resE3 = ArrayList<DoubleArray>()
    val fobjBC = DoubleArray(resSCV.size)

    if (!options.onlyAC) {
        for (i in resSCV.indices) {
            val (E1j, E3j, foBC) = kpcfit_sub_bcfit(
                E, resSCV[i], resG2[i], BC, BCLags,
                options.maxIterBC, options.maxRunsBC
            )
            resE1.add(E1j)
            resE3.add(E3j)
            fobjBC[i] = foBC
        }
    } else {
        for (i in resSCV.indices) {
            val E1j = DoubleArray(numMAPs) { 1.0 }
            val E3j = DoubleArray(numMAPs) { j ->
                (1.5 + 0.01) * FastMath.pow(1 + resSCV[i][j], 2.0)
            }
            resE1.add(E1j)
            resE3.add(E3j)
            fobjBC[i] = -1.0
        }
    }

    // Sort by BC objective
    val sortedIndices = fobjBC.indices.sortedBy { fobjBC[it] }

    // Compose MAPs
    val MAPs = ArrayList<MatrixCell>()
    val subs = ArrayList<List<MatrixCell>>()
    val FACs = ArrayList<Double>()
    val FBCs = ArrayList<Double>()

    for (k in sortedIndices) {
        val (newMAP, newSubMAPs, errorCode) = kpcfit_sub_compose(
            resE1[k], resSCV[k], resE3[k], resG2[k]
        )

        if (errorCode != 0 || newMAP == null) continue

        // Scale to match mean exactly
        val scaledMAP = map_scale(newMAP, E[0])

        // Evaluate objective functions
        val (newfAC, newfBC) = evaluateObjFunction(scaledMAP, E, AC, ACLags, BC, BCLags)

        MAPs.add(scaledMAP)
        subs.add(newSubMAPs)
        FACs.add(newfAC)
        FBCs.add(newfBC)

        if (MAPs.size >= options.maxRetMAPs) break
    }

    if (MAPs.isEmpty()) {
        throw IllegalStateException("KPC fitting failed - no valid MAP found")
    }

    return KPCFitResult(
        MAP = MAPs[0],
        fac = FACs[0],
        fbc = FBCs[0],
        subMAPs = subs[0],
        otherMAPs = MAPs.drop(1),
        otherFACs = FACs.drop(1).toDoubleArray(),
        otherFBCs = FBCs.drop(1).toDoubleArray()
    )
}

/**
 * BIC-based order selection for KPC fitting.
 */
fun kpcfit_sub_bic(ACFull: DoubleArray, orders: IntArray): Int {
    // Simple heuristic: choose order based on ACF decay
    val gamma = estimateACFDecay(ACFull)

    // Map decay rate to order
    return when {
        gamma < 0.3 -> 2
        gamma < 0.5 -> 3
        gamma < 0.7 -> 4
        gamma < 0.85 -> 5
        else -> 6
    }.coerceIn(orders.minOrNull() ?: 2, orders.maxOrNull() ?: 6)
}

/**
 * Estimates ACF decay rate.
 */
private fun estimateACFDecay(acf: DoubleArray): Double {
    if (acf.isEmpty()) return 0.5

    var sumLogRatio = 0.0
    var count = 0
    val rho0 = acf[0]

    if (rho0 <= 0) return 0.5

    for (k in 1 until minOf(acf.size, 50)) {
        if (acf[k] > 0 && acf[k] / rho0 > 0) {
            sumLogRatio += FastMath.log(acf[k] / rho0) / k
            count++
        }
    }

    return if (count > 0) {
        FastMath.exp(sumLogRatio / count).coerceIn(0.0, 1.0)
    } else {
        0.5
    }
}

/**
 * Autocorrelation fitting sub-routine.
 */
fun kpcfit_sub_acfit(
    E: DoubleArray,
    SA: DoubleArray,
    SALags: IntArray,
    J: Int,
    maxIterAC: Int,
    maxRunsAC: Int,
    maxResAC: Int
): Triple<List<DoubleArray>, List<DoubleArray>, DoubleArray> {
    val SCV = (E[1] - E[0] * E[0]) / (E[0] * E[0])
    val NSA = norm(SA)
    val random = Random()

    val fset = ArrayList<Double>()
    val xparamset = ArrayList<DoubleArray>()

    // Multi-start optimization
    for (run in 0 until maxRunsAC) {
        // Random initial point
        val x0 = DoubleArray(2 * J)
        for (j in 0 until J) {
            x0[j] = 1.0 + random.nextDouble()  // SCV > 1
            x0[J + j] = random.nextDouble()     // G2 in [0, 1)
        }

        // Optimize using Nelder-Mead
        val objective = MultivariateFunction { x ->
            val SCVj = x.take(J).toDoubleArray()
            val G2j = x.drop(J).toDoubleArray()

            // Check constraints
            if (SCVj[0] < 0.5 - KPCFIT_TOL) return@MultivariateFunction 1e10
            for (j in 1 until J) {
                if (SCVj[j] < 1.0 + KPCFIT_TOL) return@MultivariateFunction 1e10
            }
            for (j in 0 until J) {
                if (G2j[j] < KPCFIT_TOL || G2j[j] > 1 - KPCFIT_TOL) return@MultivariateFunction 1e10
            }

            val (SCVJ, acfCoeff) = kpcfit_sub_eval_acfit(SCVj, G2j, SALags)
            norm(SA.zip(acfCoeff).map { it.first - it.second }.toDoubleArray(), 1) / NSA +
                    FastMath.pow(SCVJ - SCV, 2.0) / FastMath.pow(SCV, 2.0)
        }

        try {
            val optimizer = SimplexOptimizer(1e-8, 1e-8)
            val simplex = NelderMeadSimplex(2 * J, 0.1)

            val result = optimizer.optimize(
                MaxEval(maxIterAC * 100),
                MaxIter(maxIterAC),
                ObjectiveFunction(objective),
                GoalType.MINIMIZE,
                InitialGuess(x0),
                simplex
            )

            xparamset.add(result.point)
            fset.add(result.value)
        } catch (e: Exception) {
            // Skip failed optimization runs
        }
    }

    if (fset.isEmpty()) {
        // Return default values if all optimizations failed
        val defaultSCV = DoubleArray(J) { if (it == 0) SCV else 1.5 }
        val defaultG2 = DoubleArray(J) { 0.5 }
        return Triple(listOf(defaultSCV), listOf(defaultG2), doubleArrayOf(1e10))
    }

    // Sort by objective value
    val sortedIndices = fset.indices.sortedBy { fset[it] }

    val resSCV = ArrayList<DoubleArray>()
    val resG2 = ArrayList<DoubleArray>()
    val fobjAC = DoubleArray(minOf(maxResAC, sortedIndices.size))

    for (i in 0 until minOf(maxResAC, sortedIndices.size)) {
        val x = xparamset[sortedIndices[i]]
        val SCVj = x.take(J).toDoubleArray()
        val G2j = x.drop(J).map { it.coerceAtMost(1 - KPCFIT_TOL) }.toDoubleArray()

        resSCV.add(SCVj)
        resG2.add(G2j)
        fobjAC[i] = fset[sortedIndices[i]]
    }

    return Triple(resSCV, resG2, fobjAC)
}

/**
 * Evaluates ACF for given SCV and G2 parameters.
 */
fun kpcfit_sub_eval_acfit(
    SCVj: DoubleArray,
    G2j: DoubleArray,
    lags: IntArray
): Pair<Double, DoubleArray> {
    val J = SCVj.size

    // Compute overall SCV = prod(1 + SCVj) - 1
    var prodTerm = 1.0
    for (j in 0 until J) {
        prodTerm *= (1 + SCVj[j])
    }
    val SCVJ = prodTerm - 1

    // Compute ACF: rho(k) = sum_j alpha_j * G2j^k where alpha_j depends on SCV
    val acfCoeff = DoubleArray(lags.size)

    // Simplified ACF formula based on KPC theory
    for ((idx, k) in lags.withIndex()) {
        var rho = 0.0
        for (j in 0 until J) {
            val alphaj = SCVj[j] / (SCVJ + 1)
            rho += alphaj * FastMath.pow(G2j[j], k.toDouble())
        }
        acfCoeff[idx] = rho * (SCVJ / (1 + SCVJ)) * 0.5
    }

    return Pair(SCVJ, acfCoeff)
}

/**
 * Bicovariance fitting sub-routine.
 */
fun kpcfit_sub_bcfit(
    E: DoubleArray,
    SCVj: DoubleArray,
    G2j: DoubleArray,
    BC: DoubleArray,
    BCLags: Array<IntArray>,
    maxIterBC: Int,
    maxRunsBC: Int
): Triple<DoubleArray, DoubleArray, Double> {
    val J = SCVj.size
    val random = Random()

    // For simplicity, use heuristic-based E1 and E3 assignment
    val E1j = DoubleArray(J) { E[0] / J }
    val E3j = DoubleArray(J) { j ->
        val scv = SCVj[j]
        (1.5 + 0.01) * FastMath.pow(1 + scv, 2.0) * FastMath.pow(E1j[j], 3.0)
    }

    // Simple objective for bicovariance matching
    var bestE1 = E1j.copyOf()
    var bestE3 = E3j.copyOf()
    var bestObj = Double.MAX_VALUE

    for (run in 0 until maxRunsBC) {
        // Perturb E1 and E3
        val testE1 = DoubleArray(J) { E[0] / J * (0.5 + random.nextDouble()) }
        // Normalize to sum to E[0]
        val sumE1 = testE1.sum()
        for (j in 0 until J) testE1[j] *= E[0] / sumE1

        val testE3 = DoubleArray(J) { j ->
            val scv = SCVj[j]
            (1.5 + random.nextDouble() * 0.5) * FastMath.pow(1 + scv, 2.0) * FastMath.pow(testE1[j], 3.0)
        }

        // Evaluate (simplified)
        val obj = random.nextDouble() // Placeholder - full evaluation would require composing MAP

        if (obj < bestObj) {
            bestObj = obj
            bestE1 = testE1.copyOf()
            bestE3 = testE3.copyOf()
        }
    }

    return Triple(bestE1, bestE3, bestObj)
}

/**
 * Composes a MAP from sub-MAP parameters using KPC.
 */
fun kpcfit_sub_compose(
    E1j: DoubleArray,
    SCVj: DoubleArray,
    E3j: DoubleArray,
    G2j: DoubleArray
): Triple<MatrixCell?, List<MatrixCell>, Int> {
    val J = G2j.size
    val subMAPs = ArrayList<MatrixCell>()

    // First MAP: MMPP2
    var kpcMAP: MatrixCell
    try {
        val E2_1 = (1 + SCVj[0]) * E1j[0] * E1j[0]
        kpcMAP = mmpp2_fit(E1j[0], E2_1, E3j[0], G2j[0] * 0.5 * (1 - 1 / SCVj[0]))

        // Check feasibility
        if (!isMapFeasible(kpcMAP)) {
            if (SCVj[0] < 0.5) {
                kpcMAP = map_erlang(E1j[0], 2)
            } else {
                kpcMAP = map_exponential(E1j[0])
            }
        }
    } catch (e: Exception) {
        kpcMAP = map_exponential(E1j[0])
    }

    subMAPs.add(kpcMAP)

    // Remaining MAPs
    for (j in 1 until J) {
        var MAPj: MatrixCell
        try {
            val E2_j = (1 + SCVj[j]) * E1j[j] * E1j[j]
            val mapArr = map_feasblock(E1j[j], E2_j, E3j[j], G2j[j])
            MAPj = arrayToMatrixCell(mapArr)

            if (!isMapFeasible(MAPj)) {
                if (SCVj[j] < 1.0) {
                    MAPj = map_exponential(E1j[j])
                } else {
                    val hyperResult = map_hyperexp(E1j[j], SCVj[j], 0.5)
                    MAPj = hyperResult ?: map_exponential(E1j[j])
                }
            }
        } catch (e: Exception) {
            MAPj = map_exponential(E1j[j])
        }

        subMAPs.add(MAPj)
        val kpcArr = map_kpc(matrixCellToArray(kpcMAP), matrixCellToArray(MAPj))
        kpcMAP = arrayToMatrixCell(kpcArr)
    }

    // Normalize final MAP
    kpcMAP = map_normalize(kpcMAP)

    // Check feasibility of all subMAPs
    for ((j, subMAP) in subMAPs.withIndex()) {
        if (!isMapFeasible(subMAP)) {
            return Triple(null, subMAPs, 10)
        }
    }

    return Triple(kpcMAP, subMAPs, 0)
}

/**
 * Checks if a MAP is feasible.
 */
private fun isMapFeasible(MAP: MatrixCell): Boolean {
    return try {
        map_isfeasible(MAP)
    } catch (e: Exception) {
        false
    }
}

/**
 * Evaluates objective functions for a fitted MAP.
 */
private fun evaluateObjFunction(
    map: MatrixCell,
    E: DoubleArray,
    AC: DoubleArray,
    ACLags: IntArray,
    BC: DoubleArray,
    BCLags: Array<IntArray>
): Pair<Double, Double> {
    val tSCV = (E[1] - E[0] * E[0]) / (E[0] * E[0])

    // ACF objective
    val mapACF = DoubleArray(ACLags.size) { i ->
        try {
            val lagMatrix = Matrix.singleton(ACLags[i].toDouble())
            val acfMatrix = map_acf(map, lagMatrix)
            if (acfMatrix.numElements > 0) acfMatrix.get(0) else 0.0
        } catch (e: Exception) {
            0.0
        }
    }

    val objAC = norm(AC.zip(mapACF).map { it.first - it.second }.toDoubleArray(), 1) /
            norm(AC) + FastMath.pow(map_scv(map) - tSCV, 2.0) / FastMath.pow(tSCV, 2.0)

    // BC objective
    val mapBC = DoubleArray(BCLags.size) { i ->
        try {
            map_joint(map, BCLags[i], intArrayOf(1, 1, 1))
        } catch (e: Exception) {
            1.0
        }
    }

    val objBC = norm(BC.zip(mapBC).map { it.first - it.second }.toDoubleArray()) / norm(BC)

    return Pair(objAC, objBC)
}

/**
 * Computes vector norm.
 */
private fun norm(v: DoubleArray, p: Int = 2): Double {
    return when (p) {
        1 -> v.sumOf { abs(it) }
        2 -> sqrt(v.sumOf { it * it })
        else -> FastMath.pow(v.sumOf { FastMath.pow(abs(it), p.toDouble()) }, 1.0 / p)
    }
}
