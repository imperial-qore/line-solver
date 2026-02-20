package jline.lib.kpctoolbox.kpcfit

import jline.api.mam.*
import jline.lib.kpctoolbox.basic.logspacei
import jline.lib.kpctoolbox.mmpp.mmpp2_fit3
import jline.lib.kpctoolbox.trace.*
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import org.apache.commons.math3.analysis.MultivariateFunction
import org.apache.commons.math3.linear.MatrixUtils
import org.apache.commons.math3.linear.QRDecomposition
import org.apache.commons.math3.optim.*
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.NelderMeadSimplex
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.SimplexOptimizer
import org.apache.commons.math3.util.FastMath
import java.util.Random
import kotlin.math.abs
import kotlin.math.ceil
import kotlin.math.ln
import kotlin.math.sqrt

/**
 * KPC-Toolbox fitting functions.
 * Kronecker Product Composition (KPC) method for fitting MAPs to traces.
 *
 * Ported from MATLAB: matlab/lib/kpctoolbox/kpcfit/
 */

// Global tolerance for KPC fitting (Fix 4: changed from 1e-6 to 1e-10 to match MATLAB)
const val KPCFIT_TOL = 1e-10

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

    // Fix 5: use ceil() on n/nMinSupportAC to match MATLAB: ceil(n/nMinSupportAC)
    val defaultACLags = acLags ?: logspacei(1.0, ceil(n.toDouble() / nMinSupportAC).coerceAtLeast(1.0), 500)
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
    // Fix 5: use ceil() here too, matching MATLAB: 1:ceil(length(S)/nMinSupportAC)
    val ACFull = trace_acf(S, (1..ceil(n.toDouble() / nMinSupportAC).toInt().coerceAtLeast(1)).toList().toIntArray())

    // Fix 5: AC truncation -- MATLAB finds posmax = ACLags(find(abs(AC)<1e-6, 1))
    // then deletes lags STRICTLY GREATER than posmax, keeping the boundary lag.
    // Old JAR used take(cutIdx) which excluded the boundary. Fix: take(cutIdx + 1).
    var cutIdx = validACLags.size
    for (i in AC.indices) {
        if (abs(AC[i]) < 1e-6) {
            cutIdx = i + 1  // include the boundary lag itself
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
 *
 * Fix 2: Replaced heuristic with proper BIC regression algorithm matching MATLAB kpcfit_sub_bic.m.
 * For each candidate order, fits AR(p) model via OLS regression, computes residual SSE,
 * then applies Schwarz Bayesian Criterion: SBC = n*log(SSE) - n*log(n) + log(n)*order.
 * Returns log2 of the order that minimizes SBC.
 *
 * @param ACFull Full autocorrelation sequence (SA in MATLAB)
 * @param orders Candidate orders (powers of 2)
 * @return Nstar = log2(best order)
 */
fun kpcfit_sub_bic(ACFull: DoubleArray, orders: IntArray): Int {
    val nlags = ACFull.size
    val ordermax = orders.maxOrNull() ?: 2

    // Find first lag where ACF < 1e-6 (MATLAB: nlagsend = find(SA < 1e-6))
    var nlagsend = nlags
    for (i in ACFull.indices) {
        if (ACFull[i] < 1e-6) {
            // MATLAB: nlagsend = nlagsend(1)-1 + ordermax
            // nlagsend(1) is 1-based index; i is 0-based, so i+1 is 1-based
            nlagsend = i + ordermax
            break
        }
    }

    val NLAGSMAX = 10000

    // Build SAlags: 1-based lag indices for the regression
    var effectiveOrders = orders
    var SAlags: IntArray

    if (nlagsend > NLAGSMAX) {
        // Use log-spaced lags for very long sequences
        val logLags = logspacei(1.0, (nlagsend - ordermax).toDouble(), NLAGSMAX)
        SAlags = logLags.distinct().toIntArray()
    } else {
        if (nlagsend > ordermax) {
            SAlags = (1..(nlagsend - ordermax)).toList().toIntArray()
        } else {
            val effectiveOrdermax = nlagsend - 2
            if (effectiveOrdermax < 1) return 1
            SAlags = (1..(nlagsend - effectiveOrdermax).coerceAtLeast(1)).toList().toIntArray()
            effectiveOrders = orders.filter { it <= effectiveOrdermax }.toIntArray()
            if (effectiveOrders.isEmpty()) return 1
        }
    }

    val nSamples = SAlags.size
    if (nSamples < 2) return 1

    // Y = SA(SAlags_y) where SAlags_y = SAlags
    // ACFull is 0-based: ACFull[k] = ACF at lag k+1
    // MATLAB SA is 1-based: SA(k) = ACF at lag k
    // So SA(SAlags[i]) corresponds to ACFull[SAlags[i]-1]
    val Y = DoubleArray(nSamples)
    for (i in 0 until nSamples) {
        val idx = SAlags[i] - 1  // convert 1-based lag to 0-based array index
        if (idx >= 0 && idx < nlags) {
            Y[i] = ACFull[idx]
        }
    }

    // For each candidate order, compute SBC
    val SBC = DoubleArray(effectiveOrders.size) { Double.MAX_VALUE }

    for (j in effectiveOrders.indices) {
        val order = effectiveOrders[j]

        // Build X matrix: X(:, col) = SA(SAlags_y + col) for col = 1..order
        val X = Array(nSamples) { DoubleArray(order) }
        var validMatrix = true
        for (col in 0 until order) {
            for (row in 0 until nSamples) {
                // MATLAB: lags_x = SAlags_y + i (where i goes from 1 to order)
                val lagIdx = SAlags[row] + (col + 1) - 1  // (lag value + offset) converted to 0-based
                if (lagIdx >= 0 && lagIdx < nlags) {
                    X[row][col] = ACFull[lagIdx]
                } else {
                    validMatrix = false
                }
            }
        }

        if (!validMatrix) continue

        val resid = regressResiduals(Y, X)
        if (resid != null) {
            var sse = 0.0
            for (r in resid) sse += r * r
            if (sse > 0) {
                // SBC = n*log(SSE) - n*log(n) + log(n)*order
                SBC[j] = nSamples * ln(sse) - nSamples * ln(nSamples.toDouble()) + ln(nSamples.toDouble()) * order
            }
        }
    }

    // Find order with minimum SBC
    var bestIdx = 0
    var bestSBC = SBC[0]
    for (j in 1 until SBC.size) {
        if (SBC[j] < bestSBC) {
            bestSBC = SBC[j]
            bestIdx = j
        }
    }

    // Return Nstar = log2(bestOrder), matching MATLAB: Nstar = log2(orders(ind(1)))
    val bestOrder = effectiveOrders[bestIdx]
    return FastMath.round(FastMath.log(2.0, bestOrder.toDouble())).toInt()
}

/**
 * Performs OLS regression Y = X*b (no intercept) and returns residuals r = Y - X*b.
 * Matches MATLAB regress(Y, X).
 *
 * @param Y Response vector (nSamples)
 * @param X Predictor matrix (nSamples x order)
 * @return Residual vector, or null if regression fails
 */
private fun regressResiduals(Y: DoubleArray, X: Array<DoubleArray>): DoubleArray? {
    val n = Y.size
    if (n == 0 || X.isEmpty() || X[0].isEmpty()) return null
    val p = X[0].size

    return try {
        val xMatrix = MatrixUtils.createRealMatrix(n, p)
        for (i in 0 until n) {
            for (j in 0 until p) {
                xMatrix.setEntry(i, j, X[i][j])
            }
        }
        val yVector = MatrixUtils.createRealVector(Y)

        // Solve via QR decomposition: b = (X'X)^{-1} X'Y
        val qr = QRDecomposition(xMatrix)
        val b = qr.solver.solve(yVector)

        // Compute residuals: r = Y - X*b
        val fitted = xMatrix.operate(b)
        val residuals = DoubleArray(n)
        for (i in 0 until n) {
            residuals[i] = Y[i] - fitted.getEntry(i)
        }
        residuals
    } catch (e: Exception) {
        null
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
 *
 * Fix 1: Replaced simplified additive formula with exact iterative algorithm from MATLAB
 * kpcfit_sub_eval_iat2.m. The key difference is the cross-term SCVj_1*acfCoeff*(1+X)
 * in the numerator of the iterative update.
 *
 * MATLAB algorithm:
 *   SCVj = SCV(1);
 *   acfCoeff = 0.5*(1-1/SCV(1))*G2(1).^acfLag;
 *   for j=2:J
 *       SCVj_1 = SCVj;
 *       SCVj = (1+SCVj)*(1+SCV(j))/2 - 1;
 *       r0j = 0.5*(1-1/SCV(j));
 *       X = SCV(j)*r0j*G2(j).^acfLag;
 *       acfCoeff = (X + SCVj_1*acfCoeff.*(1+X)) / SCVj;
 *   end
 */
fun kpcfit_sub_eval_acfit(
    SCVj: DoubleArray,
    G2j: DoubleArray,
    lags: IntArray
): Pair<Double, DoubleArray> {
    val J = SCVj.size

    // Initialize with first MAP
    var SCVcum = SCVj[0]

    // acfCoeff = 0.5*(1-1/SCV(1))*G2(1).^acfLag
    val acfCoeff = DoubleArray(lags.size) { idx ->
        0.5 * (1.0 - 1.0 / SCVj[0]) * FastMath.pow(G2j[0], lags[idx].toDouble())
    }

    // Iterative composition for j=2..J (0-based: j=1..J-1)
    for (j in 1 until J) {
        val SCVj_1 = SCVcum
        SCVcum = (1.0 + SCVcum) * (1.0 + SCVj[j]) / 2.0 - 1.0
        val r0j = 0.5 * (1.0 - 1.0 / SCVj[j])

        for (idx in lags.indices) {
            val X = SCVj[j] * r0j * FastMath.pow(G2j[j], lags[idx].toDouble())
            acfCoeff[idx] = (X + SCVj_1 * acfCoeff[idx] * (1.0 + X)) / SCVcum
        }
    }

    return Pair(SCVcum, acfCoeff)
}

/**
 * Bicovariance fitting sub-routine.
 *
 * Fix 3: Replaced random stub with proper optimization matching MATLAB kpcfit_sub_bcfit.m.
 * Uses Nelder-Mead optimization (approximating MATLAB fmincon) with:
 * - Proper initialization: E1j = E(1)^(1/J) * ones(1,J)
 * - Actual bicovariance objective via kpcfit_sub_compose + map_joint
 * - Multiple random restarts with best-of tracking
 * - Constraint enforcement via penalty terms
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
    val NumMAPs = SCVj.size
    val TOL = 1e-9
    val EPSTOL = 10 * TOL

    // Normalizing constant for bicovariances
    val NBC = norm(BC)
    if (NBC < 1e-30) {
        // Degenerate case
        val E1j = DoubleArray(NumMAPs) { FastMath.pow(E[0], 1.0 / NumMAPs) }
        val E3j = DoubleArray(NumMAPs) { j ->
            val E2j = (1 + SCVj[j]) * E1j[j] * E1j[j]
            1.501 * E2j * E2j / E1j[j]
        }
        return Triple(E1j, E3j, 0.0)
    }

    // Truncate bicorrelations when lags have negative diffs (matching MATLAB)
    val validIndices = ArrayList<Int>()
    for (index in BCLags.indices) {
        val lags = BCLags[index]
        var hasNegDiff = false
        for (i in 1 until lags.size) {
            if (lags[i] - lags[i - 1] < 0) {
                hasNegDiff = true
                break
            }
        }
        if (!hasNegDiff) {
            validIndices.add(index)
        }
    }
    val filteredBC = DoubleArray(validIndices.size) { BC[validIndices[it]] }
    val filteredBCLags = Array(validIndices.size) { BCLags[validIndices[it]] }
    val filteredNBC = norm(filteredBC)
    if (filteredNBC < 1e-30) {
        val E1j = DoubleArray(NumMAPs) { FastMath.pow(E[0], 1.0 / NumMAPs) }
        val E3j = DoubleArray(NumMAPs) { j ->
            val E2j = (1 + SCVj[j]) * E1j[j] * E1j[j]
            1.501 * E2j * E2j / E1j[j]
        }
        return Triple(E1j, E3j, 0.0)
    }

    // Initialize: E1j = E(1)^(1/NumMAPs) * ones(1,NumMAPs)
    val E1jBase = FastMath.pow(E[0], 1.0 / NumMAPs)
    val E1j0 = DoubleArray(NumMAPs) { E1jBase }

    val random = Random()
    var t = E[0] * random.nextDouble()

    // E2j for each MAP (used for E3j initialization)
    val E2j0 = DoubleArray(NumMAPs) { j -> (1 + SCVj[j]) * E1j0[j] * E1j0[j] }
    val E3j0 = DoubleArray(NumMAPs) { j -> (1.5 + t) * E2j0[j] * E2j0[j] / E1j0[j] }

    // x = [E1j(1), ..., E1j(NumMAPs), E3j(1), ..., E3j(NumMAPs)]
    val x0base = DoubleArray(2 * NumMAPs)
    for (j in 0 until NumMAPs) {
        x0base[j] = E1j0[j]
        x0base[NumMAPs + j] = E3j0[j]
    }

    var fBest = Double.MAX_VALUE
    var xBest: DoubleArray? = null
    var fold = 0.0

    // Helper to extract E1j, E3j from optimization vector x
    fun xtopar(x: DoubleArray): Pair<DoubleArray, DoubleArray> {
        val e1 = DoubleArray(NumMAPs) { x[it] }
        val e3 = DoubleArray(NumMAPs) { x[NumMAPs + it] }
        // E1j(1) = E(1) / prod(E1j(2:end))
        var prodE1 = 1.0
        for (j in 1 until NumMAPs) {
            prodE1 *= e1[j]
        }
        e1[0] = if (prodE1 > 0) E[0] / prodE1 else E1jBase
        return Pair(e1, e3)
    }

    // Objective function matching MATLAB objfun
    fun objfun(x: DoubleArray): Double {
        val (e1, e3) = xtopar(x)

        // Check positivity
        for (j in 0 until NumMAPs) {
            if (e1[j] <= EPSTOL || e3[j] <= EPSTOL) return maxOf(2 * fold, 1e6)
        }

        // Check constraints (penalty for violations, matching MATLAB nnlcon)
        val e2 = DoubleArray(NumMAPs) { j -> (1 + SCVj[j]) * e1[j] * e1[j] }
        for (j in 1 until NumMAPs) {
            if ((2 + EPSTOL) * e1[j] * e1[j] > e2[j]) return maxOf(2 * fold, 1e6)
            if ((1.5 + EPSTOL) * e2[j] * e2[j] / e1[j] > e3[j]) return maxOf(2 * fold, 1e6)
        }
        if (SCVj[0] > 1) {
            if ((2 + EPSTOL) * e1[0] * e1[0] > e2[0]) return maxOf(2 * fold, 1e6)
            if ((1.5 + EPSTOL) * e2[0] * e2[0] / e1[0] > e3[0]) return maxOf(2 * fold, 1e6)
        }

        // Compose MAP
        val (compMAP, _, err) = kpcfit_sub_compose(e1, SCVj, e3, G2j)
        if (err != 0 || compMAP == null) {
            return maxOf(2 * fold, 1e6)
        }

        // Scale MAP to match E(1)
        val scaledMAP = map_scale(compMAP, E[0])

        // Compute bicovariance values
        val BCj = DoubleArray(filteredBCLags.size)
        for (indexL in filteredBCLags.indices) {
            BCj[indexL] = try {
                map_joint(scaledMAP, filteredBCLags[indexL], intArrayOf(1, 1, 1))
            } catch (ex: Exception) {
                return maxOf(2 * fold, 1e6)
            }
        }

        val f = norm(DoubleArray(filteredBC.size) { filteredBC[it] - BCj[it] }) / filteredNBC
        return if (f.isNaN()) {
            2 * fold
        } else {
            fold = f
            f
        }
    }

    for (ind in 0 until maxRunsBC) {
        // Build x0 for this run (matching MATLAB restart logic)
        val x0: DoubleArray
        if (ind == 0) {
            x0 = x0base.copyOf()
        } else {
            x0 = DoubleArray(2 * NumMAPs)
            // MATLAB: x0 = x0base .* [(0.25+1.75*rand(NumMAPs,1)); ones(NumMAPs,1)]'
            for (j in 0 until NumMAPs) {
                x0[j] = x0base[j] * (0.25 + 1.75 * random.nextDouble())
            }
            // MATLAB: t = rand*E(1); E2j = (1+SCVj)*x0^2; x0(NumMAPs+j) = (3/2+t)*E2j^2/x0(j)
            t = random.nextDouble() * E[0]
            for (j in 0 until NumMAPs) {
                val e2jLocal = (1 + SCVj[j]) * x0[j] * x0[j]
                x0[NumMAPs + j] = (1.5 + t) * e2jLocal * e2jLocal / x0[j]
            }
        }

        // Ensure all values are positive
        for (j in 0 until 2 * NumMAPs) {
            if (x0[j] <= EPSTOL) x0[j] = EPSTOL
        }

        try {
            val objective = MultivariateFunction { pt -> objfun(pt) }
            val optimizer = SimplexOptimizer(TOL, TOL)
            val simplex = NelderMeadSimplex(2 * NumMAPs, 0.01)

            val result = optimizer.optimize(
                MaxEval(maxIterBC * 1000),
                MaxIter(maxIterBC),
                ObjectiveFunction(objective),
                GoalType.MINIMIZE,
                InitialGuess(x0),
                simplex
            )

            if (result.value < fBest) {
                fBest = result.value
                xBest = result.point
            }
        } catch (e: Exception) {
            // Skip failed optimization runs
        }
    }

    // Extract best parameters
    if (xBest != null) {
        val (e1Best, e3Best) = xtopar(xBest!!)
        // MATLAB: E1j(1)=E(1)/prod(E1j(2:end)) -- already done in xtopar
        return Triple(e1Best, e3Best, fBest)
    }

    // Fallback: return initial values
    val (e1Fall, e3Fall) = xtopar(x0base)
    return Triple(e1Fall, e3Fall, Double.MAX_VALUE)
}

/**
 * Composes a MAP from sub-MAP parameters using KPC.
 *
 * Fix 6: Added map2_fit fallback cascade for infeasible first MAP with SCV >= 0.5,
 * matching MATLAB kpcfit_sub_compose.m. Also added full fallback cascade for j>=2 MAPs
 * including eigenvalue normalization step.
 */
fun kpcfit_sub_compose(
    E1j: DoubleArray,
    SCVj: DoubleArray,
    E3j: DoubleArray,
    G2j: DoubleArray
): Triple<MatrixCell?, List<MatrixCell>, Int> {
    val J = G2j.size
    val subMAPs = ArrayList<MatrixCell>()

    // First MAP: MMPP2 using mmpp2_fit3 (matches MATLAB exactly)
    var kpcMAP: MatrixCell
    try {
        val E2_1 = (1 + SCVj[0]) * E1j[0] * E1j[0]
        kpcMAP = mmpp2_fit3(E1j[0], E2_1, E3j[0], G2j[0])

        // Check feasibility (MATLAB: any(imag(D0)>1e-4) || any(imag(D1)>1e-4) || !feasible)
        val hasImagOrNaN = hasNaNOrInfinite(kpcMAP)
        if (hasImagOrNaN || !isMapFeasible(kpcMAP)) {
            if (SCVj[0] < 0.5) {
                kpcMAP = map_erlang(E1j[0], 2)
            } else {
                // Fix 6: map2_fit fallback cascade matching MATLAB
                // Try map2_fit(E1, E2, -1, G2) first
                val fitResult1 = map2_fit(E1j[0], E2_1, -1.0, G2j[0])
                if (fitResult1.MAP != null && fitResult1.MAP.size() > 0 && fitResult1.error.toInt() == 0) {
                    kpcMAP = fitResult1.MAP
                } else {
                    // Then try map2_fit(E1, E2, -1, 0)
                    val fitResult2 = map2_fit(E1j[0], E2_1, -1.0, 0.0)
                    if (fitResult2.MAP != null && fitResult2.MAP.size() > 0 && fitResult2.error.toInt() == 0) {
                        kpcMAP = fitResult2.MAP
                    } else {
                        return Triple(null, subMAPs, 1)
                    }
                }
            }
        }
    } catch (e: Exception) {
        // Fallback: try map2_fit cascade before giving up
        try {
            val E2_1 = (1 + SCVj[0]) * E1j[0] * E1j[0]
            if (SCVj[0] < 0.5) {
                kpcMAP = map_erlang(E1j[0], 2)
            } else {
                val fitResult1 = map2_fit(E1j[0], E2_1, -1.0, G2j[0])
                if (fitResult1.MAP != null && fitResult1.MAP.size() > 0 && fitResult1.error.toInt() == 0) {
                    kpcMAP = fitResult1.MAP
                } else {
                    val fitResult2 = map2_fit(E1j[0], E2_1, -1.0, 0.0)
                    if (fitResult2.MAP != null && fitResult2.MAP.size() > 0 && fitResult2.error.toInt() == 0) {
                        kpcMAP = fitResult2.MAP
                    } else {
                        return Triple(null, ArrayList(), 1)
                    }
                }
            }
        } catch (e2: Exception) {
            return Triple(null, ArrayList(), 1)
        }
    }

    subMAPs.add(kpcMAP)

    // Remaining MAPs (j=2..J in MATLAB, j=1..J-1 in 0-based)
    for (j in 1 until J) {
        var MAPj: MatrixCell?
        try {
            val E2_j = (1 + SCVj[j]) * E1j[j] * E1j[j]
            val mapArr = map_feasblock(E1j[j], E2_j, E3j[j], G2j[j])
            MAPj = arrayToMatrixCell(mapArr)

            val feasible = try {
                isMapFeasible(MAPj)
            } catch (ex: Exception) {
                false
            }

            if (!feasible) {
                if (SCVj[j] < 1.0) {
                    // Low variability: use exponential (MATLAB: map2_exponential)
                    MAPj = map_exponential(E1j[j])
                } else {
                    // Fix 6: Full fallback cascade matching MATLAB for j>=2 MAPs
                    // Try map2_fit(E1, E2, -1, G2)
                    val fitResult1 = map2_fit(E1j[j], E2_j, -1.0, G2j[j])
                    if (fitResult1.error.toInt() != 0 || fitResult1.MAP == null || fitResult1.MAP.size() == 0) {
                        // Try map2_fit(E1, E2, -1, 0)
                        val fitResult2 = map2_fit(E1j[j], E2_j, -1.0, 0.0)
                        if (fitResult2.error.toInt() != 0 || fitResult2.MAP == null || fitResult2.MAP.size() == 0) {
                            return Triple(null, subMAPs, 5)
                        }
                        MAPj = fitResult2.MAP
                    } else {
                        MAPj = fitResult1.MAP
                    }
                    // MATLAB eigenvalue normalization:
                    // lambda = eig(-inv(MAPj{1})); D0 = diag(-1./lambda);
                    // P = map_embedded(MAPj); p = min(eig(P));
                    // D1 = -D0*[p,1-p;p,1-p]; MAPj = map_normalize({D0,D1})
                    try {
                        val D0orig = MAPj!!.get(0)
                        val n = D0orig.numRows
                        if (n == 2) {
                            // Compute -inv(D0)
                            val negInvD0 = D0orig.inv()
                            for (ii in 0 until n) {
                                for (jj in 0 until n) {
                                    negInvD0.set(ii, jj, -negInvD0.get(ii, jj))
                                }
                            }
                            val eigenValues = negInvD0.eig()
                            if (eigenValues != null && eigenValues.size >= n) {
                                val D0new = Matrix(n, n)
                                for (ii in 0 until n) {
                                    D0new.set(ii, ii, -1.0 / eigenValues[ii].real)
                                }
                                val P = map_embedded(MAPj)
                                val pEig = P.eig()
                                if (pEig != null && pEig.size >= n) {
                                    var pMinReal = pEig[0].real
                                    for (ii in 1 until pEig.size) {
                                        if (pEig[ii].real < pMinReal) pMinReal = pEig[ii].real
                                    }
                                    val D1new = Matrix(n, n)
                                    for (ii in 0 until n) {
                                        D1new.set(ii, 0, -D0new.get(ii, ii) * pMinReal)
                                        D1new.set(ii, 1, -D0new.get(ii, ii) * (1 - pMinReal))
                                    }
                                    val newCell = MatrixCell(2)
                                    newCell[0] = D0new
                                    newCell[1] = D1new
                                    MAPj = map_normalize(newCell)
                                }
                            }
                        }
                    } catch (ex: Exception) {
                        // If eigenvalue normalization fails, keep the map2_fit result
                    }
                }
            }
        } catch (e: Exception) {
            MAPj = map_exponential(E1j[j])
        }

        // MATLAB: if isempty(MAPj), subMAPs{j}=map_exponential(E1j(j))
        if (MAPj == null || MAPj.size() == 0) {
            MAPj = map_exponential(E1j[j])
        }
        subMAPs.add(MAPj)

        val kpcArr = map_kpc(matrixCellToArray(kpcMAP), matrixCellToArray(MAPj))
        kpcMAP = arrayToMatrixCell(kpcArr)
    }

    // Normalize final MAP
    kpcMAP = map_normalize(kpcMAP)

    // Check feasibility of all subMAPs
    var error = 0
    for (subMAP in subMAPs) {
        if (!isMapFeasible(subMAP)) {
            error = 10
        }
    }

    return Triple(kpcMAP, subMAPs, error)
}

/**
 * Checks if a MAP has NaN or Infinite entries (corresponding to MATLAB imaginary check).
 */
private fun hasNaNOrInfinite(MAP: MatrixCell): Boolean {
    try {
        val D0 = MAP.get(0)
        val D1 = MAP.get(1)
        for (i in 0 until D0.numRows) {
            for (j in 0 until D0.numCols) {
                if (D0.get(i, j).isNaN() || D0.get(i, j).isInfinite()) return true
            }
        }
        for (i in 0 until D1.numRows) {
            for (j in 0 until D1.numCols) {
                if (D1.get(i, j).isNaN() || D1.get(i, j).isInfinite()) return true
            }
        }
    } catch (e: Exception) {
        return true
    }
    return false
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

// =============================================================================
// PH Fitting Functions
// =============================================================================

/**
 * Computes the characteristic polynomial coefficients for a hyper-exponential
 * distribution from its moments, using the Prony method.
 *
 * Ported from MATLAB: kpcfit_hyper_charpoly.m
 *
 * @param E Array of moments E[X], E[X^2], ..., E[X^(2n-1)]
 * @param n Number of exponential phases
 * @return Coefficients of the characteristic polynomial (n+1 elements)
 */
fun kpcfit_hyper_charpoly(E: DoubleArray, n: Int): DoubleArray {
    // MATLAB: E=[1,E(:)'] -- prepend 1 to make it 1-indexed
    // Ep[0]=1, Ep[1]=E[0], Ep[2]=E[1], ...
    val Ep = DoubleArray(1 + E.size)
    Ep[0] = 1.0
    for (i in E.indices) {
        Ep[i + 1] = E[i]
    }

    // f = factorial(0:(2*n-1))
    val f = DoubleArray(2 * n)
    for (i in 0 until 2 * n) {
        f[i] = factorial(i)
    }

    // Build A matrix (n+1 rows, n+1 cols)
    // MATLAB: A(i,:) = E((n+i):-1:i)./f((n+i):-1:i) for i=1..n
    val A = Array(n + 1) { DoubleArray(n + 1) }
    for (i in 1..n) {
        // MATLAB indices: (n+i):-1:i => n+i, n+i-1, ..., i (total n+1 elements)
        for (col in 0..n) {
            val epIdx = (n + i) - col  // MATLAB 1-based: (n+i)-col => 0-based: (n+i)-col-1+1 = (n+i)-col
            val fIdx = (n + i) - col   // same indexing
            // MATLAB E is 1-based; Ep is 0-based version of prepended E
            // MATLAB f is 1-based: f(k) = factorial(k-1) for k=1..2n
            // but we defined f as 0-based: f[k] = factorial(k)
            // MATLAB: E((n+i)-col) / f((n+i)-col)
            // E is 1-based in MATLAB with prepended 1, so E(k) = Ep[k-1]
            // f is 1-based: f(k) = factorial(k-1), so f[k-1] in our array
            val epIndex = epIdx - 1  // convert to 0-based
            val fIndex = fIdx - 1    // convert to 0-based
            if (epIndex >= 0 && epIndex < Ep.size && fIndex >= 0 && fIndex < f.size) {
                A[i - 1][col] = Ep[epIndex] / f[fIndex]
            }
        }
    }
    // MATLAB: A(n+1,:) = zeros(1,n+1); A(n+1,1) = 1
    for (col in 0..n) {
        A[n][col] = 0.0
    }
    A[n][0] = 1.0

    // b = zeros(n+1,1); b(n+1) = 1
    val b = DoubleArray(n + 1)
    b[n] = 1.0

    // m = A\b  (solve linear system)
    val aMatrix = MatrixUtils.createRealMatrix(n + 1, n + 1)
    for (i in 0..n) {
        for (j in 0..n) {
            aMatrix.setEntry(i, j, A[i][j])
        }
    }
    val bVector = MatrixUtils.createRealVector(b)
    val qr = QRDecomposition(aMatrix)
    val mVector = qr.getSolver().solve(bVector)

    val m = DoubleArray(n + 1)
    for (i in 0..n) {
        m[i] = mVector.getEntry(i)
    }
    return m
}

/**
 * Fits a hyper-exponential PH distribution using Prony's method.
 *
 * Ported from MATLAB: kpcfit_ph_prony.m
 *
 * @param E Array of moments E[X], E[X^2], ..., E[X^(2n-1)]
 * @param n Number of exponential phases
 * @return MatrixCell with D0 and D1 matrices of the PH distribution
 */
fun kpcfit_ph_prony(E: DoubleArray, n: Int): MatrixCell {
    // f = factorial(0:(2*n-1))
    val f = DoubleArray(2 * n)
    for (i in 0 until 2 * n) {
        f[i] = factorial(i)
    }

    // m = kpcfit_hyper_charpoly(E,n)
    val m = kpcfit_hyper_charpoly(E, n)

    // theta = roots(m)
    // MATLAB roots() takes descending coefficients: m(1)*x^n + m(2)*x^(n-1) + ... + m(n+1)
    // Maths.roots() takes ascending coefficients: a[0] + a[1]*x + ... + a[n]*x^n
    // So we need to reverse m for Maths.roots
    val mReversed = DoubleArray(m.size)
    for (i in m.indices) {
        mReversed[i] = m[m.size - 1 - i]
    }
    val thetaComplex = jline.util.Maths.roots(mReversed)

    // Take real parts of theta (for hyper-exp, roots should be real and positive)
    val theta = DoubleArray(n)
    for (i in 0 until n) {
        theta[i] = thetaComplex[i].real
    }

    // C(i,1:n) = f(i+1)*theta.^i for i=1..n
    // MATLAB f is 1-based: f(i+1) = factorial(i), our f is 0-based: f[i] = factorial(i)
    val C = Array(n) { DoubleArray(n) }
    for (i in 1..n) {
        for (j in 0 until n) {
            C[i - 1][j] = f[i] * FastMath.pow(theta[j], i.toDouble())
        }
    }

    // M = (C\E(1:n)')' -- solve C*M' = E(1:n)', then transpose
    // E(1:n) in MATLAB is 1-based, so E(1)..E(n) = our E[0]..E[n-1]
    val cMatrix = MatrixUtils.createRealMatrix(n, n)
    for (i in 0 until n) {
        for (j in 0 until n) {
            cMatrix.setEntry(i, j, C[i][j])
        }
    }
    val eVector = MatrixUtils.createRealVector(DoubleArray(n) { E[it] })
    val qr = QRDecomposition(cMatrix)
    val mSolve = qr.getSolver().solve(eVector)

    // M is a 1xn row vector (entry probabilities)
    val M = DoubleArray(n) { mSolve.getEntry(it) }

    // PH{1} = diag(-1./theta)
    val D0 = Matrix(n, n)
    for (i in 0 until n) {
        D0.set(i, i, -1.0 / theta[i])
    }

    // PH{2} = -diag(-1./theta)*ones(n,1)*M
    // = diag(1./theta)*ones(n,1)*M
    // ones(n,1)*M is an nxn matrix where each row is M
    // diag(1./theta) * (ones(n,1)*M) = matrix where row i = (1/theta[i]) * M
    val D1 = Matrix(n, n)
    for (i in 0 until n) {
        for (j in 0 until n) {
            D1.set(i, j, (1.0 / theta[i]) * M[j])
        }
    }

    val PH = MatrixCell(2)
    PH[0] = D0
    PH[1] = D1
    return PH
}

/**
 * Options for PH distribution fitting.
 */
data class KPCFitPhOptions(
    var verbose: Boolean = true,
    var runs: Int = 5,
    var minNumStates: Int = 2,
    var maxNumStates: Int = 32,
    var minExactMom: Int = 3
)

/**
 * Factory function for KPCFitPhOptions with validation.
 * Validates options (power-of-2 checks, range checks, moment sufficiency checks).
 *
 * Ported from MATLAB: kpcfit_ph_options.m
 *
 * @param E Array of moments
 * @param verbose Whether to print verbose output
 * @param runs Maximum number of optimization runs
 * @param minNumStates Minimum number of states (must be power of 2)
 * @param maxNumStates Maximum number of states (must be power of 2)
 * @param minExactMom Minimum number of moments to be fitted exactly
 * @return Validated KPCFitPhOptions
 */
fun kpcfit_ph_options(
    E: DoubleArray,
    verbose: Boolean = true,
    runs: Int = 5,
    minNumStates: Int = 2,
    maxNumStates: Int = 32,
    minExactMom: Int = 3
): KPCFitPhOptions {
    val options = KPCFitPhOptions(verbose, runs, minNumStates, maxNumStates, minExactMom)

    // Range checks
    if (options.runs < 1) options.runs = 1
    if (options.minNumStates < 2) options.minNumStates = 2
    if (options.maxNumStates < 2) options.maxNumStates = 2
    if (options.minExactMom < 1) options.minExactMom = 1
    if (options.minExactMom > E.size) options.minExactMom = E.size

    // Power-of-2 checks for MinNumStates
    if (options.minNumStates and (options.minNumStates - 1) != 0) {
        // Not a power of 2
        options.minNumStates = nextPowerOf2(options.minNumStates)
        if (options.verbose) {
            System.err.println("Warning: MinNumStates not a power of 2, fixed to ${options.minNumStates}.")
        }
    }

    // Power-of-2 checks for MaxNumStates
    if (options.maxNumStates and (options.maxNumStates - 1) != 0) {
        options.maxNumStates = nextPowerOf2(options.maxNumStates)
        if (options.verbose) {
            System.err.println("Warning: MaxNumStates not a power of 2, fixed to ${options.maxNumStates}.")
        }
    }

    // Check moment sufficiency
    if (2 * options.maxNumStates - 1 > E.size) {
        throw IllegalArgumentException(
            "MaxNumStates of ${options.maxNumStates} requires at least ${2 * options.maxNumStates - 1} moments, but only ${E.size} provided."
        )
    }

    // Ensure MinNumStates <= MaxNumStates
    if (options.minNumStates > options.maxNumStates) {
        if (options.verbose) {
            System.err.println("Warning: MaxNumStates < MinNumStates, fixed.")
        }
        val tmp = options.minNumStates
        options.minNumStates = options.maxNumStates
        options.maxNumStates = tmp
    }

    return options
}

/**
 * Returns the next power of 2 >= n.
 */
private fun nextPowerOf2(n: Int): Int {
    var v = n - 1
    v = v or (v shr 1)
    v = v or (v shr 2)
    v = v or (v shr 4)
    v = v or (v shr 8)
    v = v or (v shr 16)
    return v + 1
}

/**
 * Computes factorial as a Double.
 */
private fun factorial(n: Int): Double {
    if (n < 0) return 1.0
    if (n <= 1) return 1.0
    var result = 1.0
    for (i in 2..n) {
        result *= i.toDouble()
    }
    return result
}

/**
 * Exact PH fitting using Prony's method (SCV > 1), PH(2) fitting (SCV < 1),
 * or exponential (SCV == 1).
 *
 * Ported from MATLAB: kpcfit_ph_exact.m
 *
 * @param E Array of moments (normalized to mean=1)
 * @param options Fitting options
 * @return List of feasible PH distributions as MatrixCells
 */
fun kpcfit_ph_exact(E: DoubleArray, options: KPCFitPhOptions): List<MatrixCell> {
    val phExact = ArrayList<MatrixCell>()
    val SCV = (E[1] - E[0] * E[0]) / (E[0] * E[0])

    if (SCV > 1.0) {
        // Higher variability than exponential: Prony's method for hyper-exponentials
        if (options.verbose) {
            println("kpcfit_ph: HIGHER variability than an exponential (var/mean^2 = $SCV)")
            println()
            println("kpcfit_ph: starting exact hyper-exponential fitting method (Prony's method)")
        }
        for (n in 2..options.maxNumStates) {
            if (E.size < 2 * n - 1) {
                if (options.verbose) {
                    println("kpcfit_ph: not enough moments given in input to fit hyper-exp($n)")
                }
                break
            }
            try {
                val PH = kpcfit_ph_prony(E, n)
                if (map_isfeasible(PH)) {
                    phExact.add(PH)
                    if (options.verbose) {
                        println("\t\t\thyper-exp($n): feasible, matched exactly ${2 * n - 1} moments. result saved.")
                    }
                } else {
                    if (options.verbose) {
                        println("\t\t\thyper-exp($n): infeasible to fit exactly.")
                    }
                    break
                }
            } catch (e: Exception) {
                if (options.verbose) {
                    println("\t\t\thyper-exp($n): infeasible to fit exactly.")
                }
                break
            }
        }
    } else if (SCV < 1.0) {
        // Lower variability than exponential
        if (options.verbose) {
            println("kpcfit_ph: LOWER variability than an exponential (var/mean^2 = $SCV)")
            println()
        }
        // Find minimum n such that 1/n <= SCV
        var n = 1
        while (1.0 / n > SCV) {
            n++
        }
        if (options.verbose) {
            println("kpcfit_ph: exact fitting of E[X^2] requires at least $n states")
        }
        if (options.minExactMom >= 2 && options.maxNumStates < n) {
            if (options.verbose) {
                println("kpcfit_ph: impossible to fit exactly E[X^2] with MaxNumStates = ${options.maxNumStates}, increasing to MaxNumStates = $n.")
            }
            return phExact
        }

        if (n == 2) {
            if (options.verbose) {
                println("kpcfit_ph: attempting PH(2) fitting method")
            }
            // Try map2_fit with 3 moments
            val fitResult = map2_fit(E[0], E[1], E[2], 0.0)
            if (fitResult.MAP == null || fitResult.MAP.size() == 0) {
                // Fall back to 2 moments
                val fitResult2 = map2_fit(E[0], E[1], -1.0, 0.0)
                if (fitResult2.MAP != null && fitResult2.MAP.size() > 0 && map_isfeasible(fitResult2.MAP)) {
                    phExact.add(fitResult2.MAP)
                    if (options.verbose) {
                        println("\t\t\tph(2): feasible, matched exactly 2 moments. result saved.")
                    }
                } else {
                    throw IllegalStateException("anomalous set of moments, please check.")
                }
            } else {
                phExact.add(fitResult.MAP)
                if (options.verbose) {
                    println("\t\t\tph(2): feasible, matched exactly 3 moments. result saved.")
                }
            }
        } else if (abs(SCV - 1.0 / n) < KPCFIT_TOL) {
            // Erlang distribution
            val ERL = map_erlang(E[0], n)
            // Check if moments match
            var momDist = 0.0
            for (k in 1..E.size) {
                val diff = E[k - 1] - map_moment(ERL, k)
                momDist += diff * diff
            }
            if (sqrt(momDist) < KPCFIT_TOL) {
                if (options.verbose) {
                    println("kpcfit_ph: erlang moment set. fitted erlang-$n. result saved.")
                }
                phExact.add(ERL)
            }
        } else {
            // APH fitting
            val maxorder = options.maxNumStates
            if (options.verbose) {
                println("kpcfit_ph: fitting APH distribution (best effort, max order = $maxorder).")
            }
            try {
                val PH = aph_fit(E[0], E[1], E[2], maxorder)
                if (map_isfeasible(PH)) {
                    var aphMatched = 0
                    val phSize = PH[0].numRows
                    for (k in 1..(2 * phSize - 1)) {
                        if (k <= E.size) {
                            val momK = map_moment(PH, k)
                            if (abs(E[k - 1] - momK) < KPCFIT_TOL * momK) {
                                aphMatched++
                            }
                        }
                    }
                    if (options.verbose) {
                        println("\t\t\t      aph($phSize): feasible, matched exactly $aphMatched moments. result saved.")
                    }
                    phExact.add(PH)
                } else {
                    if (options.verbose) {
                        println("kpcfit_ph: cannot fit APH distribution.")
                    }
                }
            } catch (e: Exception) {
                if (options.verbose) {
                    println("kpcfit_ph: cannot fit APH distribution.")
                }
            }
        }
    } else {
        // SCV == 1: exponential
        if (options.verbose) {
            println("kpcfit_ph: SAME variability as an exponential (var/mean^2 = $SCV)")
            println()
        }
        val EXP = MatrixCell(2)
        EXP[0] = Matrix(1, 1)
        EXP[0].set(0, 0, -1.0 / E[0])
        EXP[1] = Matrix(1, 1)
        EXP[1].set(0, 0, 1.0 / E[0])
        // Check if moments match
        var momDist = 0.0
        for (k in 1..E.size) {
            val diff = E[k - 1] - map_moment(EXP, k)
            momDist += diff * diff
        }
        if (sqrt(momDist) < KPCFIT_TOL) {
            if (options.verbose) {
                println("kpcfit_ph: exponential moment set. fitted exponential. result saved.")
            }
        }
        phExact.add(EXP)
    }

    return phExact
}

/**
 * Optimization-based search for PH distributions using KPC composition.
 * Works in log-moment space with penalty-based constraint enforcement.
 *
 * Ported from MATLAB: kpcfit_ph_search.m
 *
 * @param E Array of moments (normalized to mean=1)
 * @param J Number of PH(2) sub-distributions to compose (total states = 2^J)
 * @param options Fitting options
 * @param x0 Initial point (optional)
 * @param maxAphOrder Maximum APH order for SCV < 1 case (optional)
 * @return Triple of (composed MAP, objective value, optimization variable x)
 */
fun kpcfit_ph_search(
    E: DoubleArray,
    J: Int,
    options: KPCFitPhOptions,
    x0: DoubleArray? = null,
    maxAphOrder: Int? = null
): Triple<MatrixCell, Double, DoubleArray> {
    val SCV = (E[1] - E[0] * E[0]) / (E[0] * E[0])
    val actualMaxAphOrder = maxAphOrder ?: ceil(1.0 / SCV).toInt().coerceAtLeast(2)

    val order = IntArray(J) { 2 }  // Current version assumes all PH(2)s
    val K = maxOf(2 * order.max()!! - 1, E.size)
    val F = DoubleArray(K)
    for (i in 0 until K) {
        F[i] = factorial(i + 1)
    }
    val logF = DoubleArray(K) { ln(F[it]) }

    // Weight vector: w = (log(E(end)).^(1/length(E))).^-(1:length(E))
    val logElast = ln(E[E.size - 1])
    val base = FastMath.pow(logElast, 1.0 / E.size)
    val w = DoubleArray(E.size) { k -> FastMath.pow(base, -(k + 1).toDouble()) }
    val logE = DoubleArray(E.size) { ln(E[it]) }

    // Initialize sub-PHs
    val PHs = ArrayList<MatrixCell>()
    for (j in 0 until J) {
        val fbResult = map_feasblock(Math.random(), 1000 * Math.random(), -1.0, 0.0)
        PHs.add(arrayToMatrixCell(fbResult))
    }

    // Build initial logEtable (J x K)
    val logEtable = Array(J) { DoubleArray(K) }
    for (k in 0 until K) {
        for (j in 0 until J) {
            val mom = map_moment(PHs[j], k + 1)
            logEtable[j][k] = if (mom > 0) ln(mom) else -200.0
        }
    }

    // Use initial point if provided
    val x0Flat: DoubleArray
    if (x0 != null && x0.isNotEmpty()) {
        x0Flat = x0.copyOf()
        // Populate logEtable from x0
        for (j in 0 until J) {
            for (k in 0 until K) {
                val idx = j * K + k
                if (idx < x0Flat.size) {
                    logEtable[j][k] = x0Flat[idx]
                }
            }
        }
    } else {
        // Flatten logEtable as initial point
        x0Flat = DoubleArray(J * K)
        for (j in 0 until J) {
            for (k in 0 until K) {
                x0Flat[j * K + k] = logEtable[j][k]
            }
        }
    }

    // Stagnation detection
    var stagnVal = 0.0
    var stagnIter = 0

    // Objective function (with penalty for constraints)
    val objective = MultivariateFunction { x ->
        // Check for NaN
        if (x.any { it.isNaN() }) return@MultivariateFunction 1e10

        // Reshape x to logEtable
        val logEt = Array(J) { j -> DoubleArray(K) { k -> x[j * K + k] } }

        // Reconstruct PHs and update logEtable based on SCV
        if (SCV < 1.0) {
            // First PH: use aph_fit
            try {
                val phj = map_scale(
                    aph_fit(FastMath.exp(logEt[0][0]), FastMath.exp(logEt[0][1]), FastMath.exp(logEt[0][2]), actualMaxAphOrder),
                    FastMath.exp(logEt[0][0])
                )
                for (k in 0 until K) {
                    val mom = map_moment(phj, k + 1)
                    logEt[0][k] = if (mom > 0) ln(mom) else -200.0
                }
            } catch (e: Exception) {
                return@MultivariateFunction 1e10
            }
            // Remaining PHs: use map_feasblock
            for (j in 1 until J) {
                try {
                    val phj = map_scale(
                        arrayToMatrixCell(map_feasblock(FastMath.exp(logEt[j][0]), FastMath.exp(logEt[j][1]), FastMath.exp(logEt[j][2]), 0.0)),
                        FastMath.exp(logEt[j][0])
                    )
                    for (k in 0 until K) {
                        val mom = map_moment(phj, k + 1)
                        logEt[j][k] = if (mom > 0) ln(mom) else -200.0
                    }
                } catch (e: Exception) {
                    return@MultivariateFunction 1e10
                }
            }
        } else {
            // SCV >= 1: use characteristic polynomial method for j>=2
            for (j in 1 until J) {
                try {
                    // Compute characteristic polynomial coefficients
                    val expLogEt = DoubleArray(K) { k -> FastMath.exp(logEt[j][k]) / F[k] }
                    val m = kpcfit_hyper_charpoly(expLogEt, 2)

                    if (m.any { it.isNaN() }) {
                        // Reset to exponential moments
                        val expPH = MatrixCell(2)
                        expPH[0] = Matrix(1, 1); expPH[0].set(0, 0, -1.0)
                        expPH[1] = Matrix(1, 1); expPH[1].set(0, 0, 1.0)
                        for (k in 3 until K) {
                            val idxStart = k - 1
                            val idxEnd = k - (2 * order[j] - 2)
                            if (idxEnd >= 0 && idxStart < K) {
                                var sum = 0.0
                                for (idx in idxStart downTo maxOf(idxEnd, 0)) {
                                    val Ej = map_moment(expPH, idx + 1)
                                    val mIdx = idxStart - idx + 1
                                    if (mIdx < m.size) {
                                        sum += m[mIdx] * (Ej / F[idx])
                                    }
                                }
                                val valK = -F[k] * sum
                                logEt[j][k] = if (valK > 0) ln(valK) else -200.0
                            }
                        }
                    } else {
                        for (k in 3 until K) {
                            val idxStart = k - 1
                            val idxEnd = k - (2 * order[j] - 2)
                            if (idxEnd >= 0 && idxStart < K) {
                                var sum = 0.0
                                for (idx in idxStart downTo maxOf(idxEnd, 0)) {
                                    val Ej = FastMath.exp(logEt[j][idx])
                                    if (Ej.isNaN()) {
                                        // Reset to exponential
                                        val expMom = map_moment(map_exponential(1.0), idx + 1)
                                        val mIdx = idxStart - idx + 1
                                        if (mIdx < m.size) {
                                            sum += m[mIdx] * (expMom / F[idx])
                                        }
                                    } else {
                                        val mIdx = idxStart - idx + 1
                                        if (mIdx < m.size) {
                                            sum += m[mIdx] * (Ej / F[idx])
                                        }
                                    }
                                }
                                val valK = -F[k] * sum
                                logEt[j][k] = if (valK > 0) ln(valK) else -200.0
                            }
                        }
                    }
                } catch (e: Exception) {
                    // Keep existing logEt values
                }
            }
        }

        // Compute composed log moments: logEcur = sum(logEtable,1) - (J-1)*log(F)
        val logEcur = DoubleArray(K) { k ->
            var sum = 0.0
            for (j in 0 until J) {
                sum += logEt[j][k]
            }
            sum - (J - 1) * logF[k]
        }

        // Objective: f = w * abs( log(E) - logEcur )'
        var f = 0.0
        for (k in 0 until E.size) {
            f += w[k] * abs(logE[k] - logEcur[k])
        }

        // Penalty for constraints
        var penalty = 0.0

        // Equality constraints: first MinExactMom moments matched
        for (k in 0 until minOf(options.minExactMom, E.size)) {
            val violation = w[k] * abs(logE[k] - logEcur[k])
            penalty += 1000.0 * violation * violation
        }

        // Inequality constraints
        if (SCV < 1.0) {
            // APH constraints for first PH
            val le2 = logEt[0][1]
            val le1 = logEt[0][0]
            val le3 = logEt[0][2]
            val n0 = actualMaxAphOrder
            // log n2 >= log n - log(n-1)
            val c1 = -(le2 - ln(n0.toDouble()) + ln((n0 - 1).toDouble()))
            if (c1 > 0) penalty += 1000.0 * c1 * c1
            // log e3 - 2*log e2 + log e1 >= log(n+1) - log(n)
            val c2 = -(le3 - 2.0 * le2 + le1 - ln((n0 + 1).toDouble()) + ln(n0.toDouble()))
            if (c2 > 0) penalty += 1000.0 * c2 * c2
        } else {
            // Hyper-exp constraints for first PH
            val c1 = -(logEt[0][1] - 2.0 * logEt[0][0] - ln(2.0)) - 10 * KPCFIT_TOL
            if (c1 > 0) penalty += 1000.0 * c1 * c1
            val c2 = -(logEt[0][2] - (ln(1.5) + 2.0 * logEt[0][1] - logEt[0][0]))
            if (c2 > 0) penalty += 1000.0 * c2 * c2
        }

        // Constraints for j>=2 (always hyper-exp)
        for (j in 1 until J) {
            val c1 = -(logEt[j][1] - 2.0 * logEt[j][0] - ln(2.0)) - 10 * KPCFIT_TOL
            if (c1 > 0) penalty += 1000.0 * c1 * c1
            val c2 = -(logEt[j][2] - (ln(1.5) + 2.0 * logEt[j][1] - logEt[j][0]))
            if (c2 > 0) penalty += 1000.0 * c2 * c2
        }

        val total = f + penalty
        if (total.isNaN()) 1e10 else total
    }

    // Bounds: -200 to 200 for all variables
    val lowerBounds = DoubleArray(J * K) { -200.0 }
    val upperBounds = DoubleArray(J * K) { 200.0 }

    // Clamp initial point to bounds
    val x0Clamped = DoubleArray(J * K) { i ->
        x0Flat.getOrElse(i) { 0.0 }.coerceIn(lowerBounds[i], upperBounds[i])
    }

    // Optimize using Nelder-Mead
    var bestX = x0Clamped.copyOf()
    var bestScore = 1e10

    try {
        val optimizer = SimplexOptimizer(KPCFIT_TOL, KPCFIT_TOL)
        val simplex = NelderMeadSimplex(J * K, 0.5)

        val result = optimizer.optimize(
            MaxEval(200 * 100),
            MaxIter(200),
            ObjectiveFunction(objective),
            GoalType.MINIMIZE,
            InitialGuess(x0Clamped),
            simplex
        )

        bestX = result.point
        bestScore = result.value
    } catch (e: Exception) {
        // Use initial point if optimization fails
    }

    // Reconstruct PHs from best solution
    val logEtBest = Array(J) { j -> DoubleArray(K) { k -> bestX[j * K + k] } }
    val subPHs = ArrayList<MatrixCell>()

    if (SCV < 1.0) {
        try {
            val ph1 = map_scale(
                aph_fit(FastMath.exp(logEtBest[0][0]), FastMath.exp(logEtBest[0][1]), FastMath.exp(logEtBest[0][2]), actualMaxAphOrder),
                FastMath.exp(logEtBest[0][0])
            )
            subPHs.add(ph1)
        } catch (e: Exception) {
            subPHs.add(map_exponential(FastMath.exp(logEtBest[0][0])))
        }
    } else {
        try {
            val fb = map_feasblock(FastMath.exp(logEtBest[0][0]), FastMath.exp(logEtBest[0][1]), FastMath.exp(logEtBest[0][2]), 0.0)
            subPHs.add(arrayToMatrixCell(fb))
        } catch (e: Exception) {
            subPHs.add(map_exponential(FastMath.exp(logEtBest[0][0])))
        }
    }

    for (j in 1 until J) {
        try {
            val fb = map_feasblock(FastMath.exp(logEtBest[j][0]), FastMath.exp(logEtBest[j][1]), FastMath.exp(logEtBest[j][2]), 0.0)
            subPHs.add(arrayToMatrixCell(fb))
        } catch (e: Exception) {
            subPHs.add(map_exponential(FastMath.exp(logEtBest[j][0])))
        }
    }

    // Compose all sub-PHs via KPC
    var kpcPH = subPHs[0]
    for (j in 1 until J) {
        val kpcArr = map_kpc(matrixCellToArray(kpcPH), matrixCellToArray(subPHs[j]))
        kpcPH = arrayToMatrixCell(kpcArr)
    }

    return Triple(kpcPH, bestScore, bestX)
}

/**
 * Fits a hyperexponential distribution with n states using constrained optimization
 * in parameter space with weighted log-moment objective.
 *
 * Ported from MATLAB: psfit_hyperexp_weighted.m
 *
 * @param E Array of moments
 * @param n Number of states (exponential phases)
 * @return Triple of (PH distribution, approximate moments, optimization variables)
 */
fun psfit_hyperexp_weighted(E: DoubleArray, n: Int): Triple<MatrixCell, DoubleArray, DoubleArray> {
    // Scale moments to mean=1
    val Eunscaled = E.copyOf()
    val Escale = DoubleArray(E.size) { k -> FastMath.pow(E[0], (k + 1).toDouble()) }
    val Enormalized = DoubleArray(E.size) { k -> E[k] / Escale[k] }

    // Weight vector
    val logElast = ln(Enormalized[Enormalized.size - 1])
    val base = FastMath.pow(logElast, 1.0 / Enormalized.size)
    val w = DoubleArray(Enormalized.size) { k -> FastMath.pow(base, -(k + 1).toDouble()) }

    val TOL = 1e-8
    val EPSTOL = 100 * TOL
    val MAXITER = 100
    val random = Random()

    // Helper: convert x to alpha and T
    // alpha = x[0:n], T = diag(-1/x[n:2n])
    fun topar(x: DoubleArray): Pair<DoubleArray, Matrix> {
        val alpha = DoubleArray(n) { x[it] }
        val T = Matrix(n, n)
        for (i in 0 until n) {
            T.set(i, i, -1.0 / x[n + i])
        }
        return Pair(alpha, T)
    }

    // Compute approximate moments from alpha and T
    fun computeMoments(alpha: DoubleArray, T: Matrix, numMoments: Int): DoubleArray {
        val Eapx = DoubleArray(numMoments)
        val negInvT = T.inv()
        for (i in 0 until n) {
            for (j in 0 until n) {
                negInvT.set(i, j, -negInvT.get(i, j))
            }
        }
        val onesVec = Matrix(n, 1)
        for (i in 0 until n) onesVec.set(i, 0, 1.0)
        val alphaMatrix = Matrix(1, n)
        for (i in 0 until n) alphaMatrix.set(0, i, alpha[i])

        var negInvTpow = Matrix.eye(n)
        for (j in 1..numMoments) {
            negInvTpow = negInvTpow.mult(negInvT)
            // E[X^j] = factorial(j) * alpha * (-inv(T))^j * ones(n,1)
            val temp = alphaMatrix.mult(negInvTpow).mult(onesVec)
            Eapx[j - 1] = factorial(j) * temp.get(0, 0)
        }
        return Eapx
    }

    // Objective function with penalty
    val objective = MultivariateFunction { x ->
        // Check bounds
        for (i in 0 until 2 * n) {
            if (x[i] < EPSTOL) return@MultivariateFunction 1e10
        }

        val (alpha, T) = topar(x)

        // Check alpha >= 0
        for (i in 0 until n) {
            if (alpha[i] < 0) return@MultivariateFunction 1e10
        }

        val Eapx = try {
            computeMoments(alpha, T, Enormalized.size)
        } catch (e: Exception) {
            return@MultivariateFunction 1e10
        }

        // Check for non-positive moments
        for (mom in Eapx) {
            if (mom <= 0 || mom.isNaN() || mom.isInfinite()) return@MultivariateFunction 1e10
        }

        // Objective: w * abs(log(E) - log(Eapx))
        var f = 0.0
        for (k in Enormalized.indices) {
            f += w[k] * abs(ln(Enormalized[k]) - ln(Eapx[k]))
        }

        // Penalty for constraints
        var penalty = 0.0

        // Equality constraints: first 3 moments matched
        for (k in 0 until minOf(3, Enormalized.size)) {
            val ceq = w[k] * abs(ln(Enormalized[k]) - ln(Eapx[k]))
            penalty += 10000.0 * ceq * ceq
        }

        // Equality constraint: sum(alpha) = 1
        val alphaSum = alpha.sum()
        penalty += 10000.0 * (alphaSum - 1.0) * (alphaSum - 1.0)

        // Inequality constraints: alpha >= 0 (already checked above)

        if (f.isNaN()) 1e10 else f + penalty
    }

    var bestX: DoubleArray? = null
    var bestF = Double.MAX_VALUE

    // Multiple random restarts
    for (run in 0 until 5) {
        val x0 = DoubleArray(2 * n)
        for (i in 0 until n) {
            x0[i] = random.nextDouble()
        }
        // Normalize alpha to sum to 1
        val alphaSum = x0.take(n).sum()
        for (i in 0 until n) {
            x0[i] /= alphaSum
        }
        for (i in n until 2 * n) {
            x0[i] = random.nextDouble()
        }
        // Enforce lower bound
        for (i in 0 until 2 * n) {
            if (x0[i] < EPSTOL) x0[i] = EPSTOL
        }

        try {
            val optimizer = SimplexOptimizer(TOL, TOL)
            val simplex = NelderMeadSimplex(2 * n, 0.01)

            val result = optimizer.optimize(
                MaxEval(MAXITER * 200),
                MaxIter(MAXITER),
                ObjectiveFunction(objective),
                GoalType.MINIMIZE,
                InitialGuess(x0),
                simplex
            )

            if (result.value < bestF) {
                bestF = result.value
                bestX = result.point
            }
        } catch (e: Exception) {
            // Skip failed runs
        }
    }

    if (bestX == null) {
        // Return exponential as fallback
        val expPH = map_exponential(Eunscaled[0])
        val Eapx = DoubleArray(E.size) { k -> map_moment(expPH, k + 1) }
        return Triple(expPH, Eapx, DoubleArray(2 * n))
    }

    val (alpha, T) = topar(bestX!!)
    val Eapx = computeMoments(alpha, T, Enormalized.size)

    // Build MAP = {T, -T*ones(n,1)*alpha}
    val onesVec = Matrix(n, 1)
    for (i in 0 until n) onesVec.set(i, 0, 1.0)
    val alphaMatrix = Matrix(1, n)
    for (i in 0 until n) alphaMatrix.set(0, i, alpha[i])
    val negTones = T.mult(onesVec)
    for (i in 0 until n) negTones.set(i, 0, -negTones.get(i, 0))
    val D1 = negTones.mult(alphaMatrix)
    val MAP = MatrixCell(2)
    MAP[0] = T
    MAP[1] = D1

    // Scale to original mean and normalize
    val scaledMAP = map_scale(MAP, Eunscaled[0])
    val normalizedMAP = map_normalize(scaledMAP)

    // Compute final approximate moments (unscaled)
    val EapxUnscaled = DoubleArray(E.size) { k -> Eapx[k] * Escale[k] }

    return Triple(normalizedMAP, EapxUnscaled, bestX!!)
}

/**
 * Fits an acyclic PH (APH) distribution with n states using constrained optimization
 * in parameter space with weighted log-moment objective.
 *
 * The difference from psfit_hyperexp_weighted is the T matrix structure:
 * T = diag(-1/x[n:2n]) + diag(1/x[n:2n-1], 1)  (upper bidiagonal)
 *
 * Ported from MATLAB: psfit_aph_weighted.m
 *
 * @param E Array of moments
 * @param n Number of states
 * @return Triple of (PH distribution, approximate moments, optimization variables)
 */
fun psfit_aph_weighted(E: DoubleArray, n: Int): Triple<MatrixCell, DoubleArray, DoubleArray> {
    // Scale moments to mean=1
    val Eunscaled = E.copyOf()
    val Escale = DoubleArray(E.size) { k -> FastMath.pow(E[0], (k + 1).toDouble()) }
    val Enormalized = DoubleArray(E.size) { k -> E[k] / Escale[k] }

    // Weight vector
    val logElast = ln(Enormalized[Enormalized.size - 1])
    val base = FastMath.pow(logElast, 1.0 / Enormalized.size)
    val w = DoubleArray(Enormalized.size) { k -> FastMath.pow(base, -(k + 1).toDouble()) }

    val TOL = 1e-8
    val EPSTOL = 100 * TOL
    val MAXITER = 100
    val random = Random()

    // Helper: convert x to alpha and T (APH structure: upper bidiagonal)
    // alpha = x[0:n], T = diag(-1/x[n:2n]) + diag(1/x[n:2n-1], 1)
    fun topar(x: DoubleArray): Pair<DoubleArray, Matrix> {
        val alpha = DoubleArray(n) { x[it] }
        val T = Matrix(n, n)
        for (i in 0 until n) {
            T.set(i, i, -1.0 / x[n + i])
        }
        // Upper diagonal: diag(1/x[n:2n-1], 1)
        for (i in 0 until n - 1) {
            T.set(i, i + 1, 1.0 / x[n + i])
        }
        return Pair(alpha, T)
    }

    // Compute approximate moments from alpha and T
    fun computeMoments(alpha: DoubleArray, T: Matrix, numMoments: Int): DoubleArray {
        val Eapx = DoubleArray(numMoments)
        val negInvT = T.inv()
        for (i in 0 until n) {
            for (j in 0 until n) {
                negInvT.set(i, j, -negInvT.get(i, j))
            }
        }
        val onesVec = Matrix(n, 1)
        for (i in 0 until n) onesVec.set(i, 0, 1.0)
        val alphaMatrix = Matrix(1, n)
        for (i in 0 until n) alphaMatrix.set(0, i, alpha[i])

        var negInvTpow = Matrix.eye(n)
        for (j in 1..numMoments) {
            negInvTpow = negInvTpow.mult(negInvT)
            val temp = alphaMatrix.mult(negInvTpow).mult(onesVec)
            Eapx[j - 1] = factorial(j) * temp.get(0, 0)
        }
        return Eapx
    }

    // Objective function with penalty
    val objective = MultivariateFunction { x ->
        for (i in 0 until 2 * n) {
            if (x[i] < EPSTOL) return@MultivariateFunction 1e10
        }

        val (alpha, T) = topar(x)
        for (i in 0 until n) {
            if (alpha[i] < 0) return@MultivariateFunction 1e10
        }

        val Eapx = try {
            computeMoments(alpha, T, Enormalized.size)
        } catch (e: Exception) {
            return@MultivariateFunction 1e10
        }

        for (mom in Eapx) {
            if (mom <= 0 || mom.isNaN() || mom.isInfinite()) return@MultivariateFunction 1e10
        }

        var f = 0.0
        for (k in Enormalized.indices) {
            f += w[k] * abs(ln(Enormalized[k]) - ln(Eapx[k]))
        }

        var penalty = 0.0
        // Equality constraints: first 3 moments
        for (k in 0 until minOf(3, Enormalized.size)) {
            val ceq = w[k] * abs(ln(Enormalized[k]) - ln(Eapx[k]))
            penalty += 10000.0 * ceq * ceq
        }
        // sum(alpha) = 1
        val alphaSum = alpha.sum()
        penalty += 10000.0 * (alphaSum - 1.0) * (alphaSum - 1.0)

        if (f.isNaN()) 1e10 else f + penalty
    }

    var bestX: DoubleArray? = null
    var bestF = Double.MAX_VALUE

    for (run in 0 until 5) {
        val x0 = DoubleArray(2 * n)
        for (i in 0 until n) {
            x0[i] = random.nextDouble()
        }
        val alphaSum = x0.take(n).sum()
        for (i in 0 until n) {
            x0[i] /= alphaSum
        }
        for (i in n until 2 * n) {
            x0[i] = random.nextDouble()
        }
        for (i in 0 until 2 * n) {
            if (x0[i] < EPSTOL) x0[i] = EPSTOL
        }

        try {
            val optimizer = SimplexOptimizer(TOL, TOL)
            val simplex = NelderMeadSimplex(2 * n, 0.01)

            val result = optimizer.optimize(
                MaxEval(MAXITER * 200),
                MaxIter(MAXITER),
                ObjectiveFunction(objective),
                GoalType.MINIMIZE,
                InitialGuess(x0),
                simplex
            )

            if (result.value < bestF) {
                bestF = result.value
                bestX = result.point
            }
        } catch (e: Exception) {
            // Skip failed runs
        }
    }

    if (bestX == null) {
        val expPH = map_exponential(Eunscaled[0])
        val Eapx = DoubleArray(E.size) { k -> map_moment(expPH, k + 1) }
        return Triple(expPH, Eapx, DoubleArray(2 * n))
    }

    val (alpha, T) = topar(bestX!!)
    val Eapx = computeMoments(alpha, T, Enormalized.size)

    // Build MAP = {T, -T*ones(n,1)*alpha}
    val onesVec = Matrix(n, 1)
    for (i in 0 until n) onesVec.set(i, 0, 1.0)
    val alphaMatrix = Matrix(1, n)
    for (i in 0 until n) alphaMatrix.set(0, i, alpha[i])
    val negTones = T.mult(onesVec)
    for (i in 0 until n) negTones.set(i, 0, -negTones.get(i, 0))
    val D1 = negTones.mult(alphaMatrix)
    val MAP = MatrixCell(2)
    MAP[0] = T
    MAP[1] = D1

    val scaledMAP = map_scale(MAP, Eunscaled[0])
    val normalizedMAP = map_normalize(scaledMAP)

    val EapxUnscaled = DoubleArray(E.size) { k -> Eapx[k] * Escale[k] }

    return Triple(normalizedMAP, EapxUnscaled, bestX!!)
}

/**
 * Automatic PH distribution fitting from moments.
 * Combines exact fitting, KPC-based approximate fitting (moment space),
 * and parameter-space fitting (hyperexp or APH).
 *
 * Ported from MATLAB: kpcfit_ph_auto.m
 *
 * @param E Array of moments E[X], E[X^2], ..., E[X^K]
 * @param options Fitting options
 * @return Array of rows where each row is: [PH (MatrixCell), dist (Double), x (DoubleArray?), method (String)]
 */
fun kpcfit_ph_auto(E: DoubleArray, options: KPCFitPhOptions): Array<Array<Any?>> {
    if (options.verbose) {
        println("kpcfit_ph: PH fitting")
        println("kpcfit_ph: type \"help kpcfit_ph_options\" for options")
    }

    // Scale moments such that E[X]=1
    val Eunscaled = E.copyOf()
    val Escale = DoubleArray(E.size) { k -> FastMath.pow(E[0], (k + 1).toDouble()) }
    val Enormalized = DoubleArray(E.size) { k -> E[k] / Escale[k] }

    // Weight function for distance computation
    val logElast10 = FastMath.log10(Enormalized[Enormalized.size - 1])
    val base10 = FastMath.pow(logElast10, 1.0 / Enormalized.size)
    fun distFun(Eref: DoubleArray, Eapx: DoubleArray): Double {
        val wLocal = DoubleArray(Eref.size) { k -> FastMath.pow(base10, -(k + 1).toDouble()) }
        var dist = 0.0
        for (k in Eref.indices) {
            val logRef = if (Eref[k] > 0) FastMath.log10(Eref[k]) else -300.0
            val logApx = if (Eapx[k] > 0) FastMath.log10(Eapx[k]) else -300.0
            dist += wLocal[k] * abs(logRef - logApx)
        }
        return dist
    }

    // ---- EXACT FITTING ----
    if (options.verbose) {
        println()
        println("kpcfit_ph: starting exact fitting methods")
    }
    val phExact = kpcfit_ph_exact(Enormalized, options)

    // ---- APPROXIMATE FITTING -- MOMENT SPACE (KPC-based) ----
    if (options.verbose) {
        println()
        println("kpcfit_ph: starting approximate fitting method -- moment space (KPC-based)")
        println("\t\t\tRUN\t\tORD\tDIST\t\tRESULT")
    }
    val phApxMS = ArrayList<Array<Any?>>()  // Each: [MAP, dist, x]

    val minJ = maxOf(2, ceil(FastMath.log(2.0, options.minNumStates.toDouble())).toInt())
    val maxJ = ceil(FastMath.log(2.0, options.maxNumStates.toDouble())).toInt()

    for (J in minJ..maxJ) {
        var mindist = Double.MAX_VALUE
        var best = 0
        var bestIdx = 0
        val randomruns = ceil(2.0 * options.runs / 3.0).toInt()

        for (run in 1..options.runs) {
            if (options.verbose) {
                if (run > randomruns) {
                    print("\n\t\t\t${run}r\t")
                } else {
                    print("\n\t\t\t$run\t")
                }
                print("\t${1 shl J}")
            }

            val x0Input: DoubleArray?
            if (run > randomruns && best > 0 && bestIdx > 0 && bestIdx <= phApxMS.size) {
                x0Input = phApxMS[bestIdx - 1][2] as? DoubleArray
            } else {
                x0Input = null
            }

            try {
                val (APX, _, x) = kpcfit_ph_search(Enormalized, J, options, x0Input)
                if (map_isfeasible(APX)) {
                    // Compute distance
                    val apxMoments = DoubleArray(Enormalized.size) { k -> map_moment(APX, k + 1) }
                    val dist = distFun(Enormalized, apxMoments)

                    if (options.verbose) {
                        print("\t${String.format("%6.6f", dist)}")
                    }

                    if (dist < mindist) {
                        if (best == 0) {
                            phApxMS.add(arrayOf<Any?>(APX, dist, x))
                            if (options.verbose) {
                                val apxSize = APX[0].numRows
                                print("\to PH($apxSize) saved.")
                            }
                        } else {
                            phApxMS[phApxMS.size - 1] = arrayOf<Any?>(APX, dist, x)
                            if (options.verbose) {
                                val apxSize = APX[0].numRows
                                print("\t+ PH($apxSize) updated.")
                            }
                        }
                        best = run
                        bestIdx = phApxMS.size
                        mindist = dist
                    }
                } else {
                    if (options.verbose) {
                        print("\t\t\t\tx infeasible result.")
                    }
                }
            } catch (e: Exception) {
                if (options.verbose) {
                    print("\t\t\t\tx optimization solver failed.")
                }
            }
        }
    }

    // ---- APPROXIMATE FITTING -- PARAMETER SPACE ----
    val SCV = (Enormalized[1] - Enormalized[0] * Enormalized[0]) / (Enormalized[0] * Enormalized[0])
    val phApxPS = ArrayList<Array<Any?>>()

    if (SCV >= 1.0) {
        if (options.verbose) {
            println()
            println()
            println("kpcfit_ph: starting approximate fitting method -- hyper-exponential parameter space (non-KPC-based)")
            println("\t\t\tRUN\t\tORD\tDIST\t\tRESULT")
        }

        for (J in minJ..maxJ) {
            var mindist = Double.MAX_VALUE
            var best = 0
            var bestIdx = 0

            for (run in 1..options.runs) {
                if (options.verbose) {
                    print("\n\t\t\t$run\t")
                    print("\t${1 shl J}")
                }
                try {
                    val (APX, _, x) = psfit_hyperexp_weighted(Enormalized, 1 shl J)
                    if (map_isfeasible(APX)) {
                        val apxMoments = DoubleArray(Enormalized.size) { k -> map_moment(APX, k + 1) }
                        val dist = distFun(Enormalized, apxMoments)
                        if (options.verbose) {
                            print("\t${String.format("%6.6f", dist)}")
                        }
                        if (dist < mindist) {
                            if (best == 0) {
                                phApxPS.add(arrayOf<Any?>(APX, dist, x))
                                if (options.verbose) {
                                    print("\to PH(${APX[0].numRows}) saved.")
                                }
                            } else {
                                phApxPS[phApxPS.size - 1] = arrayOf<Any?>(APX, dist, x)
                                if (options.verbose) {
                                    print("\t+ PH(${APX[0].numRows}) updated.")
                                }
                            }
                            best = run
                            bestIdx = phApxPS.size
                            mindist = dist
                        }
                    } else {
                        if (options.verbose) {
                            print("\t\t\t\tinfeasible result.")
                        }
                    }
                } catch (e: Exception) {
                    if (options.verbose) {
                        print("\t\t\t\tx optimization solver failed.")
                    }
                }
            }
        }
    } else {
        if (options.verbose) {
            println()
            println()
            println("kpcfit_ph: starting approximate fitting method -- aph parameter space (non-KPC-based)")
            println("\t\t\tRUN\t\tORD\tDIST\t\tRESULT")
        }

        for (J in minJ..maxJ) {
            var mindist = Double.MAX_VALUE
            var best = 0
            var bestIdx = 0

            for (run in 1..options.runs) {
                if (options.verbose) {
                    print("\n\t\t\t$run\t")
                    print("\t${1 shl J}")
                }
                try {
                    val (APX, _, x) = psfit_aph_weighted(Enormalized, 1 shl J)
                    if (map_isfeasible(APX)) {
                        val apxMoments = DoubleArray(Enormalized.size) { k -> map_moment(APX, k + 1) }
                        val dist = distFun(Enormalized, apxMoments)
                        if (options.verbose) {
                            print("\t${String.format("%6.6f", dist)}")
                        }
                        if (dist < mindist) {
                            if (best == 0) {
                                phApxPS.add(arrayOf<Any?>(APX, dist, x))
                                if (options.verbose) {
                                    print("\to PH(${APX[0].numRows}) saved.")
                                }
                            } else {
                                phApxPS[phApxPS.size - 1] = arrayOf<Any?>(APX, dist, x)
                                if (options.verbose) {
                                    print("\t+ PH(${APX[0].numRows}) updated.")
                                }
                            }
                            best = run
                            bestIdx = phApxPS.size
                            mindist = dist
                        }
                    } else {
                        if (options.verbose) {
                            print("\t\t\t\tinfeasible result.")
                        }
                    }
                } catch (e: Exception) {
                    if (options.verbose) {
                        print("\t\t\t\tx optimization solver failed.")
                    }
                }
            }
        }
    }

    // ---- ASSEMBLE RESULTS ----
    // Ensure all returned results match at least the mean (scale to original E[X])
    val nExactResults = phExact.size
    val nApproxMS = phApxMS.size

    val results = ArrayList<Array<Any?>>()

    // Exact results
    for (j in 0 until nExactResults) {
        val scaledPH = map_scale(phExact[j], Eunscaled[0])
        val apxMoments = DoubleArray(Enormalized.size) { k -> map_moment(phExact[j], k + 1) }
        val dist = distFun(Enormalized, apxMoments)
        results.add(arrayOf<Any?>(scaledPH, dist, null, "exact"))
    }

    // Approximate moment space results
    for (j in 0 until nApproxMS) {
        val apxMAP = phApxMS[j][0] as MatrixCell
        val scaledPH = map_scale(apxMAP, Eunscaled[0])
        val apxMoments = DoubleArray(Enormalized.size) { k -> map_moment(apxMAP, k + 1) }
        val dist = distFun(Enormalized, apxMoments)
        results.add(arrayOf<Any?>(scaledPH, dist, phApxMS[j][2], "approx_moment_space"))
    }

    // Approximate parameter space results
    for (j in 0 until phApxPS.size) {
        val apxMAP = phApxPS[j][0] as MatrixCell
        val scaledPH = map_scale(apxMAP, Eunscaled[0])
        val apxMoments = DoubleArray(Enormalized.size) { k -> map_moment(apxMAP, k + 1) }
        val dist = distFun(Enormalized, apxMoments)
        results.add(arrayOf<Any?>(scaledPH, dist, phApxPS[j][2], "approx_param_space"))
    }

    if (options.verbose) {
        println()
        println()
        println("Returned $nExactResults PH distribution(s) by exact moment matching.")
        println("Returned $nApproxMS PH distribution(s) by approximate moment matching (moment space).")
        println("Returned ${phApxPS.size} PH distribution(s) by approximate moment matching (parameter space).")
        println()
    }

    return results.toTypedArray()
}

/**
 * Returns a formatted summary string for PH fitting results.
 * Analogous to the MATLAB kpcfit_ph_summary script (without plotting).
 *
 * @param PH Results array from kpcfit_ph_auto: each row is [PH, dist, x, method]
 * @param E Original moments array
 * @return Formatted summary string
 */
fun kpcfit_ph_summary(PH: Array<Array<Any?>>, E: DoubleArray): String {
    val sb = StringBuilder()

    // Categorize results
    val exactEntries = ArrayList<Int>()
    val apxMSEntries = ArrayList<Int>()
    val apxPSEntries = ArrayList<Int>()
    var bestExactDist = Double.MAX_VALUE
    var bestExactPos = -1
    var bestApxMSDist = Double.MAX_VALUE
    var bestApxMSPos = -1
    var bestApxPSDist = Double.MAX_VALUE
    var bestApxPSPos = -1

    for (i in PH.indices) {
        val method = PH[i][3] as String
        val dist = PH[i][1] as Double
        val phMAP = PH[i][0] as MatrixCell

        when {
            method.equals("exact", ignoreCase = true) -> {
                exactEntries.add(i)
                if (dist < bestExactDist) {
                    bestExactDist = dist
                    bestExactPos = i
                }
            }
            method.equals("approx_moment_space", ignoreCase = true) -> {
                apxMSEntries.add(i)
                if (dist < bestApxMSDist) {
                    bestApxMSDist = dist
                    bestApxMSPos = i
                }
            }
            method.equals("approx_param_space", ignoreCase = true) -> {
                apxPSEntries.add(i)
                if (dist < bestApxPSDist) {
                    bestApxPSDist = dist
                    bestApxPSPos = i
                }
            }
        }
    }

    // Header
    sb.append("Log10\tTRACE        ")
    for (j in exactEntries.indices) {
        sb.append("\tEXMM (idx=${j + 1})   ")
    }
    for (j in apxMSEntries.indices) {
        sb.append("\tAPXMS (idx=${exactEntries.size + j + 1})  ")
    }
    for (j in apxPSEntries.indices) {
        sb.append("\tAPXPS (idx=${exactEntries.size + apxMSEntries.size + j + 1})  ")
    }
    sb.append("\n")

    // Moment rows
    for (k in 1..E.size) {
        if (k == 1) {
            sb.append("E[X]")
        } else {
            sb.append("E[X^$k]")
        }
        sb.append("\t${String.format("%13.6f", FastMath.log10(E[k - 1]))}")

        for (j in exactEntries) {
            val phMAP = PH[j][0] as MatrixCell
            val mom = map_moment(phMAP, k)
            sb.append("\t${String.format("%13.6f", FastMath.log10(mom))}")
        }
        for (j in apxMSEntries) {
            val phMAP = PH[j][0] as MatrixCell
            val mom = map_moment(phMAP, k)
            sb.append("\t${String.format("%13.6f", FastMath.log10(mom + Double.MIN_VALUE))}")
        }
        for (j in apxPSEntries) {
            val phMAP = PH[j][0] as MatrixCell
            val mom = map_moment(phMAP, k)
            sb.append("\t${String.format("%13.6f", FastMath.log10(mom + Double.MIN_VALUE))}")
        }
        sb.append("\n")
    }

    // Summary
    sb.append("\n")
    if (bestExactPos >= 0) {
        sb.append("Best exact moment match: idx=${bestExactPos + 1}, dist=$bestExactDist, ")
        val phMAP = PH[bestExactPos][0] as MatrixCell
        sb.append("order=${phMAP[0].numRows}\n")
    }
    if (bestApxMSPos >= 0) {
        sb.append("Best approx moment space: idx=${bestApxMSPos + 1}, dist=$bestApxMSDist, ")
        val phMAP = PH[bestApxMSPos][0] as MatrixCell
        sb.append("order=${phMAP[0].numRows}\n")
    }
    if (bestApxPSPos >= 0) {
        sb.append("Best approx param space: idx=${bestApxPSPos + 1}, dist=$bestApxPSDist, ")
        val phMAP = PH[bestApxPSPos][0] as MatrixCell
        sb.append("order=${phMAP[0].numRows}\n")
    }

    return sb.toString()
}

/**
 * Creates an exponential MAP with given mean.
 * Convenience wrapper around map_exponential.
 */
private fun map_exponential(mean: Double): MatrixCell {
    return jline.api.mam.map_exponential(mean)
}
