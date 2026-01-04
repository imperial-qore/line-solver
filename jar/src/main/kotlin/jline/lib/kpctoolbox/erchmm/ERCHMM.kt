package jline.lib.kpctoolbox.erchmm

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import org.apache.commons.math3.linear.LUDecomposition
import org.apache.commons.math3.linear.MatrixUtils
import org.apache.commons.math3.util.FastMath
import kotlin.math.abs
import kotlin.math.ln

/**
 * Extended Renewal Continuous-time Hidden Markov Model (ER-CHMM) functions.
 *
 * The EM algorithm is based on:
 * - Okamura et al. (2008). An EM algorithm for a Superposition of Markovian Arrival Processes.
 * - Horvath & Okamura (2013). A Fast EM algorithm for Fitting Marked Markov Arrival Processes.
 * - Horvath et al. (2018). Parallel Algorithms for Fitting Markov Arrival Processes.
 *
 * Ported from MATLAB: matlab/lib/kpctoolbox/erchmm/
 */

/**
 * Result class for ER-CHMM EM fitting.
 */
data class ERCHMMFitResult(
    val MAP: MatrixCell,
    val logLikelihood: Double,
    val orders: IntArray
)

/**
 * Fits an Extended Renewal Continuous-time Hidden Markov Model to a trace using EM algorithm.
 *
 * @param trace Array of inter-arrival times
 * @param orders Either a single integer (sum of orders) or array of Erlang branch orders
 * @param iterMax Maximum number of EM iterations (default 300)
 * @param iterTol Convergence tolerance (default 1e-7)
 * @param verbose Whether to print progress (default false)
 * @return ERCHMMFitResult containing fitted MAP, log-likelihood, and best orders
 */
fun erchmm_emfit(
    trace: DoubleArray,
    orders: IntArray,
    iterMax: Int = 300,
    iterTol: Double = 1e-7,
    verbose: Boolean = false
): ERCHMMFitResult {

    // If orders is a single value, try all combinations of Erlang branch orders
    if (orders.size == 1) {
        val totalOrder = orders[0]
        var bestMAP: MatrixCell? = null
        var bestLogLi = Double.NEGATIVE_INFINITY
        var bestOrders: IntArray = intArrayOf()

        // Try different numbers of branches (2 to totalOrder)
        for (numBranches in 2..totalOrder) {
            val allCombinations = allErlangCombinations(numBranches, totalOrder)

            for (combination in allCombinations) {
                if (verbose) {
                    println("Calculating with orders ${combination.joinToString(",")}")
                }

                val result = erchmm_emfit(trace, combination, iterMax, iterTol, verbose)

                if (result.logLikelihood > bestLogLi) {
                    bestMAP = result.MAP
                    bestLogLi = result.logLikelihood
                    bestOrders = combination
                }
            }
        }

        if (verbose) {
            println("Best solution: log-likelihood=$bestLogLi, orders=${bestOrders.joinToString(",")}")
        }

        return ERCHMMFitResult(
            bestMAP ?: MatrixCell(2),
            bestLogLi,
            bestOrders
        )
    }

    // Main EM algorithm for given orders
    val M = orders.size
    val K = trace.size

    // Initialize pi and lambda parameters such that the mean is matched
    val piV = DoubleArray(M) { 1.0 / M }
    val lambda = DoubleArray(M) { i -> orders[i].toDouble() * (i + 1) }

    val traceMean = trace.sum() / K
    var piMean = 0.0
    for (i in 0 until M) {
        piMean += piV[i] / (i + 1)
    }

    // Match the mean
    for (i in 0 until M) {
        lambda[i] = lambda[i] * piMean / traceMean
    }

    // Initialize the transition matrix
    val T = Array(M) { piV.copyOf() }

    // Initialize matrices for EM algorithm
    val F = Array(M) { DoubleArray(K) }  // Matrix for densities
    val aLikelihoods = Array(K) { DoubleArray(M) }  // Forward likelihoods
    val bLikelihoods = Array(M) { DoubleArray(K) }  // Backward likelihoods
    val aLikelihoodsScale = DoubleArray(K)
    val bLikelihoodsScale = DoubleArray(K)

    var ologli = 1.0
    var logL = 0.0
    var steps = 1

    // Main EM iteration loop
    while (abs((ologli - logL) / logL) > iterTol && steps < iterMax) {
        ologli = logL

        // E-step: Compute branch densities
        for (i in 0 until M) {
            for (k in 0 until K) {
                val lt = lambda[i] * trace[k]
                val order = orders[i]
                F[i][k] = (FastMath.pow(lt, (order - 1).toDouble()) / factorial(order - 1) * lambda[i]) *
                        FastMath.exp(-lambda[i] * trace[k])
            }
        }

        // Compute forward likelihood vectors
        val prevPi = piV.copyOf()
        var scaledPrev = 0.0

        for (k in 0 until K) {
            // prevPi = prevPi * diag(F(:,k)) * T
            val temp = DoubleArray(M)
            for (j in 0 until M) {
                var sum = 0.0
                for (i in 0 until M) {
                    sum += prevPi[i] * F[i][k] * T[i][j]
                }
                temp[j] = sum
            }

            val scale = log2(temp.sum())
            val scaleFactor = FastMath.pow(2.0, -scale)
            for (j in 0 until M) {
                temp[j] *= scaleFactor
                prevPi[j] = temp[j]
            }

            aLikelihoodsScale[k] = scaledPrev + scale
            aLikelihoods[k] = temp.copyOf()
            scaledPrev = aLikelihoodsScale[k]
        }

        // Forward likelihoods shifted
        val aForwardLikelihoods = Array(K) { k ->
            if (k == 0) piV.copyOf() else aLikelihoods[k - 1].copyOf()
        }
        val aScaleV = DoubleArray(K) { k ->
            if (k == 0) 0.0 else aLikelihoodsScale[k - 1]
        }

        // Compute backward likelihood vectors
        val nextB = DoubleArray(M) { 1.0 }
        scaledPrev = 0.0

        for (k in K - 1 downTo 0) {
            // nextB = diag(F(:,k)) * T * nextB
            val temp = DoubleArray(M)
            for (i in 0 until M) {
                var sum = 0.0
                for (j in 0 until M) {
                    sum += F[i][k] * T[i][j] * nextB[j]
                }
                temp[i] = sum
            }

            val scale = log2(temp.sum())
            val scaleFactor = FastMath.pow(2.0, -scale)
            for (i in 0 until M) {
                temp[i] *= scaleFactor
                nextB[i] = temp[i]
            }

            bLikelihoodsScale[k] = scaledPrev + scale
            for (i in 0 until M) {
                bLikelihoods[i][k] = nextB[i]
            }
            scaledPrev = bLikelihoodsScale[k]
        }

        // Backward likelihoods shifted
        val bBackwardLikelihoods = Array(M) { i ->
            DoubleArray(K) { k ->
                if (k == K - 1) 1.0 else bLikelihoods[i][k + 1]
            }
        }
        val bScaleV = DoubleArray(K) { k ->
            if (k == K - 1) 0.0 else bLikelihoodsScale[k + 1]
        }

        // Calculate likelihood
        var likelihoodValue = 0.0
        for (i in 0 until M) {
            likelihoodValue += piV[i] * bLikelihoods[i][0]
        }

        logL = (ln(likelihoodValue) + bLikelihoodsScale[0] * ln(2.0)) / K
        val iLikelihood = 1.0 / likelihoodValue

        // M-step: Update parameters

        // Calculate likelihoods multiplied: a[k] .* b[k]
        val likelihoodsMultiplied = Array(K) { k ->
            DoubleArray(M) { m -> aForwardLikelihoods[k][m] * bLikelihoods[m][k] }
        }

        // Normalize each row
        for (k in 0 until K) {
            val sumRow = likelihoodsMultiplied[k].sum()
            if (sumRow > 0) {
                for (m in 0 until M) {
                    likelihoodsMultiplied[k][m] /= sumRow
                }
            }
        }

        // Compute numerator and denominator for lambda estimations
        val numeratorEstimation = DoubleArray(M)
        val denominatorEstimation = DoubleArray(M)

        for (m in 0 until M) {
            for (k in 0 until K) {
                numeratorEstimation[m] += likelihoodsMultiplied[k][m]
                denominatorEstimation[m] += trace[k] * likelihoodsMultiplied[k][m]
            }
        }

        // Update pi and lambda
        for (m in 0 until M) {
            piV[m] = numeratorEstimation[m] / K
            if (denominatorEstimation[m] > 0) {
                lambda[m] = orders[m] * numeratorEstimation[m] / denominatorEstimation[m]
            }
        }

        // Compute multiplication of forward likelihood and densities for T estimation
        val densMutALikelihood = Array(K) { k ->
            DoubleArray(M) { m -> aForwardLikelihoods[k][m] * F[m][k] }
        }

        // Scale factors
        val summedLm = DoubleArray(K) { k ->
            iLikelihood * FastMath.pow(2.0, aScaleV[k] + bScaleV[k] - bLikelihoodsScale[0])
        }

        for (k in 0 until K) {
            for (m in 0 until M) {
                densMutALikelihood[k][m] *= summedLm[k]
            }
        }

        // Update T: T = (densMutALikelihood' * bBackwardLikelihoods') .* T
        for (i in 0 until M) {
            for (j in 0 until M) {
                var sum = 0.0
                for (k in 0 until K) {
                    sum += densMutALikelihood[k][i] * bBackwardLikelihoods[j][k]
                }
                T[i][j] = sum * T[i][j]
            }
        }

        // Normalize T rows
        for (i in 0 until M) {
            val rowSum = T[i].sum()
            if (rowSum > 0) {
                for (j in 0 until M) {
                    T[i][j] /= rowSum
                }
            }
        }

        steps++

        if (verbose && steps % 50 == 0) {
            println("Num of iterations: $steps, log-likelihood: $logL")
        }
    }

    if (verbose) {
        println("Num of iterations: $steps, log-likelihood: $logL")
        println("EM algorithm terminated. (orders=${orders.joinToString(",")})")
    }

    // Finalize D0 and D1 matrices
    val D0 = generateD0FromErlangs(lambda, orders)
    val D1 = generateD1FromErlangs(lambda, orders, T)

    val MAP = MatrixCell(2)
    MAP[0] = D0
    MAP[1] = D1

    return ERCHMMFitResult(MAP, logL, orders)
}

/**
 * Generates all unique combinations of Erlang branch orders that sum to a given total.
 */
private fun allErlangCombinations(branches: Int, sumErlangs: Int): List<IntArray> {
    if (branches == 1) {
        return listOf(intArrayOf(sumErlangs))
    }

    val result = mutableListOf<IntArray>()

    for (k1 in 1..sumErlangs - branches + 1) {
        val subCombinations = allErlangCombinations(branches - 1, sumErlangs - k1)

        for (sub in subCombinations) {
            val combined = (sub.toList() + k1).sorted().toIntArray()

            // Check if this combination already exists
            val exists = result.any { existing ->
                existing.contentEquals(combined)
            }

            if (!exists) {
                result.add(combined)
            }
        }
    }

    return result
}

/**
 * Generates D0 matrix from Erlang representation.
 */
private fun generateD0FromErlangs(lambda: DoubleArray, orders: IntArray): Matrix {
    val n = orders.sum()
    val D0 = Matrix(n, n)

    var startIdx = 0
    for (i in lambda.indices) {
        val order = orders[i]
        val lam = lambda[i]

        for (j in 0 until order) {
            D0.set(startIdx + j, startIdx + j, -lam)
            if (j < order - 1) {
                D0.set(startIdx + j, startIdx + j + 1, lam)
            }
        }

        startIdx += order
    }

    return D0
}

/**
 * Generates D1 matrix from Erlang representation.
 */
private fun generateD1FromErlangs(lambda: DoubleArray, orders: IntArray, T: Array<DoubleArray>): Matrix {
    val n = orders.sum()
    val M = orders.size
    val D1 = Matrix(n, n)

    // Compute indices
    val indicesTo = IntArray(M)
    val indicesFrom = IntArray(M)

    var cumSum = 0
    for (i in 0 until M) {
        indicesTo[i] = cumSum
        cumSum += orders[i]
        indicesFrom[i] = cumSum - 1
    }

    // D1(indicesFrom, indicesTo) = diag(lambda) * T
    for (i in 0 until M) {
        for (j in 0 until M) {
            D1.set(indicesFrom[i], indicesTo[j], lambda[i] * T[i][j])
        }
    }

    return D1
}

/**
 * Computes factorial.
 */
private fun factorial(n: Int): Double {
    if (n <= 1) return 1.0
    var result = 1.0
    for (i in 2..n) {
        result *= i.toDouble()
    }
    return result
}

/**
 * Computes log base 2.
 */
private fun log2(x: Double): Double {
    return if (x > 0) ln(x) / ln(2.0) else 0.0
}

/**
 * Simplified version of erchmm_emfit that returns only the MAP.
 */
fun erchmm_emfit_simple(
    trace: DoubleArray,
    orders: IntArray,
    iterMax: Int = 300,
    iterTol: Double = 1e-7
): MatrixCell {
    return erchmm_emfit(trace, orders, iterMax, iterTol, false).MAP
}
