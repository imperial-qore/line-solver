/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
/**
 * @file Single trace variance computation
 * 
 * Computes sample variance and related statistics for empirical trace data. 
 * Essential for characterizing variability in arrival/service processes and 
 * parameterizing stochastic models from measurement data.
 * 
 * @since LINE 3.0
 */
package jline.api.trace

import jline.util.matrix.Matrix
import kotlin.collections.ArrayList
import kotlin.math.*
import kotlin.random.Random

/**
 * APIs for statistical analysis of trace data.
 */

/**
 * Computes the variance of the trace data.
 *
 * @param trace the array containing the trace data
 * @return the variance of the trace data
 */
fun trace_var(trace: DoubleArray): Double {
    var e1 = 0.0
    var e2 = 0.0
    for (i in trace.indices) {
        e2 += trace[i] * trace[i]
        e1 += trace[i]
    }
    return e2 / trace.size - (e1 / trace.size) * (e1 / trace.size)
}

/**
 * Computes the squared coefficient of variation for trace S
 *
 * @param S the trace data
 * @return the squared coefficient of variation
 */
fun trace_scv(S: DoubleArray): Double {
    return trace_var(S) / (trace_mean(S) * trace_mean(S))
}

/**
 * Computes the autocorrelation function for trace S at the specified lags
 *
 * @param S the trace data  
 * @param lags the lag values to compute (default is arrayOf(1))
 * @return the autocorrelation values at the specified lags
 */
fun trace_acf(S: DoubleArray, lags: IntArray = intArrayOf(1)): DoubleArray {
    val maxLag = lags.maxOrNull() ?: 1
    
    if (maxLag > S.size - 2) {
        // Truncate lags that are too large
        val filteredLags = lags.filter { it <= S.size - 2 && it > 0 }.toIntArray()
        if (filteredLags.isEmpty()) {
            return doubleArrayOf()
        }
        return trace_acf(S, filteredLags)
    }
    
    val mean = trace_mean(S)
    val centeredS = S.map { it - mean }.toDoubleArray()
    
    // Compute autocovariance using FFT-based method similar to autocov.m
    val autocov = autocov(centeredS)
    val rho = DoubleArray(lags.size)
    
    for (i in lags.indices) {
        rho[i] = autocov[lags[i]] / autocov[0]
    }
    
    return rho
}

/**
 * Computes autocovariance using FFT-based method (Java 8 compatible)
 */
private fun autocov(X: DoubleArray): DoubleArray {
    val M = X.size
    val acv = DoubleArray(M - 1)
    
    // Simple direct computation for Java 8 compatibility
    for (p in 0 until M - 1) {
        var sum = 0.0
        var count = 0
        for (i in 0 until M - p) {
            sum += X[i] * X[i + p]
            count++
        }
        acv[p] = sum / count
    }
    
    return acv
}

/**
 * Estimates the auto-correlation decay rate of a trace
 *
 * @param T the trace data
 * @param limit maximum lag considered (default 1000)
 * @return array containing [GAMMA, RHO0, RESIDUALS]
 */
fun trace_gamma(T: DoubleArray, limit: Int = 1000): DoubleArray {
    val M1 = trace_mean(T)
    val M2 = T.map { it * it }.average()
    
    val maxLag = minOf(limit, T.size - 1)
    val lag = (1..maxLag).toList().toIntArray()
    val rho = trace_acf(T, lag)
    
    val VAR = M2 - M1 * M1
    val SCV = VAR / (M1 * M1)
    val RHO0 = 0.5 * (1.0 - 1.0 / SCV)
    
    // Simple exponential fit - in real implementation would use nonlinear fitting
    var bestGamma = 0.99
    var minResiduals = Double.MAX_VALUE
    
    // Grid search for gamma
    for (gamma in 990..999) {
        val g = gamma / 1000.0
        var residuals = 0.0
        for (i in lag.indices) {
            val expected = RHO0 * g.pow(lag[i])
            residuals += (rho[i] - expected) * (rho[i] - expected)
        }
        if (residuals < minResiduals) {
            minResiduals = residuals
            bestGamma = g
        }
    }
    
    return doubleArrayOf(bestGamma, RHO0, minResiduals)
}

/**
 * Computes the counting process of S
 *
 * @param S the inter-arrival times
 * @param scale the time scale
 * @return the counts after "scale" units of time from an arrival
 */
fun trace_iat2counts(S: DoubleArray, scale: Double): IntArray {
    val n = S.size
    val CS = DoubleArray(n + 1) // cumulative sum with 0 at start
    CS[0] = 0.0
    for (i in 1..n) {
        CS[i] = CS[i-1] + S[i-1]
    }
    
    val C = ArrayList<Int>()
    
    for (i in 0 until n-1) {
        var cur = i
        while (cur + 1 < n && CS[cur + 1] - CS[i] <= scale) {
            cur++
        }
        C.add(cur - i)
        if (cur == n - 1) {
            break
        }
    }
    
    return C.toIntArray()
}

/**
 * Computes the Index of Dispersion for Intervals of a trace
 *
 * @param S the trace data
 * @param kset the set of k values to compute IDI for
 * @param option aggregation option (null, "aggregate", "aggregate-mix")  
 * @param n aggregation parameter
 * @return pair of (IDI values, support values)
 */
fun trace_idi(S: DoubleArray, kset: IntArray, option: String? = null, n: Int = 1): Pair<DoubleArray, IntArray> {
    val IDIk = ArrayList<Double>()
    val support = ArrayList<Int>()
    
    for (k in kset) {
        when (option) {
            null -> {
                val supportVal = S.size - k - 1
                support.add(supportVal)
                val Sk = DoubleArray(S.size - k)
                for (t in 0 until S.size - k - 1) {
                    var sum = 0.0
                    for (j in t until t + k) {
                        sum += S[j]
                    }
                    Sk[t] = sum
                }
                val variance = trace_var(Sk)
                val mean = trace_mean(Sk)
                IDIk.add(k * variance / (mean * mean))
            }
            "aggregate" -> {
                val keff = k / n
                val supportVal = S.size / keff
                support.add(supportVal)
                val Sk = DoubleArray(S.size - keff)
                for (t in 0 until S.size - keff - 1) {
                    var sum = 0.0
                    for (j in t until t + keff) {
                        sum += S[j]
                    }
                    Sk[t] = sum
                }
                val variance = trace_var(Sk)
                val mean = trace_mean(Sk)
                IDIk.add(k * variance / (mean * mean))
            }
            // Simplified version of "aggregate-mix" - would need more complex implementation
        }
    }
    
    return Pair(IDIk.toDoubleArray(), support.toIntArray())
}

/**
 * Computes the Index of Dispersion for Counts (asymptotically equal to IDI)
 *
 * @param S the trace data
 * @return the IDC value
 */
fun trace_idc(S: DoubleArray): Double {
    val kset = intArrayOf(minOf(1000, S.size / 30))
    val result = trace_idi(S, kset)
    return result.first[0]
}

/**
 * Computes the probability mass function of discrete data
 *
 * @param X the discrete data
 * @return pair of (pmf values, unique values)
 */
fun trace_pmf(X: IntArray): Pair<DoubleArray, IntArray> {
    val uniqueValues = X.distinct().sorted()
    val pmf = DoubleArray(uniqueValues.size)
    
    for (i in uniqueValues.indices) {
        val count = X.count { it == uniqueValues[i] }
        pmf[i] = count.toDouble() / X.size
    }
    
    return Pair(pmf, uniqueValues.toIntArray())
}

/**
 * Shuffles the trace data randomly
 *
 * @param S the trace data
 * @return the shuffled trace
 */
fun trace_shuffle(S: DoubleArray): DoubleArray {
    val result = S.copyOf()
    val random = Random.Default
    for (i in result.indices.reversed()) {
        val j = random.nextInt(i + 1)
        val temp = result[i]
        result[i] = result[j]
        result[j] = temp
    }
    return result
}

/**
 * Computes the bicovariance of the trace
 *
 * @param S the trace data
 * @param GRID the grid of lag values
 * @return pair of (bicovariance values, lag combinations)
 */
fun trace_bicov(S: DoubleArray, GRID: IntArray): Pair<DoubleArray, Array<IntArray>> {
    val BiCovLags = ArrayList<IntArray>()
    val BiCov = ArrayList<Double>()
    
    for (i in GRID) {
        for (j in GRID) {
            BiCovLags.add(intArrayOf(1, i, j))
        }
    }
    
    for (lags in BiCovLags) {
        val jointMoment = trace_joint(S, lags, intArrayOf(1, 1, 1))
        BiCov.add(jointMoment)
    }
    
    return Pair(BiCov.toDoubleArray(), BiCovLags.toTypedArray())
}

/**
 * Computes the counts in each bin with specified timescale
 *
 * @param S the inter-arrival times
 * @param scale the bin timescale
 * @return pair of (counts per bin, bin membership for each element)
 */
fun trace_iat2bins(S: DoubleArray, scale: Double): Pair<IntArray, IntArray> {
    val n = S.size
    val CS = DoubleArray(n + 1)
    CS[0] = 0.0
    for (i in 1..n) {
        CS[i] = CS[i-1] + S[i-1]
    }
    
    val bins = ceil((CS[n] - CS[0]) / scale).toInt()
    val C = IntArray(bins)
    val bC = ArrayList<Int>()
    
    var cur = 0
    var last = 0
    
    for (i in 0 until bins) {
        if (cur >= n - 1) break
        
        while (cur < n - 1 && CS[cur + 1] <= (i + 1) * scale) {
            cur++
        }
        C[i] = cur - last
        for (j in 0 until cur - last) {
            bC.add(i)
        }
        last = cur
    }
    
    return Pair(C, bC.toIntArray())
}

/**
 * Computes joint moments E[X^{k_1}_{i} X^{k_2}_{i+j} ...] for a trace
 *
 * @param S the trace data
 * @param lag the lag values (cumulative)
 * @param order the moment orders
 * @return the joint moment value
 */
fun trace_joint(S: DoubleArray, lag: IntArray, order: IntArray): Double {
    val sortedLag = lag.sorted().toIntArray()
    val K = sortedLag.size
    val adjustedLag = IntArray(K)
    val baseLag = sortedLag[0]
    for (i in 0 until K) {
        adjustedLag[i] = sortedLag[i] - baseLag
    }
    
    val maxLag = adjustedLag.maxOrNull() ?: 0
    val validLength = S.size - maxLag
    
    if (validLength <= 0) return 0.0
    
    var sum = 0.0
    for (i in 0 until validLength) {
        var product = 1.0
        for (j in 0 until order.size) {
            val idx = i + adjustedLag[minOf(j, adjustedLag.size - 1)]
            if (idx < S.size) {
                product *= S[idx].pow(order[j])
            }
        }
        sum += product
    }
    
    return sum / validLength
}

/**
 * Computes comprehensive summary statistics for a trace
 *
 * @param m the trace data
 * @return array containing [MEAN,SCV,MAD,SKEW,KURT,Q25,Q50,Q75,P95,MIN,MAX,IQR,ACF1-4,IDC_SCV_RATIO]
 */
fun trace_summary(m: DoubleArray): DoubleArray {
    val mean = trace_mean(m)
    val scv = trace_scv(m)
    val sortedM = m.sorted().toDoubleArray()
    
    // Compute percentiles
    val q25 = percentile(sortedM, 25.0)
    val q50 = percentile(sortedM, 50.0) // median
    val q75 = percentile(sortedM, 75.0)
    val p95 = percentile(sortedM, 95.0)
    val min = sortedM[0]
    val max = sortedM[sortedM.size - 1]
    val iqr = q75 - q25
    
    // MAD (median absolute deviation)
    val mad = m.map { abs(it - q50) }.sorted()[m.size / 2]
    
    // Skewness and kurtosis (simplified)
    val variance = trace_var(m)
    val std = sqrt(variance)
    var skew = 0.0
    var kurt = 0.0
    
    for (value in m) {
        val z = (value - mean) / std
        skew += z * z * z
        kurt += z * z * z * z
    }
    skew /= m.size
    kurt = kurt / m.size - 3.0 // excess kurtosis
    
    // ACF for lags 1-4
    val acf = trace_acf(m, intArrayOf(1, 2, 3, 4))
    
    // IDC/SCV ratio
    val idc = trace_idc(m)
    val idcScvRatio = idc / scv
    
    return doubleArrayOf(mean, scv, mad, skew, kurt, q25, q50, q75, p95, min, max, iqr, 
                        acf[0], acf[1], acf[2], acf[3], idcScvRatio)
}

/**
 * Computes percentile of sorted data
 */
private fun percentile(sortedData: DoubleArray, p: Double): Double {
    val index = (p / 100.0) * (sortedData.size - 1)
    val lower = floor(index).toInt()
    val upper = ceil(index).toInt()
    
    if (lower == upper) {
        return sortedData[lower]
    }
    
    val weight = index - lower
    return sortedData[lower] * (1 - weight) + sortedData[upper] * weight
}
/**
 * Trace Var algorithms
 */
@Suppress("unused")
class TraceVarAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}