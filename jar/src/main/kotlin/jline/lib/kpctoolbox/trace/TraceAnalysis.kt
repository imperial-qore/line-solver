package jline.lib.kpctoolbox.trace

import org.apache.commons.math3.stat.StatUtils
import org.apache.commons.math3.stat.descriptive.moment.Kurtosis
import org.apache.commons.math3.stat.descriptive.moment.Skewness
import org.apache.commons.math3.stat.descriptive.rank.Percentile
import org.apache.commons.math3.transform.DftNormalization
import org.apache.commons.math3.transform.FastFourierTransformer
import org.apache.commons.math3.transform.TransformType
import org.apache.commons.math3.util.FastMath
import java.util.Random
import kotlin.math.ceil

/**
 * Trace analysis functions for the KPC-Toolbox.
 *
 * Ported from MATLAB: matlab/lib/kpctoolbox/trace/
 */

/**
 * Computes the mean of a trace.
 *
 * @param S Input trace
 * @return Mean value
 */
fun trace_mean(S: DoubleArray): Double {
    return StatUtils.mean(S)
}

/**
 * Computes the variance of a trace.
 *
 * @param S Input trace
 * @return Variance
 */
fun trace_var(S: DoubleArray): Double {
    return StatUtils.variance(S)
}

/**
 * Computes the squared coefficient of variation (SCV) of a trace.
 *
 * @param S Input trace
 * @return Squared coefficient of variation (variance / mean^2)
 */
fun trace_scv(S: DoubleArray): Double {
    val mean = trace_mean(S)
    val variance = trace_var(S)
    return variance / (mean * mean)
}

/**
 * Computes the autocovariance of a trace using FFT.
 *
 * @param X Input trace (column vector conceptually)
 * @return Autocovariance vector
 */
fun autocov(X: DoubleArray): DoubleArray {
    val M = X.size
    if (M < 2) {
        throw IllegalArgumentException("Trace too short")
    }

    val mean = StatUtils.mean(X)

    // Pad X - mean(X) with zeros
    val padLen = nextPowerOf2(2 * M)
    val xPad = DoubleArray(padLen)
    for (i in 0 until M) {
        xPad[i] = X[i] - mean
    }

    // FFT
    val fft = FastFourierTransformer(DftNormalization.STANDARD)
    val xHat = fft.transform(xPad, TransformType.FORWARD)

    // Compute |xHat|^2 (autoconvolution in frequency domain)
    val product = DoubleArray(padLen)
    for (i in 0 until padLen) {
        product[i] = xHat[i].real * xHat[i].real + xHat[i].imaginary * xHat[i].imaginary
    }

    // Inverse FFT
    val result = fft.transform(
        product.map { org.apache.commons.math3.complex.Complex(it, 0.0) }.toTypedArray(),
        TransformType.INVERSE
    )

    // Extract and normalize autocovariance
    val acv = DoubleArray(M - 1)
    for (i in 0 until M - 1) {
        acv[i] = result[i].real / (M - i)
    }

    return acv
}

/**
 * Helper function to find next power of 2.
 */
private fun nextPowerOf2(n: Int): Int {
    var power = 1
    while (power < n) {
        power *= 2
    }
    return power
}

/**
 * Computes the autocorrelation function (ACF) for a trace at specified lags.
 *
 * @param S Input trace
 * @param lags Lags at which to compute ACF (1-based in MATLAB, converted to 0-based internally)
 * @return Autocorrelation values at specified lags
 */
fun trace_acf(S: DoubleArray, lags: IntArray): DoubleArray {
    if (lags.isEmpty()) {
        return doubleArrayOf()
    }

    val n = S.size
    val maxLag = lags.maxOrNull() ?: 1

    // Validate lags
    if (maxLag > n - 2) {
        // Truncate lags
        val validLags = lags.filter { it <= n - 2 }.toIntArray()
        if (validLags.isEmpty()) {
            return doubleArrayOf()
        }
        return trace_acf(S, validLags)
    }

    // Compute autocovariance using FFT
    val acv = autocov(S)

    // Compute ACF: rho(k) = acv(k) / acv(0)
    val rho = DoubleArray(lags.size)
    val variance = acv[0]

    if (variance == 0.0) {
        return DoubleArray(lags.size) { 0.0 }
    }

    for ((idx, lag) in lags.withIndex()) {
        if (lag > 0 && lag < acv.size) {
            rho[idx] = acv[lag] / variance
        } else if (lag == 0) {
            rho[idx] = 1.0
        }
    }

    return rho
}

/**
 * Computes the autocorrelation function at lag 1.
 *
 * @param S Input trace
 * @return Lag-1 autocorrelation
 */
fun trace_acf(S: DoubleArray): Double {
    return trace_acf(S, intArrayOf(1))[0]
}

/**
 * Computes the bias-corrected skewness of a trace.
 *
 * @param S Input trace
 * @return Skewness value
 */
fun trace_skew(S: DoubleArray): Double {
    val filtered = S.filter { !it.isNaN() }.toDoubleArray()
    val n = filtered.size

    if (n < 3) {
        return Double.NaN
    }

    val mean = StatUtils.mean(filtered)
    var s2 = 0.0
    var m3 = 0.0

    for (x in filtered) {
        val res = x - mean
        s2 += res * res
        m3 += res * res * res
    }

    s2 /= n
    m3 /= n

    // Bias correction factor
    val correction = FastMath.sqrt((n - 1.0) / n) * n / (n - 2.0)

    return m3 * correction / FastMath.pow(s2, 1.5)
}

/**
 * Computes the joint moment E[X^{k_1}_{i} X^{k_2}_{i+lag_1} ...] for a trace.
 *
 * @param S Input trace
 * @param lags Array of lag values (cumulative)
 * @param orders Array of moment orders for each position
 * @return Joint moment value
 */
fun trace_joint(S: DoubleArray, lags: IntArray, orders: IntArray): Double {
    if (lags.size != orders.size) {
        throw IllegalArgumentException("lags and orders must have same length")
    }

    // Compute cumulative lags and normalize
    val cumLags = IntArray(lags.size)
    var sum = 0
    for (i in lags.indices) {
        sum += lags[i]
        cumLags[i] = sum
    }

    // Shift to start from 0
    val minLag = cumLags.minOrNull() ?: 0
    for (i in cumLags.indices) {
        cumLags[i] -= minLag
    }

    val maxLag = cumLags.maxOrNull() ?: 0
    val n = S.size

    if (maxLag >= n) {
        throw IllegalArgumentException("Lag exceeds trace length")
    }

    // Compute joint moment
    var total = 0.0
    val validCount = n - maxLag

    for (t in 0 until validCount) {
        var product = 1.0
        for (i in cumLags.indices) {
            product *= FastMath.pow(S[t + cumLags[i]], orders[i].toDouble())
        }
        total += product
    }

    return total / validCount
}

/**
 * Computes the bicovariance of a trace for a grid of lags.
 *
 * @param S Input trace
 * @param grid Lag values for bicovariance grid
 * @return Pair of (bicovariance values, lag pairs)
 */
fun trace_bicov(S: DoubleArray, grid: IntArray): Pair<DoubleArray, Array<IntArray>> {
    val lagPairs = ArrayList<IntArray>()
    for (i in grid) {
        for (j in grid) {
            lagPairs.add(intArrayOf(1, i, j))
        }
    }

    val bicov = DoubleArray(lagPairs.size)
    for ((idx, lagPair) in lagPairs.withIndex()) {
        bicov[idx] = trace_joint(S, lagPair, intArrayOf(1, 1, 1))
    }

    return Pair(bicov, lagPairs.toTypedArray())
}

/**
 * Computes the index of dispersion for intervals (IDI) of a trace.
 *
 * @param S Input trace (inter-arrival times)
 * @param kset Array of k values at which to compute IDI
 * @return Array of IDI values
 */
fun trace_idi(S: DoubleArray, kset: IntArray): DoubleArray {
    val idiValues = DoubleArray(kset.size)

    for ((idx, k) in kset.withIndex()) {
        if (k >= S.size) {
            idiValues[idx] = Double.NaN
            continue
        }

        // Compute sums of k consecutive values
        val numSums = S.size - k
        if (numSums <= 0) {
            idiValues[idx] = Double.NaN
            continue
        }

        val Sk = DoubleArray(numSums)
        for (t in 0 until numSums) {
            var sum = 0.0
            for (i in 0 until k) {
                sum += S[t + i]
            }
            Sk[t] = sum
        }

        val meanSk = StatUtils.mean(Sk)
        val varSk = StatUtils.variance(Sk)

        if (meanSk != 0.0) {
            idiValues[idx] = k * varSk / (meanSk * meanSk)
        } else {
            idiValues[idx] = Double.NaN
        }
    }

    return idiValues
}

/**
 * Computes the index of dispersion for intervals with default k values.
 *
 * @param S Input trace
 * @return IDI value at k = min(1000, trace_length/30)
 */
fun trace_idi(S: DoubleArray): Double {
    val k = minOf(1000, S.size / 30)
    if (k <= 0) return Double.NaN
    return trace_idi(S, intArrayOf(k))[0]
}

/**
 * Computes the index of dispersion for counts (IDC) of a trace.
 * Asymptotically, IDI and IDC are equal.
 *
 * @param S Input trace
 * @return IDC value
 */
fun trace_idc(S: DoubleArray): Double {
    return trace_idi(S)
}

/**
 * Estimates the autocorrelation decay rate of a trace.
 *
 * @param T Input trace
 * @param limit Maximum lag considered (default: 1000)
 * @return Autocorrelation decay rate gamma
 */
fun trace_gamma(T: DoubleArray, limit: Int = 1000): Double {
    val M1 = trace_mean(T)
    val M2 = StatUtils.mean(T.map { it * it }.toDoubleArray())

    val maxLag = minOf(limit, T.size - 2)
    if (maxLag <= 0) return Double.NaN

    val lags = (1..maxLag).toList().toIntArray()
    val rho = trace_acf(T, lags)

    val variance = M2 - M1 * M1
    val scv = variance / (M1 * M1)
    val rho0 = 0.5 * (1 - 1 / scv)

    // Simple least squares fit for geometric decay: rho(k) = rho0 * gamma^k
    // Taking log: log(rho(k)) = log(rho0) + k * log(gamma)
    // We solve for gamma by minimizing sum of squared residuals

    // Filter positive rho values for log fitting
    val validPairs = lags.zip(rho.toList())
        .filter { it.second > 0 && it.second < 1 }

    if (validPairs.isEmpty()) {
        return 0.5 // Default fallback
    }

    // Use simple averaging of log(rho_k / rho_0) / k
    var sumLogRatio = 0.0
    var count = 0
    for ((k, rhoK) in validPairs) {
        if (rho0 > 0 && rhoK / rho0 > 0) {
            sumLogRatio += FastMath.log(rhoK / rho0) / k
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
 * Shuffles a trace randomly.
 *
 * @param S Input trace
 * @return Shuffled trace
 */
fun trace_shuffle(S: DoubleArray): DoubleArray {
    val result = S.copyOf()
    val random = Random()
    for (i in result.size - 1 downTo 1) {
        val j = random.nextInt(i + 1)
        val temp = result[i]
        result[i] = result[j]
        result[j] = temp
    }
    return result
}

/**
 * Computes the counting process from inter-arrival times.
 * Returns the number of arrivals within 'scale' time units from each arrival.
 *
 * @param S Inter-arrival times
 * @param scale Time window size
 * @return Counts array
 */
fun trace_iat2counts(S: DoubleArray, scale: Double): IntArray {
    val n = S.size
    if (n <= 1) return intArrayOf()

    val CS = DoubleArray(n)
    CS[0] = S[0]
    for (i in 1 until n) {
        CS[i] = CS[i - 1] + S[i]
    }

    val counts = ArrayList<Int>()
    for (i in 0 until n - 1) {
        var cur = i
        while (cur + 1 < n && CS[cur + 1] - CS[i] <= scale) {
            cur++
        }
        counts.add(cur - i)
        if (cur >= n - 1) break
    }

    return counts.toIntArray()
}

/**
 * Computes binned counts from inter-arrival times.
 *
 * @param S Inter-arrival times
 * @param scale Bin width
 * @return Pair of (counts per bin, bin membership for each event)
 */
fun trace_iat2bins(S: DoubleArray, scale: Double): Pair<IntArray, IntArray> {
    val n = S.size
    if (n == 0) return Pair(intArrayOf(), intArrayOf())

    val CS = DoubleArray(n)
    CS[0] = S[0]
    for (i in 1 until n) {
        CS[i] = CS[i - 1] + S[i]
    }

    val totalTime = CS[n - 1] - CS[0]
    val numBins = ceil(totalTime / scale).toInt()
    if (numBins <= 0) return Pair(intArrayOf(), intArrayOf())

    val counts = IntArray(numBins)
    val binMembership = ArrayList<Int>()

    var cur = 0
    var last = 0

    for (binIdx in 0 until numBins) {
        val binEnd = (binIdx + 1) * scale
        while (cur + 1 < n && CS[cur + 1] <= binEnd) {
            cur++
        }
        counts[binIdx] = cur - last
        for (j in 0 until cur - last) {
            binMembership.add(binIdx + 1) // 1-based bin index
        }
        last = cur
    }

    return Pair(counts, binMembership.toIntArray())
}

/**
 * Computes the probability mass function of discrete values in a trace.
 *
 * @param X Input trace with discrete values
 * @return Pair of (PMF values, unique values)
 */
fun trace_pmf(X: DoubleArray): Pair<DoubleArray, DoubleArray> {
    val sorted = X.sorted()
    val unique = sorted.distinct()
    val n = X.size.toDouble()

    val pmf = DoubleArray(unique.size)
    for ((idx, value) in unique.withIndex()) {
        pmf[idx] = X.count { it == value } / n
    }

    return Pair(pmf, unique.toDoubleArray())
}

/**
 * Summary statistics for a trace.
 */
data class TraceSummary(
    val mean: Double,
    val scv: Double,
    val mad: Double,
    val skewness: Double,
    val kurtosis: Double,
    val quartiles: DoubleArray,
    val percentile95: Double,
    val min: Double,
    val max: Double,
    val iqr: Double,
    val acf: DoubleArray,
    val idc: Double,
    val length: Int
) {
    override fun equals(other: Any?): Boolean {
        if (this === other) return true
        if (javaClass != other?.javaClass) return false
        other as TraceSummary
        if (mean != other.mean) return false
        if (scv != other.scv) return false
        return true
    }

    override fun hashCode(): Int {
        var result = mean.hashCode()
        result = 31 * result + scv.hashCode()
        return result
    }
}

/**
 * Computes comprehensive summary statistics for a trace.
 *
 * @param m Input trace
 * @return TraceSummary containing all computed statistics
 */
fun trace_summary(m: DoubleArray): TraceSummary {
    val mean = trace_mean(m)
    val scv = trace_scv(m)
    val skewness = Skewness().evaluate(m)
    val kurtosis = Kurtosis().evaluate(m)

    val percentile = Percentile()
    percentile.data = m
    val q25 = percentile.evaluate(25.0)
    val q50 = percentile.evaluate(50.0)
    val q75 = percentile.evaluate(75.0)
    val p95 = percentile.evaluate(95.0)
    val iqr = q75 - q25

    // Median absolute deviation
    val median = q50
    val absDeviations = m.map { FastMath.abs(it - median) }.toDoubleArray()
    val mad = Percentile().evaluate(absDeviations, 50.0)

    val acf = trace_acf(m, (1..10).toList().toIntArray())
    val idc = trace_idc(m)

    return TraceSummary(
        mean = mean,
        scv = scv,
        mad = mad,
        skewness = skewness,
        kurtosis = kurtosis,
        quartiles = doubleArrayOf(q25, q50, q75),
        percentile95 = p95,
        min = m.minOrNull() ?: Double.NaN,
        max = m.maxOrNull() ?: Double.NaN,
        iqr = iqr,
        acf = acf,
        idc = idc,
        length = m.size
    )
}

/**
 * Computes the mean of multiple traces (multi-trace mean).
 *
 * @param traces List of traces
 * @return Array of means for each trace
 */
fun mtrace_mean(traces: List<DoubleArray>): DoubleArray {
    return traces.map { trace_mean(it) }.toDoubleArray()
}
