package jline.lib.kpctoolbox.aph

import jline.api.mam.map_normalize
import jline.api.mam.map_renewal
import jline.api.mam.map_scale
import jline.api.mam.map_exponential
import jline.lib.kpctoolbox.mc.dtmc_solve
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import org.apache.commons.math3.util.FastMath
import java.util.Random

/**
 * Acyclic Phase-Type (APH) distribution functions.
 *
 * Ported from MATLAB: matlab/lib/kpctoolbox/aph/
 */

/**
 * Convolution patterns for APH simplification.
 */
enum class ConvolutionPattern {
    SEQUENCE,   // Pattern 1: Sequential structure
    PARALLEL,   // Pattern 2: Parallel structure
    BRANCH      // Pattern 3: Branch structure
}

/**
 * Simplifies/combines two APH distributions using a specified pattern.
 *
 * @param a1 Initial probability vector of first distribution
 * @param T1 Rate matrix of first distribution
 * @param a2 Initial probability vector of second distribution
 * @param T2 Rate matrix of second distribution
 * @param p1 Branch probability for first distribution (for BRANCH pattern)
 * @param p2 Branch probability for second distribution (for BRANCH pattern)
 * @param pattern Convolution pattern to use
 * @return Pair of (combined alpha, combined T)
 */
fun aph_simplify(
    a1: DoubleArray,
    T1: Matrix,
    a2: DoubleArray,
    T2: Matrix,
    p1: Double = 1.0,
    p2: Double = 1.0,
    pattern: ConvolutionPattern
): Pair<DoubleArray, Matrix> {
    val order1 = a1.size
    val order2 = a2.size

    when (pattern) {
        ConvolutionPattern.SEQUENCE -> {
            // Sequential structure: first complete first distribution, then second
            val n = order1 + order2

            // Compute sum of a1 (probability of starting in first distribution)
            var a1Sum = 0.0
            for (v in a1) a1Sum += v

            // alpha = [a1, (1-sum(a1))*a2]
            val alpha = DoubleArray(n)
            for (i in 0 until order1) alpha[i] = a1[i]
            for (i in 0 until order2) alpha[order1 + i] = (1.0 - a1Sum) * a2[i]

            // T = [T1, (-T1*e1)*a2; 0, T2]
            val T = Matrix(n, n)

            // Copy T1
            for (i in 0 until order1) {
                for (j in 0 until order1) {
                    T.set(i, j, T1.get(i, j))
                }
            }

            // Compute -T1*e1 (exit rates from T1)
            val exitRates = DoubleArray(order1)
            for (i in 0 until order1) {
                var sum = 0.0
                for (j in 0 until order1) {
                    sum += T1.get(i, j)
                }
                exitRates[i] = -sum
            }

            // (-T1*e1)*a2 block
            for (i in 0 until order1) {
                for (j in 0 until order2) {
                    T.set(i, order1 + j, exitRates[i] * a2[j])
                }
            }

            // Copy T2
            for (i in 0 until order2) {
                for (j in 0 until order2) {
                    T.set(order1 + i, order1 + j, T2.get(i, j))
                }
            }

            return Pair(alpha, T)
        }

        ConvolutionPattern.PARALLEL -> {
            // Parallel structure: minimum of two distributions
            val n = order1 * order2 + order1 + order2

            // Compute sums
            var a1Sum = 0.0
            var a2Sum = 0.0
            for (v in a1) a1Sum += v
            for (v in a2) a2Sum += v

            // alpha = [kron(a1,a2), (1-sum(a2))*a1, (1-sum(a1))*a2]
            val alpha = DoubleArray(n)
            var idx = 0

            // Kronecker product part
            for (i in 0 until order1) {
                for (j in 0 until order2) {
                    alpha[idx++] = a1[i] * a2[j]
                }
            }

            // (1-sum(a2))*a1
            for (i in 0 until order1) {
                alpha[idx++] = (1.0 - a2Sum) * a1[i]
            }

            // (1-sum(a1))*a2
            for (i in 0 until order2) {
                alpha[idx++] = (1.0 - a1Sum) * a2[i]
            }

            // Build T matrix for parallel structure
            val T = Matrix(n, n)

            // Compute exit rates
            val exit1 = DoubleArray(order1)
            val exit2 = DoubleArray(order2)
            for (i in 0 until order1) {
                var sum = 0.0
                for (j in 0 until order1) sum += T1.get(i, j)
                exit1[i] = -sum
            }
            for (i in 0 until order2) {
                var sum = 0.0
                for (j in 0 until order2) sum += T2.get(i, j)
                exit2[i] = -sum
            }

            // First block: kron(T1, I) + kron(I, T2)
            for (i in 0 until order1) {
                for (j in 0 until order2) {
                    val rowIdx = i * order2 + j
                    for (k in 0 until order1) {
                        for (l in 0 until order2) {
                            val colIdx = k * order2 + l
                            var value = 0.0
                            if (j == l) value += T1.get(i, k)
                            if (i == k) value += T2.get(j, l)
                            T.set(rowIdx, colIdx, value)
                        }
                    }
                }
            }

            // Transition to second block: kron(I, -T2*e2)
            for (i in 0 until order1) {
                for (j in 0 until order2) {
                    val rowIdx = i * order2 + j
                    val colIdx = order1 * order2 + i
                    T.set(rowIdx, colIdx, exit2[j])
                }
            }

            // Transition to third block: kron(-T1*e1, I)
            for (i in 0 until order1) {
                for (j in 0 until order2) {
                    val rowIdx = i * order2 + j
                    val colIdx = order1 * order2 + order1 + j
                    T.set(rowIdx, colIdx, exit1[i])
                }
            }

            // Second block: T1
            for (i in 0 until order1) {
                for (j in 0 until order1) {
                    T.set(order1 * order2 + i, order1 * order2 + j, T1.get(i, j))
                }
            }

            // Third block: T2
            for (i in 0 until order2) {
                for (j in 0 until order2) {
                    T.set(order1 * order2 + order1 + i, order1 * order2 + order1 + j, T2.get(i, j))
                }
            }

            return Pair(alpha, T)
        }

        ConvolutionPattern.BRANCH -> {
            // Branch structure: probabilistic choice between distributions
            val n = order1 + order2

            // alpha = [p1*a1, p2*a2]
            val alpha = DoubleArray(n)
            for (i in 0 until order1) alpha[i] = p1 * a1[i]
            for (i in 0 until order2) alpha[order1 + i] = p2 * a2[i]

            // T = [T1, 0; 0, T2]
            val T = Matrix(n, n)
            for (i in 0 until order1) {
                for (j in 0 until order1) {
                    T.set(i, j, T1.get(i, j))
                }
            }
            for (i in 0 until order2) {
                for (j in 0 until order2) {
                    T.set(order1 + i, order1 + j, T2.get(i, j))
                }
            }

            return Pair(alpha, T)
        }
    }
}

/**
 * Performs convolution on parallel structure with any number of elements.
 *
 * @param distributions List of pairs (alpha, T) for each distribution
 * @return Combined (alpha, T)
 */
fun aph_convpara(distributions: List<Pair<DoubleArray, Matrix>>): Pair<DoubleArray, Matrix> {
    if (distributions.isEmpty()) {
        throw IllegalArgumentException("Need at least one distribution")
    }
    if (distributions.size == 1) {
        return distributions[0]
    }

    var (alpha, T) = aph_simplify(
        distributions[0].first, distributions[0].second,
        distributions[1].first, distributions[1].second,
        pattern = ConvolutionPattern.PARALLEL
    )

    for (i in 2 until distributions.size) {
        val result = aph_simplify(
            alpha, T,
            distributions[i].first, distributions[i].second,
            pattern = ConvolutionPattern.SEQUENCE
        )
        alpha = result.first
        T = result.second
    }

    return Pair(alpha, T)
}

/**
 * Performs convolution on sequential structure with any number of elements.
 *
 * @param distributions List of pairs (alpha, T) for each distribution
 * @return Combined (alpha, T)
 */
fun aph_convseq(distributions: List<Pair<DoubleArray, Matrix>>): Pair<DoubleArray, Matrix> {
    if (distributions.isEmpty()) {
        throw IllegalArgumentException("Need at least one distribution")
    }
    if (distributions.size == 1) {
        return distributions[0]
    }

    var (alpha, T) = aph_simplify(
        distributions[0].first, distributions[0].second,
        distributions[1].first, distributions[1].second,
        pattern = ConvolutionPattern.SEQUENCE
    )

    for (i in 2 until distributions.size) {
        val result = aph_simplify(
            alpha, T,
            distributions[i].first, distributions[i].second,
            pattern = ConvolutionPattern.SEQUENCE
        )
        alpha = result.first
        T = result.second
    }

    return Pair(alpha, T)
}

/**
 * Generates a random APH (acyclic phase-type) distribution as a MAP.
 *
 * @param K Order of the APH distribution (default: 2)
 * @return MAP representation {D0, D1}
 */
fun aph_rand(K: Int = 2): MatrixCell {
    val random = Random()

    val D0 = Matrix(K, K)
    val D1 = Matrix(K, K)

    // Generate random upper triangular D0 (acyclic)
    for (i in 0 until K) {
        for (j in i until K) {
            D0.set(i, j, random.nextDouble())
        }
        // Clear lower triangular part
        for (j in 0 until i) {
            D0.set(i, j, 0.0)
        }
    }

    // Generate random D1
    for (i in 0 until K) {
        for (j in 0 until K) {
            D1.set(i, j, random.nextDouble())
        }
    }

    val MAP = MatrixCell(2)
    MAP[0] = D0
    MAP[1] = D1

    return map_normalize(map_renewal(MAP))
}

/**
 * Fits an APH distribution to match first three moments.
 *
 * Implementation based on: A.Bobbio, A.Horvath, M.Telek, "Matching three moments
 * with minimal acyclic phase type distributions", Stochastic Models 21:303-326, 2005.
 *
 * @param e1 First moment (mean)
 * @param e2 Second moment
 * @param e3 Third moment
 * @param nmax Maximum order to try (default: 10)
 * @return Pair of (fitted MAP, isExact flag)
 */
fun aph_fit(e1: Double, e2: Double, e3: Double, nmax: Int = 10): Pair<MatrixCell, Boolean> {
    var isExact = true

    if (e2.isInfinite() || e3.isInfinite()) {
        return Pair(map_exponential(e1), true)
    }

    val n2 = e2 / (e1 * e1)
    val n3 = e3 / (e1 * e2)

    // Find suitable order
    var n2Feas = false
    var n3UbFeas = false
    var n3LbFeas = false
    var n = 1
    var un = 0.0
    var unPrev = 0.0

    while ((!n2Feas || !n3LbFeas || !n3UbFeas) && n < nmax) {
        n++
        unPrev = un

        val pn = ((n + 1) * (n2 - 2) / (3 * n2 * (n - 1))) *
                (-2 * FastMath.sqrt(n + 1.0) / FastMath.sqrt(4.0 * (n + 1) - 3 * n * n2) - 1)
        val an = (n2 - 2) / (pn * (1 - n2) + FastMath.sqrt(pn * pn + pn * n * (n2 - 2) / (n - 1)))
        val ln = ((3 + an) * (n - 1) + 2 * an) / ((n - 1) * (1 + an * pn)) -
                (2 * an * (n + 1)) / (2 * (n - 1) + an * pn * (n * an + 2 * n - 2))
        un = (1.0 / (n * n * n2)) * (2 * (n - 2) * (n * n2 - n - 1) *
                FastMath.sqrt(1 + n * (n2 - 2) / (n - 1)) + (n + 2) * (3 * n * n2 - 2 * n - 2))

        if (n2 >= (n + 1.0) / n && n2 <= (n + 4.0) / (n + 1)) {
            n2Feas = true
            if (n3 >= ln) n3LbFeas = true
        } else if (n2 >= (n + 4.0) / (n + 1)) {
            n2Feas = true
            if (n3 >= n2 * (n + 1) / n) n3LbFeas = true
        }

        if (n2 >= (n + 1.0) / n && n2 <= n.toDouble() / (n - 1)) {
            n2Feas = true
            if (n3 <= un) n3UbFeas = true
        } else if (n2 >= n.toDouble() / (n - 1)) {
            n2Feas = true
            n3UbFeas = true
        }
    }

    var fitN2 = n2
    var fitN3 = n3

    if (!n2Feas || !n3LbFeas || !n3UbFeas || n >= nmax) {
        // Cannot match exactly, use feasible approximation
        fitN2 = (n + 1.0) / n
        fitN3 = 2 * fitN2 - 1
        isExact = false
    }

    // Fitting algorithm
    val MAP: MatrixCell

    if (fitN2 <= n.toDouble() / (n - 1) || fitN3 <= 2 * fitN2 - 1) {
        // Case 1
        val b = 2 * (4 - n * (3 * fitN2 - 4)) / (fitN2 * (4 + n - n * fitN3) +
                FastMath.sqrt(n * fitN2) * FastMath.sqrt(12 * fitN2 * fitN2 * (n + 1) +
                        16 * fitN3 * (n + 1) + fitN2 * (n * (fitN3 - 15) * (fitN3 + 1) - 8 * (fitN3 + 3))))
        val a = (b * fitN2 - 2) * (n - 1) * b / ((b - 1) * n)
        val p = (b - 1) / a
        val lambda = 1.0
        val mu = lambda * (n - 1) / a

        val alpha = DoubleArray(n)
        alpha[0] = p
        alpha[n - 1] = 1 - p

        val T = Matrix(n, n)
        for (i in 0 until n) {
            T.set(i, i, -mu)
            if (i < n - 1) T.set(i, i + 1, mu)
        }
        T.set(n - 1, n - 1, -lambda)

        // Build D0 = T and D1 = -T*e*alpha
        val D0 = T.copy()
        val D1 = Matrix(n, n)
        for (i in 0 until n) {
            var exitRate = 0.0
            for (j in 0 until n) exitRate += T.get(i, j)
            exitRate = -exitRate
            for (j in 0 until n) {
                D1.set(i, j, exitRate * alpha[j])
            }
        }

        MAP = MatrixCell(2)
        MAP[0] = D0
        MAP[1] = D1
    } else {
        // Case 2 - more complex fitting
        // Simplified implementation using approximation
        val lambda = 1.0
        val mu = lambda * n / e1

        val alpha = DoubleArray(n)
        alpha[0] = 1.0

        val T = Matrix(n, n)
        for (i in 0 until n) {
            T.set(i, i, -mu)
            if (i < n - 1) T.set(i, i + 1, mu)
        }

        val D0 = T.copy()
        val D1 = Matrix(n, n)
        for (i in 0 until n) {
            var exitRate = 0.0
            for (j in 0 until n) exitRate += T.get(i, j)
            exitRate = -exitRate
            for (j in 0 until n) {
                D1.set(i, j, exitRate * alpha[j])
            }
        }

        MAP = MatrixCell(2)
        MAP[0] = D0
        MAP[1] = D1
        isExact = false
    }

    return Pair(map_scale(map_normalize(MAP), e1), isExact)
}

/**
 * Converts a hyper-exponential PH distribution to its rate/probability form.
 *
 * @param PH Phase-type distribution as {D0, D1}
 * @return Pair of (rates, probabilities)
 */
fun ph2hyper(PH: MatrixCell): Pair<DoubleArray, DoubleArray> {
    val D0 = PH[0]
    val n = D0.numRows

    // Check if diagonal (hyper-exponential)
    for (i in 0 until n) {
        for (j in 0 until n) {
            if (i != j && FastMath.abs(D0.get(i, j)) > 1e-10) {
                throw IllegalArgumentException("The PH distribution is not hyper-exponential")
            }
        }
    }

    // Extract rates (negative diagonal of D0)
    val lambda = DoubleArray(n)
    for (i in 0 until n) {
        lambda[i] = -D0.get(i, i)
    }

    // Compute probabilities from stationary distribution of embedded DTMC
    val D1 = PH[1]

    // P = inv(-D0) * D1
    val P = Matrix(n, n)
    for (i in 0 until n) {
        for (j in 0 until n) {
            P.set(i, j, D1.get(i, j) / lambda[i])
        }
    }

    val prob = dtmc_solve(P)

    return Pair(lambda, prob)
}

/**
 * Generates random samples from a hyper-exponential distribution.
 *
 * @param rates Exponential rates
 * @param probs Selection probabilities
 * @param nSamples Number of samples
 * @return Array of samples
 */
fun hyper_rand(rates: DoubleArray, probs: DoubleArray, nSamples: Int): DoubleArray {
    val random = Random()
    val samples = DoubleArray(nSamples)

    // Compute cumulative probabilities
    val cumProbs = DoubleArray(probs.size)
    cumProbs[0] = probs[0]
    for (i in 1 until probs.size) {
        cumProbs[i] = cumProbs[i - 1] + probs[i]
    }

    for (s in 0 until nSamples) {
        // Select which exponential to use
        val u = random.nextDouble()
        var selected = 0
        for (i in cumProbs.indices) {
            if (u <= cumProbs[i]) {
                selected = i
                break
            }
        }

        // Generate exponential sample
        samples[s] = -FastMath.log(random.nextDouble()) / rates[selected]
    }

    return samples
}
