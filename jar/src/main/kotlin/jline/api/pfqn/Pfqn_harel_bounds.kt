/**
 * Harel et al. throughput bounds for closed queueing networks
 *
 * Implements the throughput bounds from:
 * A. Harel, S. Namn, J. Sturm (1999) "Simple bounds for closed queueing networks",
 * Queueing Systems 31:125-135.
 *
 * These bounds are sharper than the classical Balanced Job Bounds (BJB) of Zahorjan et al.
 *
 * Lower Bound (LB):
 *   LB = N / (A1 + (N-1) * [AN/A1]^(1/(N-1)))
 *
 * Upper Bounds (UB):
 *   UB(n) = N / (A1 + ((N-1)/(n-1)) * (n/TH(n) - A1))
 *
 * where:
 *   Ai = sum(rho_j^i) for j=1..k
 *   TH(n) = G(n-1) / G(n) is the throughput with n customers
 *   G(n) is the normalizing constant (complete symmetric function)
 *
 * @since LINE 3.0
 */
package jline.api.pfqn

import jline.io.Ret
import jline.util.matrix.Matrix
import kotlin.math.pow

/**
 * Compute Harel et al. throughput bounds for a closed single-class queueing network.
 *
 * @param rho Loading vector (k x 1 matrix), where rho_i = x_i / mu_i
 *            (visit ratio divided by service rate for queue i)
 * @param N Number of customers in the closed network
 * @param Z Think time (must be zero; throws exception if nonzero)
 * @param maxUB Maximum n for upper bound computation (default: min(N, 7))
 * @return Ret.pfqnHarelBounds containing LB, UB(n) for n=2..maxUB, and TH(n) values
 * @throws IllegalArgumentException if Z is nonzero or if parameters are invalid
 */
@JvmOverloads
fun pfqn_harel_bounds(
    rho: Matrix,
    N: Int,
    Z: Double = 0.0,
    maxUB: Int? = null
): Ret.pfqnHarelBounds {
    // Validate think time is zero
    if (Z != 0.0) {
        throw IllegalArgumentException(
            "Harel bounds are only valid for networks with zero think time. " +
            "The provided think time Z=$Z is nonzero."
        )
    }

    // Validate population
    if (N < 1) {
        throw IllegalArgumentException("Population N must be at least 1, got N=$N")
    }

    // Get number of queues
    val k = rho.numRows
    if (k < 1) {
        throw IllegalArgumentException("Loading vector rho must have at least 1 element")
    }

    // Validate all rho values are positive
    for (i in 0 until k) {
        if (rho.get(i, 0) <= 0) {
            throw IllegalArgumentException("All loading factors must be positive, rho[$i]=${rho.get(i, 0)}")
        }
    }

    // Determine maximum n for upper bounds (limited by G(n) polynomial availability)
    val effectiveMaxUB = maxUB ?: minOf(N, 7)
    if (effectiveMaxUB > 7) {
        throw IllegalArgumentException(
            "Upper bounds UB(n) can only be computed for n <= 7 due to G(n) polynomial limitations. " +
            "Requested maxUB=$effectiveMaxUB"
        )
    }

    // Compute A_i values: A_i = sum(rho_j^i) for i=1..N
    val A = computeAValues(rho, N)

    // Compute Lower Bound (LB)
    // LB = N / (A1 + (N-1) * [AN/A1]^(1/(N-1)))
    val LB = computeLowerBound(A, N)

    // Compute TH(n) for n=1..effectiveMaxUB using G(n) polynomials
    // TH(n) = G(n-1) / G(n)
    val TH = computeThroughputValues(A, effectiveMaxUB)

    // Compute Upper Bounds UB(n) for n=2..effectiveMaxUB
    // UB(n) = N / (A1 + ((N-1)/(n-1)) * (n/TH(n) - A1))
    val UB = computeUpperBounds(A, N, TH, effectiveMaxUB)

    return Ret.pfqnHarelBounds(LB, UB, TH, N, k, effectiveMaxUB)
}

/**
 * Compute A_i values: A_i = sum(rho_j^i) for i=1..maxPower
 *
 * @param rho Loading vector
 * @param maxPower Maximum power to compute (typically N)
 * @return Array of A values, A[i] = sum(rho_j^i) for i=1..maxPower (A[0] is unused)
 */
private fun computeAValues(rho: Matrix, maxPower: Int): DoubleArray {
    val k = rho.numRows
    val A = DoubleArray(maxPower + 1)

    for (i in 1..maxPower) {
        var sum = 0.0
        for (j in 0 until k) {
            sum += rho.get(j, 0).pow(i.toDouble())
        }
        A[i] = sum
    }

    return A
}

/**
 * Compute the Lower Bound (LB) from Theorem 2.8
 *
 * LB = N / (A1 + (N-1) * [AN/A1]^(1/(N-1)))
 *
 * @param A Array of A_i values
 * @param N Population size
 * @return Lower bound on throughput
 */
private fun computeLowerBound(A: DoubleArray, N: Int): Double {
    if (N == 1) {
        // Special case: N=1, LB = 1/A1
        return 1.0 / A[1]
    }

    val A1 = A[1]
    val AN = A[N]
    val exponent = 1.0 / (N - 1)

    return N.toDouble() / (A1 + (N - 1) * (AN / A1).pow(exponent))
}

/**
 * Compute G(n) using the polynomial formulas from the paper's appendix.
 *
 * The polynomials are given for n!G(n), so we divide by n! to get G(n).
 *
 * @param A Array of A_i values (A[i] = sum(rho_j^i))
 * @param n Population for which to compute G(n)
 * @return G(n) normalizing constant
 */
private fun computeG(A: DoubleArray, n: Int): Double {
    return when (n) {
        0 -> 1.0

        1 -> A[1]

        2 -> {
            // 2*G(2) = A1^2 + A2
            (A[1] * A[1] + A[2]) / 2.0
        }

        3 -> {
            // 6*G(3) = A1^3 + 3*A2*A1 + 2*A3
            val A1 = A[1]
            val A2 = A[2]
            val A3 = A[3]
            (A1 * A1 * A1 + 3 * A2 * A1 + 2 * A3) / 6.0
        }

        4 -> {
            // 24*G(4) = A1^4 + 6*A2*A1^2 + 8*A3*A1 + 3*A2^2 + 6*A4
            val A1 = A[1]
            val A2 = A[2]
            val A3 = A[3]
            val A4 = A[4]
            (A1.pow(4.0) + 6 * A2 * A1 * A1 + 8 * A3 * A1 + 3 * A2 * A2 + 6 * A4) / 24.0
        }

        5 -> {
            // 120*G(5) = A1^5 + 10*A2*A1^3 + 20*A3*A1^2 + 15*A1*A2^2 + 30*A4*A1 + 20*A2*A3 + 24*A5
            val A1 = A[1]
            val A2 = A[2]
            val A3 = A[3]
            val A4 = A[4]
            val A5 = A[5]
            (A1.pow(5.0) + 10 * A2 * A1.pow(3.0) + 20 * A3 * A1 * A1 +
             15 * A1 * A2 * A2 + 30 * A4 * A1 + 20 * A2 * A3 + 24 * A5) / 120.0
        }

        6 -> {
            // 720*G(6) = A1^6 + 15*A2*A1^4 + 40*A3*A1^3 + 45*A1^2*A2^2 + 90*A4*A1^2
            //          + 120*A1*A2*A3 + 144*A5*A1 + 15*A2^3 + 90*A2*A4 + 40*A3^2 + 120*A6
            val A1 = A[1]
            val A2 = A[2]
            val A3 = A[3]
            val A4 = A[4]
            val A5 = A[5]
            val A6 = A[6]
            (A1.pow(6.0) + 15 * A2 * A1.pow(4.0) + 40 * A3 * A1.pow(3.0) +
             45 * A1 * A1 * A2 * A2 + 90 * A4 * A1 * A1 + 120 * A1 * A2 * A3 +
             144 * A5 * A1 + 15 * A2.pow(3.0) + 90 * A2 * A4 + 40 * A3 * A3 + 120 * A6) / 720.0
        }

        7 -> {
            // 5040*G(7) = A1^7 + 21*A2*A1^5 + 70*A3*A1^4 + 105*A1^3*A2^2 + 210*A4*A1^3
            //           + 504*A5*A1^2 + 105*A1*A2^3 + 280*A1*A3^2 + 210*A2^2*A3
            //           + 504*A2*A5 + 420*A3*A4 + 840*A6*A1 + 630*A1*A2*A4
            //           + 420*A1^2*A2*A3 + 720*A7
            val A1 = A[1]
            val A2 = A[2]
            val A3 = A[3]
            val A4 = A[4]
            val A5 = A[5]
            val A6 = A[6]
            val A7 = A[7]
            (A1.pow(7.0) + 21 * A2 * A1.pow(5.0) + 70 * A3 * A1.pow(4.0) +
             105 * A1.pow(3.0) * A2 * A2 + 210 * A4 * A1.pow(3.0) + 504 * A5 * A1 * A1 +
             105 * A1 * A2.pow(3.0) + 280 * A1 * A3 * A3 + 210 * A2 * A2 * A3 +
             504 * A2 * A5 + 420 * A3 * A4 + 840 * A6 * A1 + 630 * A1 * A2 * A4 +
             420 * A1 * A1 * A2 * A3 + 720 * A7) / 5040.0
        }

        else -> {
            throw IllegalArgumentException("G(n) polynomial not available for n=$n (max supported: 7)")
        }
    }
}

/**
 * Compute throughput values TH(n) = G(n-1)/G(n) for n=1..maxN
 *
 * @param A Array of A_i values
 * @param maxN Maximum population for throughput computation
 * @return Array of TH values, TH[n] = throughput with n customers (TH[0] is unused)
 */
private fun computeThroughputValues(A: DoubleArray, maxN: Int): DoubleArray {
    val TH = DoubleArray(maxN + 1)

    for (n in 1..maxN) {
        val Gn = computeG(A, n)
        val Gn_1 = computeG(A, n - 1)
        TH[n] = Gn_1 / Gn
    }

    return TH
}

/**
 * Compute Upper Bounds UB(n) from Theorem 3.2
 *
 * UB(n) = N / (A1 + ((N-1)/(n-1)) * (n/TH(n) - A1))
 *
 * Note: UB(N) = TH(N) (exact throughput)
 *
 * @param A Array of A_i values
 * @param N Population size
 * @param TH Array of throughput values
 * @param maxN Maximum n for upper bounds
 * @return Array of UB values, UB[n] for n=2..maxN (UB[0] and UB[1] are unused)
 */
private fun computeUpperBounds(A: DoubleArray, N: Int, TH: DoubleArray, maxN: Int): DoubleArray {
    val UB = DoubleArray(maxN + 1)
    val A1 = A[1]

    for (n in 2..maxN) {
        if (n > TH.size - 1) break

        val nOverTH = n.toDouble() / TH[n]
        val denominator = A1 + ((N - 1).toDouble() / (n - 1)) * (nOverTH - A1)
        UB[n] = N.toDouble() / denominator
    }

    return UB
}

/**
 * Compute only the lower bound LB (convenience function).
 *
 * @param rho Loading vector (k x 1 matrix)
 * @param N Number of customers
 * @param Z Think time (must be zero)
 * @return Lower bound on throughput
 */
@JvmOverloads
fun pfqn_harel_lb(rho: Matrix, N: Int, Z: Double = 0.0): Double {
    if (Z != 0.0) {
        throw IllegalArgumentException(
            "Harel LB is only valid for networks with zero think time. " +
            "The provided think time Z=$Z is nonzero."
        )
    }

    if (N < 1) {
        throw IllegalArgumentException("Population N must be at least 1, got N=$N")
    }

    val A = computeAValues(rho, N)
    return computeLowerBound(A, N)
}

/**
 * Compute a specific upper bound UB(n) (convenience function).
 *
 * @param rho Loading vector (k x 1 matrix)
 * @param N Number of customers
 * @param n The n value for UB(n), must be 2 <= n <= min(N, 7)
 * @param Z Think time (must be zero)
 * @return Upper bound UB(n) on throughput
 */
@JvmOverloads
fun pfqn_harel_ub(rho: Matrix, N: Int, n: Int, Z: Double = 0.0): Double {
    if (Z != 0.0) {
        throw IllegalArgumentException(
            "Harel UB is only valid for networks with zero think time. " +
            "The provided think time Z=$Z is nonzero."
        )
    }

    if (N < 1) {
        throw IllegalArgumentException("Population N must be at least 1, got N=$N")
    }

    if (n < 2) {
        throw IllegalArgumentException("n must be at least 2 for upper bound UB(n), got n=$n")
    }

    if (n > N) {
        throw IllegalArgumentException("n cannot exceed N for UB(n), got n=$n, N=$N")
    }

    if (n > 7) {
        throw IllegalArgumentException(
            "UB(n) can only be computed for n <= 7 due to G(n) polynomial limitations, got n=$n"
        )
    }

    // Compute A values up to n
    val A = computeAValues(rho, n)

    // Compute TH(n)
    val Gn = computeG(A, n)
    val Gn_1 = computeG(A, n - 1)
    val THn = Gn_1 / Gn

    // Compute UB(n)
    val A1 = A[1]
    val nOverTH = n.toDouble() / THn
    val denominator = A1 + ((N - 1).toDouble() / (n - 1)) * (nOverTH - A1)

    return N.toDouble() / denominator
}

/**
 * Documentation marker class for Dokka
 */
@Suppress("unused")
class PfqnHarelBoundsAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}
