package jline.api.map

import jline.util.matrix.Matrix
import org.apache.commons.math3.distribution.PoissonDistribution
import org.apache.commons.math3.util.FastMath
import java.util.*
import kotlin.math.*

/**
 * MAP/M/1-PS Sojourn Time Distribution
 *
 * Computes the complementary distribution function of sojourn time in a
 * MAP/M/1 processor-sharing queue using the algorithm from:
 *
 * Masuyama, H., & Takine, T. (2003). Sojourn time distribution in a
 * MAP/M/1 processor-sharing queue. Operations Research Letters, 31(6), 406-412.
 *
 * @author LINE Development Team
 */
object MAPM1PSCdfRespT {

    /**
     * Compute complementary sojourn time CDF for MAP/M/1-PS queue
     *
     * @param C MAP C matrix (transitions without arrivals)
     * @param D MAP D matrix (transitions with arrivals)
     * @param mu Service rate (must be > 0)
     * @param x Time points for CDF evaluation
     * @param epsilon Queue length truncation parameter (default: 1e-11)
     * @param epsilonPrime Uniformization truncation parameter (default: 1e-10)
     * @return Complementary CDF values W̄(x) = Pr[W > x] at each point in x
     */
    @JvmStatic
    @JvmOverloads
    fun computeCdf(
        C: Matrix,
        D: Matrix,
        mu: Double,
        x: DoubleArray,
        epsilon: Double = 1e-11,
        epsilonPrime: Double = 1e-10
    ): DoubleArray {

        val M = C.getNumRows()
        require(M == C.getNumCols()) { "C must be square" }
        require(M == D.getNumRows() && M == D.getNumCols()) { "C and D must have same dimensions" }
        require(mu > 0) { "Service rate mu must be positive" }

        // Compute stationary distribution π of underlying Markov chain
        val pi = computeStationaryDistribution(C, D)

        // Compute mean arrival rate λ
        val e = Matrix.ones(M, 1)
        val lambda = pi.mult(D).mult(e).get(0, 0)

        // Check stability condition
        val rho = lambda / mu
        require(rho < 1.0) { "System is unstable (rho = $rho >= 1)" }

        // Compute R matrix
        val R = computeRMatrix(C, D, mu)

        // Compute π₀ = π(I - R)
        val I = Matrix.eye(M)
        val pi0 = pi.mult(I.sub(1.0, R))

        // Determine N(epsilon) - queue length truncation
        val Nepsilon = determineNEpsilon(pi0, R, D, lambda, epsilon, e)

        // Compute theta (uniformization parameter)
        val theta = (0 until M).map { abs(C.get(it, it)) }.maxOrNull() ?: 0.0

        // Allocate result array
        val result = DoubleArray(x.size)

        // Main computation loop for each x value
        for (idx in x.indices) {
            val xVal = x[idx]

            // Determine K_max for uniformization
            val thetaPlusMu = theta + mu
            val meanVal = thetaPlusMu * xVal

            val L: Int
            val Kmax: Int

            if (meanVal > 0) {
                L = max(0, floor(meanVal - 10 * sqrt(meanVal)).toInt())
                Kmax = ceil(meanVal + 10 * sqrt(meanVal)).toInt()
            } else {
                L = 0
                Kmax = 0
            }

            // Compute h_{n,k} recursively
            val h = computeHRecursive(C, D, mu, Nepsilon, Kmax, theta)

            // Compute W̄(x) = (1/λ) Σ π₀R^n D Σ Poisson * h_{n,k}
            var WBar = 0.0

            for (n in 0..Nepsilon) {
                // Compute weight = π₀ * R^n * D
                val weight = if (n == 0) {
                    pi0.mult(D)
                } else {
                    pi0.mult(Matrix.pow(R, n)).mult(D)
                }

                // Compute sum over k
                var sumK = Matrix.zeros(M, 1)
                for (k in L..Kmax) {
                    val poissonTerm = poissonPmf(k, meanVal)
                    sumK = sumK.add(1.0, h[n][k].scale(poissonTerm))
                }

                WBar += (1.0 / lambda) * weight.mult(sumK).get(0, 0)
            }

            result[idx] = WBar
        }

        return result
    }

    /**
     * Compute stationary distribution of MAP
     */
    private fun computeStationaryDistribution(C: Matrix, D: Matrix): Matrix {
        val M = C.getNumRows()
        val Q = C.add(1.0, D)

        // Solve π * Q = 0 with π * e = 1
        // Convert to Q^T * π^T = 0, replace last equation with sum = 1
        val A = Q.transpose()

        // Replace last row with normalization constraint
        for (j in 0 until M) {
            A.set(M - 1, j, 1.0)
        }

        val b = Matrix.zeros(M, 1)
        b.set(M - 1, 0, 1.0)

        // Solve system
        val piT = Matrix.zeros(M, 1)
        Matrix.solve(A, b, piT)
        return piT.transpose()
    }

    /**
     * Compute R matrix (minimal nonnegative solution of D + R(C - μI) + μR² = O)
     */
    private fun computeRMatrix(C: Matrix, D: Matrix, mu: Double): Matrix {
        val M = C.getNumRows()

        if (M == 1) {
            // Scalar case: use quadratic formula
            // μR² + (C - μ)R + D = 0
            val a = mu
            val b = C.get(0, 0) - mu
            val c = D.get(0, 0)

            val discriminant = b * b - 4 * a * c
            require(discriminant >= 0) { "No real solution for R matrix" }

            val R1 = (-b - sqrt(discriminant)) / (2 * a)
            val R2 = (-b + sqrt(discriminant)) / (2 * a)

            // Choose minimal nonnegative solution
            val Rval = when {
                R1 >= 0 && R1 < 1 -> R1
                R2 >= 0 && R2 < 1 -> R2
                else -> throw IllegalStateException("No valid R in [0,1)")
            }

            val result = Matrix(1, 1)
            result.set(0, 0, Rval)
            return result
        } else {
            // Matrix case: Newton iteration
            val I = Matrix.eye(M)
            val CminusMuI = C.sub(mu, I)

            // Initial guess: R = -D(C - μI)^{-1}
            var R = D.scale(-1.0).mult(CminusMuI.inv())

            val maxIter = 1000
            val tol = 1e-12

            for (iter in 0 until maxIter) {
                val Rold = R.copy()

                // F(R) = D + R(C - μI) + μR²
                val F = D.add(1.0, R.mult(CminusMuI)).add(mu, R.mult(R))

                // F'(R) = (C - μI) + 2μR
                val Fprime = CminusMuI.add(2 * mu, R)

                // Newton step: R_new = R_old - F(R) * inv(F'(R))
                R = R.sub(1.0, F.mult(Fprime.inv()))

                // Check convergence
                if (R.sub(1.0, Rold).elementMaxAbs() < tol) {
                    break
                }
            }

            // Ensure non-negativity
            for (i in 0 until M) {
                for (j in 0 until M) {
                    if (R.get(i, j) < 0) R.set(i, j, 0.0)
                }
            }

            return R
        }
    }

    /**
     * Determine N(epsilon) - minimum n such that truncated sum > 1 - epsilon
     */
    private fun determineNEpsilon(
        pi0: Matrix,
        R: Matrix,
        D: Matrix,
        lambda: Double,
        epsilon: Double,
        e: Matrix
    ): Int {
        val M = R.getNumRows()
        var cumsumProb = 0.0
        var Rpower = Matrix.eye(M)

        for (n in 0..1000) {
            cumsumProb += (1.0 / lambda) * pi0.mult(Rpower).mult(D).mult(e).get(0, 0)
            if (cumsumProb > 1 - epsilon) {
                return n
            }
            Rpower = Rpower.mult(R)
        }

        return 100  // Default fallback
    }

    /**
     * Compute h_{n,k} recursively using Algorithm from Theorem 1
     *
     * h_{n,0} = e for all n
     * h_{n,k+1} = 1/(θ+μ) * [nμ/(n+1) h_{n-1,k} + (θI + C)h_{n,k} + D h_{n+1,k}]
     */
    private fun computeHRecursive(
        C: Matrix,
        D: Matrix,
        mu: Double,
        N: Int,
        K: Int,
        theta: Double
    ): Array<Array<Matrix>> {

        val M = C.getNumRows()
        val I = Matrix.eye(M)
        val e = Matrix.ones(M, 1)

        val thetaPlusMu = theta + mu
        val thetaIPlusC = I.scale(theta).add(1.0, C)

        // Initialize h array: h[n][k] where n ∈ [0,N], k ∈ [0,K]
        val h = Array(N + 1) { Array(K + 1) { Matrix.zeros(M, 1) } }

        // Base case: h_{n,0} = e for all n
        for (n in 0..N) {
            h[n][0] = e.copy()
        }

        // Recursive computation: fill column by column (increasing k)
        for (k in 0 until K) {
            for (n in 0..N) {
                var term1 = Matrix.zeros(M, 1)
                val term2 = thetaIPlusC.mult(h[n][k])
                var term3 = Matrix.zeros(M, 1)

                // First term: nμ/(n+1) * h_{n-1,k}
                if (n > 0) {
                    term1 = h[n - 1][k].scale((n * mu) / (n + 1))
                }

                // Third term: D * h_{n+1,k}
                if (n < N) {
                    term3 = D.mult(h[n + 1][k])
                }

                h[n][k + 1] = term1.add(1.0, term2).add(1.0, term3).scale(1.0 / thetaPlusMu)
            }
        }

        return h
    }

    /**
     * Compute Poisson PMF avoiding overflow
     */
    private fun poissonPmf(k: Int, lambda: Double): Double {
        if (lambda == 0.0) {
            return if (k == 0) 1.0 else 0.0
        }

        // Use log-space computation to avoid overflow
        // P(k; λ) = exp(k*log(λ) - λ - log(k!))
        var logProb = k * ln(lambda) - lambda

        // Subtract log(k!)
        for (i in 1..k) {
            logProb -= ln(i.toDouble())
        }

        return exp(logProb)
    }
}
