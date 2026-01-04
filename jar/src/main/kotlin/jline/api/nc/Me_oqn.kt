package jline.api.nc

import jline.io.line_warning
import jline.util.matrix.Matrix
import kotlin.math.abs
import kotlin.math.max

/**
 * Maximum Entropy algorithm for Open Queueing Networks.
 *
 * Implements the ME algorithm from Kouvatsos (1994) "Entropy Maximisation
 * and Queueing Network Models", Section 3.2.
 *
 * @since LINE 3.0
 */

/**
 * Options for the ME OQN algorithm.
 *
 * @property tol Convergence tolerance (default: 1e-6)
 * @property maxIter Maximum number of iterations (default: 1000)
 * @property verbose Print iteration information (default: false)
 */
data class MeOqnOptions(
    val tol: Double = 1e-6,
    val maxIter: Int = 1000,
    val verbose: Boolean = false
)

/**
 * Result of the ME OQN algorithm.
 *
 * @property L Mean queue lengths [M x R matrix]
 * @property W Mean waiting times [M x R matrix]
 * @property Ca Arrival SCV at each queue [M x R matrix]
 * @property Cd Departure SCV at each queue [M x R matrix]
 * @property lambda Total arrival rates [M x R matrix]
 * @property rho Utilizations [M x R matrix]
 */
data class MeOqnResult(
    val L: Matrix,
    val W: Matrix,
    val Ca: Matrix,
    val Cd: Matrix,
    val lambda: Matrix,
    val rho: Matrix
)

/**
 * Maximum Entropy algorithm for Open Queueing Networks.
 *
 * Implements the ME algorithm from Kouvatsos (1994).
 *
 * @param M Number of queues (stations)
 * @param R Number of job classes
 * @param lambda0 External arrival rates [M x R matrix], lambda0(i,r) = λoi,r
 * @param Ca0 External arrival SCV [M x R matrix], Ca0(i,r) = Caoi,r
 * @param mu Service rates [M x R matrix], mu(i,r) = μi,r
 * @param Cs Service SCV [M x R matrix], Cs(i,r) = Csi,r
 * @param P Routing probability matrix [M x M x R], P(j,i,r) = pji,r
 *          (probability class r goes from queue j to queue i)
 * @param options Algorithm options
 * @return Result containing performance metrics
 *
 * Reference: Kouvatsos (1994), Equations 3.6 and 3.7
 */
fun me_oqn(
    M: Int,
    R: Int,
    lambda0: Matrix,
    Ca0: Matrix,
    mu: Matrix,
    Cs: Matrix,
    P: Array<Array<Matrix>>,
    options: MeOqnOptions = MeOqnOptions()
): MeOqnResult {
    // Step 1: Feedback correction
    // If pii,r > 0, apply feedback elimination transformation
    val P_eff = Array(M) { j -> Array(M) { i -> Matrix.zeros(R, 1) } }
    val mu_eff = mu.copy()
    val Cs_eff = Cs.copy()

    for (i in 0 until M) {
        for (j in 0 until M) {
            for (r in 0 until R) {
                P_eff[j][i].set(r, 0, P[j][i].get(r, 0))
            }
        }
    }

    for (i in 0 until M) {
        for (r in 0 until R) {
            val pii = P[i][i].get(r, 0)
            if (pii > 0) {
                // Feedback correction: adjust service rate and scv
                mu_eff.set(i, r, mu.get(i, r) * (1 - pii))
                // Adjusted service scv accounting for feedback
                Cs_eff.set(i, r, Cs.get(i, r) / (1 - pii) + pii / (1 - pii))
                // Remove self-loop from routing
                P_eff[i][i].set(r, 0, 0.0)
                // Renormalize routing probabilities
                var rowSum = 0.0
                for (dest in 0 until M) {
                    rowSum += P_eff[i][dest].get(r, 0)
                }
                if (rowSum > 0) {
                    for (dest in 0 until M) {
                        val newProb = P_eff[i][dest].get(r, 0) / rowSum * (1 - pii)
                        P_eff[i][dest].set(r, 0, newProb)
                    }
                }
            }
        }
    }

    // Step 2: Initialize arrival scv
    val Ca = Matrix.ones(M, R)

    // Step 3: Solve job flow balance equations using ORIGINAL P (including feedback)
    // λi,r = λoi,r + Σj λj,r * pji,r
    val lambda = Matrix.zeros(M, R)

    for (r in 0 until R) {
        // Build the system (I - P') * lambda = lambda0
        val Pr = Matrix.zeros(M, M)
        for (j in 0 until M) {
            for (i in 0 until M) {
                Pr.set(j, i, P[j][i].get(r, 0))
            }
        }
        val PrT = Pr.transpose()
        val A = Matrix.eye(M).sub(PrT)
        val lambda0_r = Matrix.extractColumn(lambda0, r, null)
        val lambda_r = Matrix.zeros(M, 1)
        Matrix.solve(A, lambda0_r, lambda_r)
        for (i in 0 until M) {
            lambda.set(i, r, lambda_r.get(i, 0))
        }
    }

    // Compute utilizations using ORIGINAL mu (not mu_eff)
    val rho = Matrix.zeros(M, R)
    for (i in 0 until M) {
        for (r in 0 until R) {
            if (mu.get(i, r) > 0) {
                rho.set(i, r, lambda.get(i, r) / mu.get(i, r))
            }
        }
    }

    // Check stability
    val rhoTotal = Matrix.zeros(M, 1)
    for (i in 0 until M) {
        var sum = 0.0
        for (r in 0 until R) {
            sum += rho.get(i, r)
        }
        rhoTotal.set(i, 0, sum)
    }

    for (i in 0 until M) {
        if (rhoTotal.get(i, 0) >= 1.0) {
            line_warning("me_oqn", "Network is unstable (utilization >= 1 at queue %d)", i)
        }
    }

    // Initialize outputs
    val L = Matrix.zeros(M, R)
    val Cd = Matrix.ones(M, R)

    // Step 4-5: Iterative computation of Ca, Cd, and L
    for (iter in 1..options.maxIter) {
        val Ca_old = Ca.copy()

        // Step 4: Apply GE-type formulae for mean queue length
        for (i in 0 until M) {
            var rho_i = 0.0
            for (r in 0 until R) {
                rho_i += rho.get(i, r)
            }

            if (rho_i > 0 && rho_i < 1) {
                for (r in 0 until R) {
                    if (rho.get(i, r) > 0) {
                        // Proportion of class r traffic
                        val prop_r = rho.get(i, r) / rho_i
                        // GE/GE/1 mean queue length formula (ME approximation)
                        val Ca_agg = Ca.get(i, r)
                        val Cs_agg = Cs_eff.get(i, r)
                        val L_total = rho_i + (rho_i * rho_i * (Ca_agg + Cs_agg)) / (2.0 * (1.0 - rho_i))
                        L.set(i, r, prop_r * L_total)
                    }
                }
            }
        }

        // Step 5a: Compute departure scv using equation (3.6)
        for (j in 0 until M) {
            var rho_j = 0.0
            for (r in 0 until R) {
                rho_j += rho.get(j, r)
            }
            for (r in 0 until R) {
                if (lambda.get(j, r) > 0) {
                    val Cd_jr = 2.0 * L.get(j, r) * (1.0 - rho_j) + Ca.get(j, r) * (1.0 - 2.0 * rho_j)
                    Cd.set(j, r, max(0.0, Cd_jr))
                }
            }
        }

        // Step 5b: Compute arrival scv using equation (3.7)
        for (i in 0 until M) {
            for (r in 0 until R) {
                if (lambda.get(i, r) > 0) {
                    var sumInv = 0.0

                    // Contribution from other queues
                    for (j in 0 until M) {
                        val pji = P_eff[j][i].get(r, 0)
                        if (pji > 0 && lambda.get(j, r) > 0) {
                            // Thinning formula for departure scv after splitting
                            val Cdji = 1.0 + pji * (Cd.get(j, r) - 1.0)
                            val weight = (lambda.get(j, r) * pji) / lambda.get(i, r)
                            sumInv += weight / (Cdji + 1.0)
                        }
                    }

                    // Contribution from external arrivals
                    if (lambda0.get(i, r) > 0) {
                        val weight0 = lambda0.get(i, r) / lambda.get(i, r)
                        sumInv += weight0 / (Ca0.get(i, r) + 1.0)
                    }

                    // Compute new arrival scv
                    if (sumInv > 0) {
                        val Ca_ir = -1.0 + 1.0 / sumInv
                        Ca.set(i, r, max(0.0, Ca_ir))
                    }
                }
            }
        }

        // Check convergence
        var delta = 0.0
        for (i in 0 until M) {
            for (r in 0 until R) {
                delta = max(delta, abs(Ca.get(i, r) - Ca_old.get(i, r)))
            }
        }

        if (options.verbose) {
            println("Iteration $iter: max delta = $delta")
        }

        if (delta < options.tol) {
            if (options.verbose) {
                println("Converged after $iter iterations")
            }
            break
        }

        if (iter == options.maxIter && delta >= options.tol) {
            line_warning("me_oqn", "Did not converge within %d iterations (delta=%f)",
                         options.maxIter, delta)
        }
    }

    // Step 6: Compute final statistics
    val W = Matrix.zeros(M, R)
    for (i in 0 until M) {
        for (r in 0 until R) {
            if (lambda.get(i, r) > 0) {
                W.set(i, r, L.get(i, r) / lambda.get(i, r))
            }
        }
    }

    return MeOqnResult(L, W, Ca, Cd, lambda, rho)
}

/**
 * ME OQN algorithms documentation marker for Dokka.
 */
@Suppress("unused")
class MeOqn {
    companion object {
        // Class documentation marker for Dokka
    }
}
