/**
 * BMAP/MAP/1 Queue Solver using M/G/1 type analysis with ETAQA.
 *
 * Solves a BMAP/MAP/1 queue where:
 * - Arrivals follow a Batch Markovian Arrival Process (BMAP)
 * - Service follows a Markovian Arrival Process (MAP) for single service
 * - Single server
 *
 * The queue is modeled as an M/G/1-type Markov chain because:
 * - BMAP arrivals can increase level by 1, 2, 3, ... (batch sizes)
 * - MAP service decreases level by exactly 1
 *
 * @since LINE 3.1.0
 */
package jline.solvers.mam.handlers

import jline.util.matrix.Matrix
import jline.lib.smc.mg1_pi_etaqa
import jline.lib.smc.mg1_g_etaqa
import jline.lib.smc.mg1_qlen_etaqa
import jline.api.mam.map_prob

/**
 * Result of BMAP/MAP/1 queue analysis.
 */
data class BMAPMAP1Result(
    /** Mean queue length E[N] */
    val meanQueueLength: Double,
    /** Server utilization rho */
    val utilization: Double,
    /** Mean response time E[R] */
    val meanResponseTime: Double,
    /** Throughput (total customer arrival rate) */
    val throughput: Double,
    /** Aggregated stationary probabilities [pi0, pi1, piStar] */
    val pi: Matrix,
    /** G matrix */
    val G: Matrix,
    /** Mean batch size */
    val meanBatchSize: Double
)

/**
 * Solves a BMAP/MAP/1 queue using M/G/1 type matrix-analytic methods.
 *
 * The BMAP is specified by matrices {D0, D1, D2, ..., DK} where:
 * - D0: transitions without arrivals
 * - Dk: transitions triggering batch size k arrivals (k >= 1)
 *
 * The MAP for service is specified by matrices (S0, S1) where:
 * - S0: transitions without service completions
 * - S1: transitions triggering service completions
 *
 * The M/G/1 type structure for BMAP/MAP/1 is:
 * ```
 *     A0 = I_ma \otimes S1           (service completion, level -1)
 *     A1 = D0 \otimes I_ms + I_ma \otimes S0  (phase changes, level 0)
 *     A_{k+1} = D_k \otimes I_ms     (batch arrival size k, level +k)
 * ```
 *
 * @param D BMAP matrices as list [D0, D1, D2, ..., DK]
 * @param S0 MAP service matrix for transitions without service completions
 * @param S1 MAP service matrix for service completions
 * @param nMoments Number of queue length moments to compute (default: 3)
 * @return BMAPMAP1Result with performance metrics
 */
fun solver_mam_bmap_map_1(D: List<Matrix>, S0: Matrix, S1: Matrix, nMoments: Int = 3): BMAPMAP1Result {
    require(D.isNotEmpty()) { "BMAP must have at least D0 matrix" }
    require(D.size >= 2) { "BMAP must have at least D0 and D1 matrices" }

    val K = D.size - 1  // Maximum batch size
    val ma = D[0].numRows  // Number of BMAP phases
    val ms = S0.numRows    // Number of MAP service phases
    val m = ma * ms        // Combined phases per level

    // Validate BMAP matrices are square and same size
    for (i in D.indices) {
        require(D[i].numRows == ma && D[i].numCols == ma) {
            "All BMAP matrices must be ${ma}x${ma}"
        }
    }
    require(S0.numRows == ms && S0.numCols == ms) { "S0 must be ${ms}x${ms}" }
    require(S1.numRows == ms && S1.numCols == ms) { "S1 must be ${ms}x${ms}" }

    // Compute arrival rate from BMAP
    val D0 = D[0]
    var D1Total = Matrix.zeros(ma, ma)
    for (k in 1..K) {
        D1Total = D1Total.add(D[k])
    }

    val piBmap = map_prob(D0, D1Total)
    val eMa = Matrix.ones(ma, 1)

    // Total customer arrival rate: sum_k (k * pi * Dk * e)
    var lambdaTotal = 0.0
    var batchRate = 0.0
    for (k in 1..K) {
        val rateK = piBmap.mult(D[k]).mult(eMa)[0, 0]
        lambdaTotal += k * rateK
        batchRate += rateK
    }

    // Mean batch size
    val meanBatchSize = if (batchRate > 0) lambdaTotal / batchRate else 0.0

    // Compute service rate from MAP
    val piMap = map_prob(S0, S1)
    val eMs = Matrix.ones(ms, 1)
    val mu = piMap.mult(S1).mult(eMs)[0, 0]

    // Utilization
    val rho = lambdaTotal / mu

    if (rho >= 1.0) {
        System.err.println("Warning: System is unstable (rho = $rho >= 1). Results may be invalid.")
    }

    // Construct M/G/1-type matrices using Kronecker products
    val I_ma = Matrix.eye(ma)
    val I_ms = Matrix.eye(ms)

    // A = [A0, A1, A2, ..., A_{K+1}] for levels >= 1
    val A = Matrix.zeros(m, m * (K + 2))

    // A0 = I_ma \otimes S1 (service completion)
    val A0 = I_ma.kron(S1)
    for (i in 0 until m) {
        for (j in 0 until m) {
            A[i, j] = A0[i, j]
        }
    }

    // A1 = D0 \otimes I_ms + I_ma \otimes S0 (phase changes only)
    val A1 = D0.kron(I_ms).add(I_ma.kron(S0))
    for (i in 0 until m) {
        for (j in 0 until m) {
            A[i, m + j] = A1[i, j]
        }
    }

    // A_{k+1} = D_k \otimes I_ms for k >= 1 (batch arrivals)
    for (k in 1..K) {
        val Ak = D[k].kron(I_ms)
        for (i in 0 until m) {
            for (j in 0 until m) {
                A[i, (k + 1) * m + j] = Ak[i, j]
            }
        }
    }

    // B = [B0, B1, B2, ..., BK] for level 0 (empty queue)
    val B = Matrix.zeros(m, m * (K + 1))

    // B0: At level 0, no service (add S1 back as self-loop)
    val B0 = D0.kron(I_ms).add(I_ma.kron(S0.add(S1)))
    for (i in 0 until m) {
        for (j in 0 until m) {
            B[i, j] = B0[i, j]
        }
    }

    // B_k = D_k \otimes I_ms for k >= 1 (batch arrivals from empty queue)
    for (k in 1..K) {
        val Bk = D[k].kron(I_ms)
        for (i in 0 until m) {
            for (j in 0 until m) {
                B[i, k * m + j] = Bk[i, j]
            }
        }
    }

    // Compute G matrix using ETAQA
    val G = mg1_g_etaqa(A)

    // Compute stationary probabilities using ETAQA
    val pi = mg1_pi_etaqa(B, A, G)

    // Compute queue length moments
    val qlenMoments = DoubleArray(nMoments)
    for (n in 1..nMoments) {
        qlenMoments[n - 1] = mg1_qlen_etaqa(B, A, pi, n)
    }

    // Performance metrics
    val meanQueueLength = qlenMoments[0]
    val meanResponseTime = if (lambdaTotal > 0) meanQueueLength / lambdaTotal else 0.0

    return BMAPMAP1Result(
        meanQueueLength = meanQueueLength,
        utilization = rho,
        meanResponseTime = meanResponseTime,
        throughput = lambdaTotal,
        pi = pi,
        G = G,
        meanBatchSize = meanBatchSize
    )
}
