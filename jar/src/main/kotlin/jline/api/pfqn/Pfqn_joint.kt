/**
 * Joint Queue-Length Probability for Product-Form Networks
 *
 * Computes the joint queue-length probability for closed product-form queueing networks
 * using normalizing constants. Supports both total queue-lengths and per-class queue-lengths.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn

import jline.api.pfqn.nc.pfqn_ca
import jline.lib.perm.Permanent
import jline.util.Maths
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath
import kotlin.math.exp
import kotlin.math.ln

/**
 * Compute the joint queue-length probability for vector n
 *
 * @param n Queue population vector (M x 1 for total queue lengths, or M x R for per-class queue lengths)
 * @param L Service demands (M x R matrix, M queues, R classes)
 * @param N Total number of jobs for each class (1 x R vector)
 * @param Z Think times for each class (1 x R vector, default: zeros)
 * @param lGN Logarithm of normalizing constant at N (optional, computed if not provided)
 * @return Joint queue-length probability
 *
 * Examples:
 * - Total queue-lengths: n is M x 1 vector
 * - Per-class queue-lengths: n is M x R matrix
 */
fun pfqn_joint(n: Matrix, L: Matrix, N: Matrix, Z: Matrix? = null, lGN: Double? = null): Double {
    val M = L.numRows
    val R = L.numCols

    // Handle default Z parameter
    val zThink = Z ?: Matrix.zeros(1, R)

    // Compute normalizing constant if not provided
    val lGn = lGN ?: run {
        val result = pfqn_ca(L, N, zThink)
        result.lG
    }

    return if (n.numCols == 1) {
        // Joint probability of total queue lengths
        computeTotalQueueLengthProb(n, L, N, zThink, lGn)
    } else if (n.numCols == R) {
        // Joint probability of per-class queue-lengths
        computePerClassQueueLengthProb(n, L, N, zThink, lGn, M, R)
    } else {
        throw IllegalArgumentException("Invalid argument to pfqn_joint: n must have 1 or $R columns")
    }
}

/**
 * Compute joint probability for total queue lengths (n is M x 1)
 */
private fun computeTotalQueueLengthProb(
    n: Matrix,
    L: Matrix,
    N: Matrix,
    Z: Matrix,
    lGn: Double
): Double {
    val hasThinkTime = Z.elementSum() > 0.0

    return if (hasThinkTime) {
        val n0 = N.elementSum() - n.elementSum()
        val Fjoint = fper(Matrix.concatRows(L, Z, null), N, Matrix.concatRows(n, Matrix.singleton(n0), null))
        val logFjoint = ln(Fjoint)
        val logFactorial = Maths.factln(n0.toInt())
        exp(logFjoint - lGn - logFactorial)
    } else {
        val Fjoint = fper(L, N, n)
        val logFjoint = ln(Fjoint)
        exp(logFjoint - lGn)
    }
}

/**
 * Compute joint probability for per-class queue lengths (n is M x R)
 */
private fun computePerClassQueueLengthProb(
    n: Matrix,
    L: Matrix,
    N: Matrix,
    Z: Matrix,
    lGn: Double,
    M: Int,
    R: Int
): Double {
    // n0 = N - sum(n,1)  (sum across rows for each class)
    val n0 = Matrix(1, R)
    for (r in 0 until R) {
        var colSum = 0.0
        for (i in 0 until M) {
            colSum += n.get(i, r)
        }
        n0.set(0, r, N.get(0, r) - colSum)
    }

    var Fjoint = 0.0

    // Handle think time contribution
    val hasThinkTime = Z.elementSum() > 0.0
    if (hasThinkTime) {
        for (r in 0 until R) {
            if (n0.get(0, r) > 0) {
                Fjoint += n0.get(0, r) * ln(Z.get(0, r))
                Fjoint -= Maths.factln(n0.get(0, r).toInt())
            }
        }
    }

    // Compute multinomial contributions for each queue
    for (i in 0 until M) {
        val nRow = Matrix(1, R)
        for (r in 0 until R) {
            nRow.set(0, r, n.get(i, r))
        }
        Fjoint += Maths.multinomialln(nRow)

        // Add service demand contributions
        for (r in 0 until R) {
            if (n.get(i, r) > 0 && L.get(i, r) > 0) {
                Fjoint += n.get(i, r) * ln(L.get(i, r))
            }
        }
    }

    return exp(Fjoint - lGn)
}

/**
 * Helper function F_per: computes permanent-based probability term
 *
 * @param L Service demands matrix
 * @param N Population vector
 * @param m State vector (queue populations)
 * @return Probability term based on permanent
 */
private fun fper(L: Matrix, N: Matrix, m: Matrix): Double {
    val M = L.numRows
    val R = L.numCols

    // Build Ak: replicate each column of L by N(r) times
    var Ak: Matrix? = null
    for (r in 0 until R) {
        val nRep = N.get(0, r).toInt()
        val col = Matrix.extractColumn(L, r, null)
        val replicatedCols = col.repmat(1, nRep)
        Ak = if (Ak == null) {
            replicatedCols
        } else {
            Matrix.concatColumns(Ak, replicatedCols, null)
        }
    }

    // Build A: for each queue i with m(i) > 0, replicate Ak(i,:) m(i) times
    var A: Matrix? = null
    for (i in 0 until M) {
        val mi = m.get(i, 0).toInt()
        if (mi > 0) {
            val rowToReplicate = Matrix.extractRows(Ak!!, i, i + 1, null)
            val replicatedRows = rowToReplicate.repmat(mi, 1)
            A = if (A == null) {
                replicatedRows
            } else {
                Matrix.concatRows(A, replicatedRows, null)
            }
        }
    }

    // Handle case where A is null (all m(i) = 0)
    if (A == null || A.numRows == 0) {
        return 1.0
    }

    // Handle non-square matrix (invalid state - probability is 0)
    if (A.numRows != A.numCols) {
        return 0.0
    }

    // Compute permanent of A
    val permanent = Permanent(A, true)
    val permValue = permanent.value

    // Divide by product of factorials
    var logProdFactorial = 0.0
    for (r in 0 until R) {
        logProdFactorial += Maths.factln(N.get(0, r).toInt())
    }

    return permValue / exp(logProdFactorial)
}
