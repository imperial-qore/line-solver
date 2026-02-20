package jline.api.mc

import jline.util.matrix.Matrix

/**
 * Solve reducible CTMCs by converting to DTMC via randomization
 *
 * @param Q Infinitesimal generator matrix (possibly reducible)
 * @param pi0 Initial state distribution (optional)
 * @param options Solution options (optional)
 * @return Pair of (pi: steady-state distribution, scc: strongly connected components info)
 */
fun ctmc_solve_reducible(
    Q: Matrix,
    pi0: Matrix? = null,
    options: Map<String, Any> = mapOf("tol" to 1e-12)
): Pair<Matrix, List<List<Int>>> {
    // Convert CTMC to DTMC via randomization
    val (P, _) = ctmc_randomization(Q)

    // Solve the reducible DTMC
    val dtmcResult = dtmc_solve_reducible(P, pi0, options)
    val pi = dtmcResult.first
    val scc = dtmcResult.second

    return Pair(pi, scc)
}

/**
 * Alternative signature that returns additional information
 */
data class CtmcSolveReducibleResult(
    val pi: Matrix,
    val pis: List<Matrix>,
    val pi0: Matrix,
    val scc: List<List<Int>>,
    val isrec: List<Boolean>
)

fun ctmc_solve_reducible_full(
    Q: Matrix,
    pi0: Matrix? = null,
    options: Map<String, Any> = mapOf("tol" to 1e-12)
): CtmcSolveReducibleResult {
    // Convert CTMC to DTMC via randomization
    val (P, _) = ctmc_randomization(Q)

    // Solve the reducible DTMC (would need full implementation)
    val result = dtmc_solve_reducible_full(P, pi0, options)

    return CtmcSolveReducibleResult(
        pi = result.pi,
        pis = result.pis,
        pi0 = result.pi0,
        scc = result.scc,
        isrec = result.isrec
    )
}
/**
 * CTMC solve reducible algorithms
 */
@Suppress("unused")
class CtmcSolveReducibleAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}
