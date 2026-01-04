package jline.solvers.ctmc.analyzers

import jline.lang.NetworkStruct
import jline.solvers.SolverOptions
import jline.solvers.ctmc.handlers.Solver_ctmc.Companion.solver_ctmc
import jline.util.matrix.Matrix

/**
 * Result class for CTMC reward computation via value iteration.
 *
 * @property valueFunction Map from reward name to value function matrix [Tmax+1 x nStates]
 * @property time Time vector scaled by uniformization rate
 * @property rewardNames List of reward names in order
 * @property stateSpace The state space matrix [nStates x nDims]
 * @property steadyState Map from reward name to steady-state expected reward
 * @property runtime Computation time in seconds
 */
data class RewardResult(
    val valueFunction: Map<String, Matrix>,
    val time: DoubleArray,
    val rewardNames: List<String>,
    val stateSpace: Matrix,
    val steadyState: Map<String, Double>,
    val runtime: Double
)

/**
 * CTMC reward analyzer using value iteration with uniformization.
 *
 * Computes cumulative rewards via the Bellman equation:
 *   V^{k+1}(s) = r(s) + sum_{s'} P(s,s') * V^k(s')
 *
 * where P is the uniformized transition probability matrix.
 */
class Solver_ctmc_reward {
    companion object {
        /**
         * Compute rewards via value iteration on uniformized CTMC.
         *
         * @param sn NetworkStruct with reward definitions in sn.reward
         * @param options SolverOptions with rewardIterations setting
         * @return RewardResult containing value functions and steady-state rewards
         */
        @JvmStatic
        fun solver_ctmc_reward(sn: NetworkStruct, options: SolverOptions): RewardResult {
            val startTime = System.nanoTime()

            // Validate rewards are defined
            if (sn.reward == null || sn.reward.isEmpty()) {
                throw IllegalStateException("No rewards defined. Use model.setReward() before calling reward analysis.")
            }

            // Get number of iterations
            val Tmax = options.rewardIterations ?: 1000

            // Get generator and state space from CTMC solver
            val ctmcResult = solver_ctmc(sn, options)
            val Q = ctmcResult.getQ()
            val stateSpace = ctmcResult.getStateSpace()
            val stateSpaceAggr = ctmcResult.getStateSpaceAggr()

            val nstates = Q.getNumRows()
            val rewards = sn.reward
            val rewardNames = rewards.keys.toList()
            val nRewards = rewardNames.size

            // Build reward vectors (one per reward definition)
            val R = mutableMapOf<String, DoubleArray>()
            for (name in rewardNames) {
                val rewardFn = rewards[name]!!
                val rv = DoubleArray(nstates)
                for (s in 0 until nstates) {
                    // Pass the aggregated state row to the reward function
                    val stateRow = stateSpaceAggr.getRow(s)
                    rv[s] = rewardFn.compute(stateRow, sn)
                }
                R[name] = rv
            }

            // Uniformization: find maximum exit rate (absolute value of diagonal)
            var q = 0.0
            for (i in 0 until nstates) {
                val diagVal = Math.abs(Q.get(i, i))
                if (diagVal > q) {
                    q = diagVal
                }
            }
            if (q == 0.0) {
                q = 1.0  // Handle absorbing states edge case
            }

            // Build transition probability matrix P = Q/q + I
            val P = Matrix(nstates, nstates)
            for (i in 0 until nstates) {
                for (j in 0 until nstates) {
                    var pij = Q.get(i, j) / q
                    if (i == j) {
                        pij += 1.0
                    }
                    P.set(i, j, pij)
                }
            }

            // Value iteration
            val V = mutableMapOf<String, Matrix>()
            for (name in rewardNames) {
                val vMatrix = Matrix(Tmax + 1, nstates)
                vMatrix.zero()
                val r = R[name]!!

                // v_prev stores V^k, we compute V^{k+1}
                var vPrev = DoubleArray(nstates) { 0.0 }

                for (k in 0 until Tmax) {
                    val vNew = DoubleArray(nstates)
                    for (s in 0 until nstates) {
                        // V^{k+1}(s) = r(s) + sum_{s'} P(s,s') * V^k(s')
                        var sum = r[s]
                        for (sp in 0 until nstates) {
                            sum += P.get(s, sp) * vPrev[sp]
                        }
                        vNew[s] = sum
                        vMatrix.set(k + 1, s, sum)
                    }
                    vPrev = vNew
                }
                V[name] = vMatrix
            }

            // Time scaling: convert iteration index to continuous time
            val time = DoubleArray(Tmax + 1) { it.toDouble() / q }

            // Compute steady-state expected rewards
            // Average reward rate = mean(V(end,:)) / t(end)
            val steadyState = mutableMapOf<String, Double>()
            val tEnd = time[Tmax]
            for (name in rewardNames) {
                val vMatrix = V[name]!!
                var sum = 0.0
                for (s in 0 until nstates) {
                    sum += vMatrix.get(Tmax, s)
                }
                val avgReward = if (tEnd > 0) (sum / nstates) / tEnd else 0.0
                steadyState[name] = avgReward
            }

            val runtime = (System.nanoTime() - startTime) / 1e9

            return RewardResult(
                valueFunction = V,
                time = time,
                rewardNames = rewardNames,
                stateSpace = stateSpaceAggr,
                steadyState = steadyState,
                runtime = runtime
            )
        }
    }
}
