package jline.solvers.mam.handlers

import jline.api.mc.ctmc_makeinfgen
import jline.api.mc.ctmc_solve
import jline.util.matrix.Matrix
import kotlin.math.abs
import kotlin.random.Random

data class INAPResult(
    val x: Matrix,
    val pi: List<Matrix>,
    val Q: List<Matrix>,
    val iter: Int
)

fun solver_mam_inap(
    rcat: RCATModel,
    tol: Double = 1e-6,
    maxiter: Int = 1000,
    seed: Int = 0,
    method: String = "inap"
): INAPResult {
    val R = rcat.R
    val AP = rcat.AP
    val N = rcat.N

    val numActions = rcat.actionMap.size
    val numProcesses = N.size

    if (numProcesses == 0) {
        return INAPResult(
            x = Matrix(0, 1),
            pi = emptyList(),
            Q = emptyList(),
            iter = 0
        )
    }

    val Aa = Array<Matrix?>(numActions) { a -> R[a][0] }
    val Pb = Array<Matrix?>(numActions) { a -> R[a][1] }
    val L = Array<Matrix?>(numProcesses) { k -> R[numActions][k] }

    val ACT = IntArray(numActions) { a -> AP.get(a, 0).toInt() }
    val PSV = IntArray(numActions) { a -> AP.get(a, 1).toInt() }

    val random = if (seed != 0) Random(seed) else Random
    val x = Matrix(numActions, 1)
    for (a in 0 until numActions) {
        x.set(a, 0, random.nextDouble())
    }

    var (pi, Q) = computeEquilibrium(x, Aa, Pb, L, ACT, PSV, numProcesses, numActions, N)

    var iter = 1
    while (iter <= maxiter) {
        val piprev = pi.map { it.copy() }

        for (a in 0 until numActions) {
            val k = ACT[a]
            val activeMatrix = Aa[a] ?: continue
            val piK = pi[k]

            if (method == "inapplus") {
                // inapplus: LAMBDA(i,j) = Aa(i,j) * pi(i)
                // x(a) = sum(LAMBDA) for non-zero entries
                var lambdaSum = 0.0
                for (i in 0 until N[k]) {
                    for (j in 0 until N[k]) {
                        val aij = activeMatrix.get(i, j)
                        if (aij > 0) {
                            val piI = piK.get(0, i)
                            lambdaSum += aij * piI
                        }
                    }
                }
                if (lambdaSum > 0) {
                    x.set(a, 0, lambdaSum)
                }
            } else {
                // inap: LAMBDA(i,j) = Aa(i,j) * pi(i) / pi(j)
                // x(a) = mean(LAMBDA) for non-zero entries
                val lambdaVec = mutableListOf<Double>()
                for (i in 0 until N[k]) {
                    for (j in 0 until N[k]) {
                        val aij = activeMatrix.get(i, j)
                        val piJ = piK.get(0, j)
                        if (aij > 0 && piJ > 0) {
                            val piI = piK.get(0, i)
                            lambdaVec.add(aij * piI / piJ)
                        }
                    }
                }
                if (lambdaVec.isNotEmpty()) {
                    x.set(a, 0, lambdaVec.average())
                }
            }
        }

        val result = computeEquilibrium(x, Aa, Pb, L, ACT, PSV, numProcesses, numActions, N)
        pi = result.first
        Q = result.second

        var maxErr = 0.0
        for (k in 0 until numProcesses) {
            var errK = 0.0
            for (j in 0 until N[k]) {
                errK += abs(pi[k].get(0, j) - piprev[k].get(0, j))
            }
            if (errK > maxErr) maxErr = errK
        }

        if (maxErr < tol) {
            return INAPResult(x, pi, Q, iter)
        }

        iter++
    }

    return INAPResult(x, pi, Q, iter)
}

private fun computeEquilibrium(
    x: Matrix,
    Aa: Array<Matrix?>,
    Pb: Array<Matrix?>,
    L: Array<Matrix?>,
    ACT: IntArray,
    PSV: IntArray,
    numProcesses: Int,
    numActions: Int,
    N: IntArray
): Pair<List<Matrix>, List<Matrix>> {
    val Q = mutableListOf<Matrix>()
    val pi = mutableListOf<Matrix>()

    for (k in 0 until numProcesses) {
        val size = N[k]

        val Lk = L[k]
        var Qk = if (Lk != null) {
            Lk.copy()
        } else {
            Matrix(size, size)
        }

        val diagVec = Matrix(1, size)
        for (i in 0 until size) {
            var rowSum = 0.0
            for (j in 0 until size) {
                if (i != j) {
                    rowSum += Qk.get(i, j)
                }
            }
            diagVec.set(0, i, -rowSum)
        }
        for (i in 0 until size) {
            Qk.set(i, i, diagVec.get(0, i))
        }

        for (c in 0 until numActions) {
            if (PSV[c] == k) {
                val pbMatrix = Pb[c] ?: continue
                val xc = x.get(c, 0)
                for (i in 0 until size) {
                    for (j in 0 until size) {
                        if (i != j) {
                            Qk.set(i, j, Qk.get(i, j) + xc * pbMatrix.get(i, j))
                        }
                    }
                }
            } else if (ACT[c] == k) {
                val aaMatrix = Aa[c] ?: continue
                for (i in 0 until size) {
                    for (j in 0 until size) {
                        if (i != j) {
                            Qk.set(i, j, Qk.get(i, j) + aaMatrix.get(i, j))
                        }
                    }
                }
            }
        }

        Qk = ctmc_makeinfgen(Qk)
        Q.add(Qk)

        val piK = ctmc_solve(Qk)
        pi.add(piK)
    }

    return Pair(pi, Q)
}
