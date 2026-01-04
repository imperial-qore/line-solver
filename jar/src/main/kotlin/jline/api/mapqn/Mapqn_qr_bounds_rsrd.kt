/**
 * @file QR Bounds via RSRD Method
 *
 * Implements Quadratic Reduction bounds using the Repetitive-Service Random-Destination
 * (RSRD) method for MAP queueing networks. This is a direct port of qrf_rsrd.m from
 * the qrf-revised repository.
 *
 * @since LINE 3.0
 */
package jline.api.mapqn

import org.apache.commons.math3.optim.linear.LinearConstraint
import org.apache.commons.math3.optim.linear.LinearConstraintSet
import org.apache.commons.math3.optim.linear.LinearObjectiveFunction
import org.apache.commons.math3.optim.linear.Relationship
import org.apache.commons.math3.optim.linear.SimplexSolver
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType
import org.apache.commons.math3.optim.MaxIter
import kotlin.math.min

/**
 * Implementation of qrf_rsrd.m linear program
 * QR Bounds with Repetitive-Service Random-Destination (RSRD)
 */
object Mapqn_qr_bounds_rsrd {

    /**
     * Create and solve the linear program for QR Bounds RSRD model
     *
     * @param params Model parameters
     * @param objectiveQueue Queue index to optimize utilization (1-based)
     * @param sense Optimization sense: "min" or "max"
     * @return Solution containing the optimal value and variable values
     */
    @JvmOverloads
    fun solve(params: Mapqn_qr_bounds_rsrd_parameters, objectiveQueue: Int, sense: String = "min"): Mapqn_solution {
        params.validate()
        require(objectiveQueue in 1..params.M) { "Objective queue must be in range 1..M" }
        require(sense == "min" || sense == "max") { "Sense must be 'min' or 'max'" }

        val M = params.M
        val N = params.N
        val F = params.F
        val K = params.K
        val r = params.r
        val mu = params.mu
        val v = params.v
        val alpha = params.alpha

        // Compute transition rates q[i][j][ki][hi][n] for n in 0..N
        // q[i,j](k,h,n) = r[i,j] * mu[i](k,h) * alpha[i](n) for j != i
        // q[i,i](k,h,n) = (v[i](k,h) + r[i,i]*mu[i](k,h)) * alpha[i](n)
        val q = Array(M) { i ->
            Array(M) { j ->
                Array(K[i]) { ki ->
                    Array(K[i]) { hi ->
                        DoubleArray(N + 1) { n ->
                            if (n == 0) {
                                0.0
                            } else {
                                val alphaVal = if (n <= alpha[i].size) alpha[i][n - 1] else 1.0
                                if (j != i) {
                                    r.get(i, j) * mu[i].get(ki, hi) * alphaVal
                                } else {
                                    (v[i].get(ki, hi) + r.get(i, i) * mu[i].get(ki, hi)) * alphaVal
                                }
                            }
                        }
                    }
                }
            }
        }

        val varIndex = mutableMapOf<String, Int>()
        var varCount = 0

        // p2 variables: p2[j,nj,kj,i,ni,hi]
        for (j in 0 until M) {
            for (nj in 0..F[j]) {
                for (kj in 0 until K[j]) {
                    for (i in 0 until M) {
                        for (ni in 0..F[i]) {
                            for (hi in 0 until K[i]) {
                                varIndex["p2_${j}_${nj}_${kj}_${i}_${ni}_${hi}"] = varCount++
                            }
                        }
                    }
                }
            }
        }

        // U and Ueff variables
        for (i in 0 until M) {
            for (ki in 0 until K[i]) {
                for (ni in 1..F[i]) {
                    varIndex["U_${i}_${ki}_${ni}"] = varCount++
                    varIndex["Ueff_${i}_${ki}_${ni}"] = varCount++
                }
            }
        }

        val nVars = varCount

        fun p2Idx(j: Int, nj: Int, kj: Int, i: Int, ni: Int, hi: Int): Int {
            return varIndex["p2_${j}_${nj}_${kj}_${i}_${ni}_${hi}"] ?: -1
        }

        fun uIdx(i: Int, ki: Int, ni: Int): Int {
            return varIndex["U_${i}_${ki}_${ni}"] ?: -1
        }

        fun ueffIdx(i: Int, ki: Int, ni: Int): Int {
            return varIndex["Ueff_${i}_${ki}_${ni}"] ?: -1
        }

        // Upper bounds for ZERO constraints
        val ub = DoubleArray(nVars) { Double.POSITIVE_INFINITY }

        for (j in 0 until M) {
            for (nj in 0..F[j]) {
                for (kj in 0 until K[j]) {
                    for (i in 0 until M) {
                        for (ni in 0..F[i]) {
                            for (hi in 0 until K[i]) {
                                val idx = p2Idx(j, nj, kj, i, ni, hi)
                                if (idx < 0) continue
                                // ZERO1: i==j, nj==ni, h<>k
                                if (i == j && nj == ni && hi != kj) ub[idx] = 0.0
                                // ZERO2: i==j, nj<>ni
                                if (i == j && nj != ni) ub[idx] = 0.0
                                // ZERO3: i<>j, nj+ni > N
                                if (i != j && nj + ni > N) ub[idx] = 0.0
                                // ZERO6: i<>j, N-nj-ni > sum of other capacities
                                if (i != j) {
                                    var sumOtherF = 0
                                    for (y in 0 until M) {
                                        if (y != i && y != j) sumOtherF += F[y]
                                    }
                                    if (N - nj - ni > sumOtherF) ub[idx] = 0.0
                                }
                            }
                        }
                    }
                }
                // ZERO7: N-nj > sum of other capacities (for diagonal)
                for (kj in 0 until K[j]) {
                    var sumOtherF = 0
                    for (y in 0 until M) {
                        if (y != j) sumOtherF += F[y]
                    }
                    if (N - nj > sumOtherF) {
                        val idx = p2Idx(j, nj, kj, j, nj, kj)
                        if (idx >= 0) ub[idx] = 0.0
                    }
                }
            }
        }

        val constraints = mutableListOf<LinearConstraint>()

        // ZERO constraints as explicit equalities
        for (idx in 0 until nVars) {
            if (ub[idx] == 0.0) {
                val coeffs = DoubleArray(nVars)
                coeffs[idx] = 1.0
                constraints.add(LinearConstraint(coeffs, Relationship.EQ, 0.0))
            }
        }

        // ONE: Normalization - sum of diagonal p2 = 1 for each queue
        for (j in 0 until M) {
            val coeffs = DoubleArray(nVars)
            for (nj in 0..F[j]) {
                for (kj in 0 until K[j]) {
                    val idx = p2Idx(j, nj, kj, j, nj, kj)
                    if (idx >= 0) coeffs[idx] = 1.0
                }
            }
            constraints.add(LinearConstraint(coeffs, Relationship.EQ, 1.0))
        }

        // SYMMETRY: p2[j,nj,k,i,ni,h] = p2[i,ni,h,j,nj,k]
        for (j in 0 until M) {
            for (nj in 0..F[j]) {
                for (kj in 0 until K[j]) {
                    for (i in j + 1 until M) {
                        for (ni in 0..F[i]) {
                            for (hi in 0 until K[i]) {
                                val idx1 = p2Idx(j, nj, kj, i, ni, hi)
                                val idx2 = p2Idx(i, ni, hi, j, nj, kj)
                                if (idx1 >= 0 && idx2 >= 0 && idx1 != idx2 &&
                                    !(ub[idx1] == 0.0 && ub[idx2] == 0.0)) {
                                    val coeffs = DoubleArray(nVars)
                                    coeffs[idx1] = 1.0
                                    coeffs[idx2] = -1.0
                                    constraints.add(LinearConstraint(coeffs, Relationship.EQ, 0.0))
                                }
                            }
                        }
                    }
                }
            }
        }

        // MARGINALS: p2[j,nj,k,j,nj,k] = sum{ni,h} p2[j,nj,k,i,ni,h] for each i != j
        for (j in 0 until M) {
            for (kj in 0 until K[j]) {
                for (nj in 0..F[j]) {
                    for (i in 0 until M) {
                        if (i == j) continue
                        val coeffs = DoubleArray(nVars)
                        val idxDiag = p2Idx(j, nj, kj, j, nj, kj)
                        if (idxDiag >= 0) coeffs[idxDiag] = 1.0
                        // Match MATLAB: iterate full range, ZERO constraints handle invalid states
                        for (ni in 0..F[i]) {
                            for (hi in 0 until K[i]) {
                                val idx = p2Idx(j, nj, kj, i, ni, hi)
                                if (idx >= 0) coeffs[idx] -= 1.0
                            }
                        }
                        constraints.add(LinearConstraint(coeffs, Relationship.EQ, 0.0))
                    }
                }
            }
        }

        // UCLASSIC: U[i,k,ni] = p2[i,ni,k,i,ni,k]
        for (i in 0 until M) {
            for (ki in 0 until K[i]) {
                for (ni in 1..F[i]) {
                    val coeffs = DoubleArray(nVars)
                    val idxU = uIdx(i, ki, ni)
                    val idxP2 = p2Idx(i, ni, ki, i, ni, ki)
                    if (idxU >= 0) coeffs[idxU] = 1.0
                    if (idxP2 >= 0) coeffs[idxP2] = -1.0
                    constraints.add(LinearConstraint(coeffs, Relationship.EQ, 0.0))
                }
            }
        }

        // UEFFS: Ueff[i,k,ni] = p2[i,ni,k,i,ni,k] - sum{j: r[i,j]>0} r[i,j]*p2[i,ni,k,j,F[j],h]
        for (i in 0 until M) {
            for (ki in 0 until K[i]) {
                for (ni in 1..F[i]) {
                    val coeffs = DoubleArray(nVars)
                    val idxUeff = ueffIdx(i, ki, ni)
                    val idxP2Diag = p2Idx(i, ni, ki, i, ni, ki)
                    if (idxUeff >= 0) coeffs[idxUeff] = 1.0
                    if (idxP2Diag >= 0) coeffs[idxP2Diag] = -1.0
                    for (j in 0 until M) {
                        if (j != i && r.get(i, j) > 0) {
                            for (hj in 0 until K[j]) {
                                val idxBlock = p2Idx(i, ni, ki, j, F[j], hj)
                                if (idxBlock >= 0) coeffs[idxBlock] = r.get(i, j)
                            }
                        }
                    }
                    constraints.add(LinearConstraint(coeffs, Relationship.EQ, 0.0))
                }
            }
        }

        // THM2: Phase balance (Theorem 1 - THM:sdeffective)
        // sum{ni} (sum{j<>i, h<>k} q[i,j,k,h,ni]*Ueff[i,k,ni] + sum{h<>k} q[i,i,k,h,ni]*U[i,k,ni])
        // = sum{ni} (sum{j<>i, h<>k} q[i,j,h,k,ni]*Ueff[i,h,ni] + sum{h<>k} q[i,i,h,k,ni]*U[i,h,ni])
        for (i in 0 until M) {
            for (ki in 0 until K[i]) {
                val coeffs = DoubleArray(nVars)
                // LHS terms
                for (ni in 1..F[i]) {
                    // Terms with Ueff (j<>i)
                    for (j in 0 until M) {
                        if (j == i) continue
                        for (hi in 0 until K[i]) {
                            if (hi == ki) continue
                            val coef = q[i][j][ki][hi][ni]
                            val idxUeff = ueffIdx(i, ki, ni)
                            if (idxUeff >= 0) coeffs[idxUeff] += coef
                        }
                    }
                    // Terms with U (j==i, h<>k)
                    for (hi in 0 until K[i]) {
                        if (hi == ki) continue
                        val coef = q[i][i][ki][hi][ni]
                        val idxP2 = p2Idx(i, ni, ki, i, ni, ki)
                        if (idxP2 >= 0) coeffs[idxP2] += coef
                    }
                }
                // RHS terms (subtract)
                for (ni in 1..F[i]) {
                    // Terms with Ueff (j<>i)
                    for (j in 0 until M) {
                        if (j == i) continue
                        for (hi in 0 until K[i]) {
                            if (hi == ki) continue
                            val coef = q[i][j][hi][ki][ni]
                            val idxUeff = ueffIdx(i, hi, ni)
                            if (idxUeff >= 0) coeffs[idxUeff] -= coef
                        }
                    }
                    // Terms with U (j==i, h<>k)
                    for (hi in 0 until K[i]) {
                        if (hi == ki) continue
                        val coef = q[i][i][hi][ki][ni]
                        val idxP2 = p2Idx(i, ni, hi, i, ni, hi)
                        if (idxP2 >= 0) coeffs[idxP2] -= coef
                    }
                }
                constraints.add(LinearConstraint(coeffs, Relationship.EQ, 0.0))
            }
        }

        // THM1: Population constraint per (j,nj,k)
        // sum{i,ni,h} ni*p2[j,nj,k,i,ni,h] = N * p2[j,nj,k,j,nj,k]
        for (j in 0 until M) {
            for (kj in 0 until K[j]) {
                for (nj in 0..F[j]) {
                    val coeffs = DoubleArray(nVars)
                    // RHS: -N * p2[j,nj,k,j,nj,k]
                    val idxDiag = p2Idx(j, nj, kj, j, nj, kj)
                    if (idxDiag >= 0) coeffs[idxDiag] = -N.toDouble()
                    // LHS: sum
                    for (i in 0 until M) {
                        for (ni in 1..F[i]) {
                            for (hi in 0 until K[i]) {
                                val idx = p2Idx(j, nj, kj, i, ni, hi)
                                if (idx >= 0) coeffs[idx] += ni.toDouble()
                            }
                        }
                    }
                    constraints.add(LinearConstraint(coeffs, Relationship.EQ, 0.0))
                }
            }
        }

        // THM3a: Marginal balance for ni in 1:F[i]-1
        for (i in 0 until M) {
            for (ni in 1 until F[i]) {
                val coeffs = DoubleArray(nVars)
                // LHS: arrivals to queue i at level ni
                for (j in 0 until M) {
                    if (j == i) continue
                    for (kj in 0 until K[j]) {
                        for (hj in 0 until K[j]) {
                            for (ui in 0 until K[i]) {
                                for (nj in 1..F[j]) {
                                    val idx = p2Idx(j, nj, kj, i, ni, ui)
                                    if (idx >= 0) coeffs[idx] += q[j][i][kj][hj][nj]
                                }
                            }
                        }
                    }
                }
                // RHS: departures from queue i at level ni+1
                for (j in 0 until M) {
                    if (j == i) continue
                    for (ki in 0 until K[i]) {
                        for (hi in 0 until K[i]) {
                            for (uj in 0 until K[j]) {
                                for (nj in 0 until F[j]) {
                                    val idx = p2Idx(i, ni + 1, ki, j, nj, uj)
                                    if (idx >= 0 && ni + 1 <= F[i]) {
                                        coeffs[idx] -= q[i][j][ki][hi][ni + 1]
                                    }
                                }
                            }
                        }
                    }
                }
                constraints.add(LinearConstraint(coeffs, Relationship.EQ, 0.0))
            }
        }

        // THM3b: Marginal balance for ni=0, per phase
        for (i in 0 until M) {
            for (ui in 0 until K[i]) {
                val coeffs = DoubleArray(nVars)
                // LHS: arrivals to queue i at ni=0
                for (j in 0 until M) {
                    if (j == i) continue
                    for (kj in 0 until K[j]) {
                        for (hj in 0 until K[j]) {
                            for (nj in 1..F[j]) {
                                val idx = p2Idx(j, nj, kj, i, 0, ui)
                                if (idx >= 0) coeffs[idx] += q[j][i][kj][hj][nj]
                            }
                        }
                    }
                }
                // RHS: departures from queue i at ni=1
                for (j in 0 until M) {
                    if (j == i) continue
                    for (ki in 0 until K[i]) {
                        for (nj in 0 until F[j]) {
                            for (hj in 0 until K[j]) {
                                val idx = p2Idx(i, 1, ki, j, nj, hj)
                                if (idx >= 0) coeffs[idx] -= q[i][j][ki][ui][1]
                            }
                        }
                    }
                }
                constraints.add(LinearConstraint(coeffs, Relationship.EQ, 0.0))
            }
        }

        // QBAL: Queue balance constraint - tightens the bounds
        for (i in 0 until M) {
            for (ki in 0 until K[i]) {
                val coeffs = DoubleArray(nVars)

                // LHS Term 1: sum{h<>k} sum{j<>i} sum{ni} sum{u} sum{nj} q[i,j,k,h,ni]*ni*p2[i,ni,k,j,nj,u]
                for (hi in 0 until K[i]) {
                    if (hi == ki) continue
                    for (j in 0 until M) {
                        if (j == i) continue
                        for (ni in 1..F[i]) {
                            for (uj in 0 until K[j]) {
                                for (nj in 0 until F[j]) {
                                    val coef = q[i][j][ki][hi][ni] * ni
                                    val idx = p2Idx(i, ni, ki, j, nj, uj)
                                    if (idx >= 0) coeffs[idx] += coef
                                }
                            }
                        }
                    }
                }

                // LHS Term 2: sum{h<>k} sum{ni} q[i,i,k,h,ni]*ni*p2[i,ni,k,i,ni,k]
                for (hi in 0 until K[i]) {
                    if (hi == ki) continue
                    for (ni in 1..F[i]) {
                        val coef = q[i][i][ki][hi][ni] * ni
                        val idx = p2Idx(i, ni, ki, i, ni, ki)
                        if (idx >= 0) coeffs[idx] += coef
                    }
                }

                // LHS Term 3: sum{j<>i} sum{h} sum{ni} sum{u} sum{nj<=min(F[j]-1,N-ni)} q[i,j,h,k,ni]*p2[i,ni,h,j,nj,u]
                for (j in 0 until M) {
                    if (j == i) continue
                    for (hi in 0 until K[i]) {
                        for (ni in 1..F[i]) {
                            for (uj in 0 until K[j]) {
                                val maxNj = min(F[j] - 1, N - ni)
                                for (nj in 0..maxNj) {
                                    val coef = q[i][j][hi][ki][ni]
                                    val idx = p2Idx(i, ni, hi, j, nj, uj)
                                    if (idx >= 0) coeffs[idx] += coef
                                }
                            }
                        }
                    }
                }

                // RHS Term 1: sum{j<>i} sum{h} sum{ni<=F[i]-1} sum{u} sum{nj>=1} q[j,i,h,u,nj]*p2[i,ni,k,j,nj,h]
                for (j in 0 until M) {
                    if (j == i) continue
                    for (hj in 0 until K[j]) {
                        for (ni in 0 until F[i]) {
                            for (uj in 0 until K[j]) {
                                for (nj in 1..F[j]) {
                                    val coef = q[j][i][hj][uj][nj]
                                    val idx = p2Idx(i, ni, ki, j, nj, hj)
                                    if (idx >= 0) coeffs[idx] -= coef
                                }
                            }
                        }
                    }
                }

                // RHS Term 2: sum{h<>k} sum{ni} q[i,i,h,k,ni]*ni*p2[i,ni,h,i,ni,h]
                for (hi in 0 until K[i]) {
                    if (hi == ki) continue
                    for (ni in 1..F[i]) {
                        val coef = q[i][i][hi][ki][ni] * ni
                        val idx = p2Idx(i, ni, hi, i, ni, hi)
                        if (idx >= 0) coeffs[idx] -= coef
                    }
                }

                // RHS Term 3: sum{h<>k} sum{j<>i} sum{ni} sum{u} sum{nj} q[i,j,h,k,ni]*ni*p2[i,ni,h,j,nj,u]
                for (hi in 0 until K[i]) {
                    if (hi == ki) continue
                    for (j in 0 until M) {
                        if (j == i) continue
                        for (ni in 1..F[i]) {
                            for (uj in 0 until K[j]) {
                                for (nj in 0 until F[j]) {
                                    val coef = q[i][j][hi][ki][ni] * ni
                                    val idx = p2Idx(i, ni, hi, j, nj, uj)
                                    if (idx >= 0) coeffs[idx] -= coef
                                }
                            }
                        }
                    }
                }

                constraints.add(LinearConstraint(coeffs, Relationship.EQ, 0.0))
            }
        }

        // THM4: Bound constraint
        // sum{t,h,nj,nt} nt*p2[j,nj,k,t,nt,h] >= N*sum{h,nj,ni} p2[j,nj,k,i,ni,h]
        for (j in 0 until M) {
            for (kj in 0 until K[j]) {
                for (i in 0 until M) {
                    val coeffs = DoubleArray(nVars)
                    // LHS: sum{t,h,nj,nt} nt*p2[j,nj,k,t,nt,h]
                    for (t in 0 until M) {
                        for (ht in 0 until K[t]) {
                            for (nj in 0..F[j]) {
                                for (nt in 1..F[t]) {
                                    val idx = p2Idx(j, nj, kj, t, nt, ht)
                                    if (idx >= 0) coeffs[idx] += nt.toDouble()
                                }
                            }
                        }
                    }
                    // RHS (negated): -N*sum{h,nj,ni} p2[j,nj,k,i,ni,h]
                    for (hi in 0 until K[i]) {
                        for (nj in 0..F[j]) {
                            for (ni in 1..F[i]) {
                                val idx = p2Idx(j, nj, kj, i, ni, hi)
                                if (idx >= 0) coeffs[idx] -= N.toDouble()
                            }
                        }
                    }
                    constraints.add(LinearConstraint(coeffs, Relationship.GEQ, 0.0))
                }
            }
        }

        // Build objective function
        val objectiveCoeffs = DoubleArray(nVars)
        val targetQueue = objectiveQueue - 1
        for (ki in 0 until K[targetQueue]) {
            for (ni in 1..F[targetQueue]) {
                val idx = p2Idx(targetQueue, ni, ki, targetQueue, ni, ki)
                if (idx >= 0) objectiveCoeffs[idx] = 1.0
            }
        }

        val goalType = if (sense == "max") GoalType.MAXIMIZE else GoalType.MINIMIZE
        val objectiveFunction = LinearObjectiveFunction(objectiveCoeffs, 0.0)

        val solver = SimplexSolver()
        val constraintSet = LinearConstraintSet(constraints)

        try {
            val solution = solver.optimize(
                objectiveFunction,
                constraintSet,
                goalType,
                MaxIter(100000)
            )

            val variables = mutableMapOf<String, Double>()

            for (i in 0 until M) {
                var totalU = 0.0
                var totalUeff = 0.0
                for (ki in 0 until K[i]) {
                    for (ni in 1..F[i]) {
                        val idxU = uIdx(i, ki, ni)
                        val idxUeff = ueffIdx(i, ki, ni)
                        if (idxU >= 0) totalU += solution.point[idxU]
                        if (idxUeff >= 0) totalUeff += solution.point[idxUeff]
                    }
                }
                variables["U_${i + 1}"] = totalU
                variables["Ueff_${i + 1}"] = totalUeff
                variables["pb_${i + 1}"] = totalU - totalUeff
            }

            return Mapqn_solution(
                objectiveValue = solution.value,
                variables = variables
            )
        } catch (e: Exception) {
            return Mapqn_solution(
                objectiveValue = Double.NaN,
                variables = emptyMap()
            )
        }
    }
}
