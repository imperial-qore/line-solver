/**
 * @file QR Bounds via BAS Method
 *
 * Implements Quadratic Reduction bounds using the Blocking-After-Service (BAS)
 * method for MAP queueing networks. This is a direct port of qrf_bas.m from
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
 * Implementation of qrf_bas.m linear program
 * QR Bounds with Blocking-After-Service (BAS)
 */
object Mapqn_qr_bounds_bas {

    /**
     * Create and solve the linear program for QR Bounds BAS model
     *
     * @param params Model parameters
     * @param objectiveQueue Queue index to optimize utilization (1-based)
     * @param sense Optimization sense: "min" or "max"
     * @return Solution containing the optimal value and variable values
     */
    @JvmOverloads
    fun solve(params: Mapqn_qr_bounds_bas_parameters, objectiveQueue: Int, sense: String = "min"): Mapqn_solution {
        params.validate()
        require(objectiveQueue in 1..params.M) { "Objective queue must be in range 1..M" }
        require(sense == "min" || sense == "max") { "Sense must be 'min' or 'max'" }

        val M = params.M
        val N = params.N
        val f = params.f - 1 // Convert to 0-based
        val F = params.F
        val K = params.K
        val mu = params.mu
        val v = params.v
        val r = params.r
        val MR = params.MR
        val BB = params.BB
        val MM = params.MM
        val ZZ = params.ZZ
        val ZM = params.ZM
        val MM1 = params.MM1

        // Compute transition rates q(i,j,k,h)
        val q = params.q()

        // Build variable indexing
        val varIndex = mutableMapOf<String, Int>()
        var varCount = 0

        // p2 variables: p2(j, nj, kj, i, ni, hi, m)
        for (j in 0 until M) {
            for (nj in 0..N) {
                for (kj in 0 until K[j]) {
                    for (i in 0 until M) {
                        for (ni in 0..N) {
                            for (hi in 0 until K[i]) {
                                for (m in 0 until MR) {
                                    varIndex["p2_${j}_${nj}_${kj}_${i}_${ni}_${hi}_${m}"] = varCount++
                                }
                            }
                        }
                    }
                }
            }
        }

        // e variables: e(i, ki)
        for (i in 0 until M) {
            for (ki in 0 until K[i]) {
                varIndex["e_${i}_${ki}"] = varCount++
            }
        }

        val nVars = varCount

        fun p2Idx(j: Int, nj: Int, kj: Int, i: Int, ni: Int, hi: Int, m: Int): Int {
            return varIndex["p2_${j}_${nj}_${kj}_${i}_${ni}_${hi}_${m}"] ?: -1
        }

        fun eIdx(i: Int, ki: Int): Int {
            return varIndex["e_${i}_${ki}"] ?: -1
        }

        // Initialize bounds
        val ub = DoubleArray(nVars) { Double.POSITIVE_INFINITY }

        // ZERO constraints - fix infeasible states
        for (j in 0 until M) {
            for (nj in 0..N) {
                for (kj in 0 until K[j]) {
                    for (i in 0 until M) {
                        for (ni in 0..N) {
                            for (hi in 0 until K[i]) {
                                for (m in 0 until MR) {
                                    val idx = p2Idx(j, nj, kj, i, ni, hi, m)
                                    if (idx < 0) continue

                                    // ZERO1: i==j, nj==ni, h<>k
                                    if (i == j && nj == ni && hi != kj) {
                                        ub[idx] = 0.0
                                    }

                                    // ZERO2: i==j, nj<>ni
                                    if (i == j && nj != ni) {
                                        ub[idx] = 0.0
                                    }

                                    // ZERO3: i<>j, nj+ni > N
                                    if (i != j && nj + ni > N) {
                                        ub[idx] = 0.0
                                    }

                                    // ZERO6: nj > F(j)
                                    if (nj > F[j]) {
                                        ub[idx] = 0.0
                                    }

                                    // ZERO5: BB(m,j)==1 and nj==0
                                    if (m >= 1 && BB.get(m, j).toInt() == 1 && nj == 0) {
                                        ub[idx] = 0.0
                                    }

                                    // ZERO7: BB(m,j)==1 and i<>j and i<>f and ni+nj+F(f)>N
                                    if (m >= 1 && BB.get(m, j).toInt() == 1 && i != j && i != f && ni + nj + F[f] > N) {
                                        ub[idx] = 0.0
                                    }

                                    // ZERO8: finite queue not at capacity in blocking config
                                    if (j == f && nj >= 1 && nj <= F[f] - 1 && m >= 1) {
                                        ub[idx] = 0.0
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        // ZERO4: For m>=1 and j<>f, p2(j,nj,k,f,nf,h,m)=0 when nf < F(f)
        for (j in 0 until M) {
            if (j == f) continue
            for (nj in 0..N) {
                for (kj in 0 until K[j]) {
                    for (m in 1 until MR) {
                        for (nf in 0 until F[f]) {
                            for (hf in 0 until K[f]) {
                                val idx = p2Idx(j, nj, kj, f, nf, hf, m)
                                if (idx >= 0) ub[idx] = 0.0
                            }
                        }
                    }
                }
            }
        }

        val constraints = mutableListOf<LinearConstraint>()

        // Add upper bound constraints for zero states
        for (idx in 0 until nVars) {
            if (ub[idx] == 0.0) {
                val coeffs = DoubleArray(nVars)
                coeffs[idx] = 1.0
                constraints.add(LinearConstraint(coeffs, Relationship.EQ, 0.0))
            }
        }

        // ONE: Normalization
        for (j in 0 until M) {
            val coeffs = DoubleArray(nVars)
            for (nj in 0..N) {
                for (kj in 0 until K[j]) {
                    for (m in 0 until MR) {
                        val idx = p2Idx(j, nj, kj, j, nj, kj, m)
                        if (idx >= 0) coeffs[idx] = 1.0
                    }
                }
            }
            constraints.add(LinearConstraint(coeffs, Relationship.EQ, 1.0))
        }

        // SYMMETRY
        for (j in 0 until M) {
            for (nj in 0..min(N, F[j])) {
                for (kj in 0 until K[j]) {
                    for (i in j + 1 until M) {
                        for (ni in 0..min(N, F[i])) {
                            if (i != j && nj + ni > N) continue
                            for (hi in 0 until K[i]) {
                                for (m in 0 until MR) {
                                    val idx1 = p2Idx(j, nj, kj, i, ni, hi, m)
                                    val idx2 = p2Idx(i, ni, hi, j, nj, kj, m)
                                    if (idx1 < 0 || idx2 < 0) continue
                                    if (ub[idx1] == 0.0 && ub[idx2] == 0.0) continue
                                    if (idx1 != idx2) {
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
        }

        // MARGINALS
        for (j in 0 until M) {
            for (kj in 0 until K[j]) {
                for (nj in 0..min(N, F[j])) {
                    for (i in 0 until M) {
                        if (i == j) continue
                        for (m in 0 until MR) {
                            val coeffs = DoubleArray(nVars)
                            val idxDiag = p2Idx(j, nj, kj, j, nj, kj, m)
                            if (idxDiag >= 0) coeffs[idxDiag] = 1.0
                            for (ni in 0..min(N - nj, F[i])) {
                                for (hi in 0 until K[i]) {
                                    val idx = p2Idx(j, nj, kj, i, ni, hi, m)
                                    if (idx >= 0) coeffs[idx] -= 1.0
                                }
                            }
                            constraints.add(LinearConstraint(coeffs, Relationship.EQ, 0.0))
                        }
                    }
                }
            }
        }

        // UEFF: e(i,ki) = sum of p2 where queue i is not blocked
        for (i in 0 until M) {
            for (ki in 0 until K[i]) {
                val coeffs = DoubleArray(nVars)
                val idxE = eIdx(i, ki)
                if (idxE >= 0) coeffs[idxE] = -1.0
                for (j in 0 until M) {
                    for (nj in 0..min(N, F[j])) {
                        for (kj in 0 until K[j]) {
                            for (m in 0 until MR) {
                                if (BB.get(m, i).toInt() == 0) {
                                    for (ni in 1..min(N, F[i])) {
                                        val idx = p2Idx(j, nj, kj, i, ni, ki, m)
                                        if (idx >= 0) coeffs[idx] += 1.0
                                    }
                                }
                            }
                        }
                    }
                }
                constraints.add(LinearConstraint(coeffs, Relationship.EQ, 0.0))
            }
        }

        // THM1: Phase balance
        for (i in 0 until M) {
            for (ki in 0 until K[i]) {
                val coeffs = DoubleArray(nVars)
                // LHS
                for (j in 0 until M) {
                    for (hi in 0 until K[i]) {
                        if (j != i || hi != ki) {
                            val idxE = eIdx(i, ki)
                            if (idxE >= 0) coeffs[idxE] += q[i][j].get(ki, hi)
                        }
                    }
                }
                // RHS (subtract)
                for (j in 0 until M) {
                    for (hi in 0 until K[i]) {
                        if (j != i || hi != ki) {
                            val idxE = eIdx(i, hi)
                            if (idxE >= 0) coeffs[idxE] -= q[i][j].get(hi, ki)
                        }
                    }
                }
                constraints.add(LinearConstraint(coeffs, Relationship.EQ, 0.0))
            }
        }

        // THM2: Population constraint
        for (j in 0 until M) {
            for (kj in 0 until K[j]) {
                for (nj in 0..F[j]) {
                    for (m in 0 until MR) {
                        val coeffs = DoubleArray(nVars)
                        val idxDiag = p2Idx(j, nj, kj, j, nj, kj, m)
                        if (idxDiag >= 0) coeffs[idxDiag] = -N.toDouble()
                        for (i in 0 until M) {
                            for (ni in 1..F[i]) {
                                for (ki in 0 until K[i]) {
                                    val idx = p2Idx(j, nj, kj, i, ni, ki, m)
                                    if (idx >= 0) coeffs[idx] += ni.toDouble()
                                }
                            }
                        }
                        constraints.add(LinearConstraint(coeffs, Relationship.EQ, 0.0))
                    }
                }
            }
        }

        // COR1: Second moment constraint
        val cor1Coeffs = DoubleArray(nVars)
        for (m in 0 until MR) {
            for (i in 0 until M) {
                for (j in 0 until M) {
                    for (nj in 1..F[j]) {
                        for (ni in 1..F[i]) {
                            for (ki in 0 until K[i]) {
                                for (kj in 0 until K[j]) {
                                    val idx = p2Idx(j, nj, kj, i, ni, ki, m)
                                    if (idx >= 0) cor1Coeffs[idx] += (ni * nj).toDouble()
                                }
                            }
                        }
                    }
                }
            }
        }
        constraints.add(LinearConstraint(cor1Coeffs, Relationship.EQ, (N * N).toDouble()))

        // THM30: Marginal balance for ni=0 (per phase), i<>f
        for (i in 0 until M) {
            if (i == f) continue
            for (ui in 0 until K[i]) {
                val coeffs = DoubleArray(nVars)
                // LHS: arrivals from j<>i,j<>f with BB(m,j)==0
                for (j in 0 until M) {
                    if (j == i || j == f) continue
                    for (nj in 1..F[j]) {
                        for (kj in 0 until K[j]) {
                            for (hj in 0 until K[j]) {
                                for (m in 0 until MR) {
                                    if (BB.get(m, j).toInt() == 0) {
                                        val idx = p2Idx(j, nj, kj, i, 0, ui, m)
                                        if (idx >= 0) coeffs[idx] += q[j][i].get(kj, hj)
                                    }
                                }
                            }
                        }
                    }
                }
                // LHS: arrivals from j==f with MM(m,0)<>i
                for (nj in 1..F[f]) {
                    for (kj in 0 until K[f]) {
                        for (hj in 0 until K[f]) {
                            for (m in 0 until MR) {
                                if (MM.get(m, 0).toInt() != i) {
                                    val idx = p2Idx(f, nj, kj, i, 0, ui, m)
                                    if (idx >= 0) coeffs[idx] += q[f][i].get(kj, hj)
                                }
                            }
                        }
                    }
                }

                // RHS: departures from i at ni=1 to j<>i,j<>f with BB(m,i)==0
                for (j in 0 until M) {
                    if (j == i || j == f) continue
                    for (nj in 0..F[j]) {
                        for (ki in 0 until K[i]) {
                            for (hj in 0 until K[j]) {
                                for (m in 0 until MR) {
                                    if (BB.get(m, i).toInt() == 0) {
                                        val idx = p2Idx(j, nj, hj, i, 1, ki, m)
                                        if (idx >= 0) coeffs[idx] -= q[i][j].get(ki, ui)
                                    }
                                }
                            }
                        }
                    }
                }
                // RHS: departures to j==f with BB(m,i)==0 and nj<F(f)
                for (nj in 0 until F[f]) {
                    for (ki in 0 until K[i]) {
                        for (hj in 0 until K[f]) {
                            for (m in 0 until MR) {
                                if (BB.get(m, i).toInt() == 0) {
                                    val idx = p2Idx(f, nj, hj, i, 1, ki, m)
                                    if (idx >= 0) coeffs[idx] -= q[i][f].get(ki, ui)
                                }
                            }
                        }
                    }
                }
                // RHS: unblocking when BB(m,i)==1 and MM(m,0)==i
                for (m in 0 until MR) {
                    if (BB.get(m, i).toInt() == 1 && MM.get(m, 0).toInt() == i) {
                        for (kf in 0 until K[f]) {
                            for (pf in 0 until K[f]) {
                                for (w in 0 until M) {
                                    if (w != f && w != i) {
                                        val idx = p2Idx(f, F[f], kf, i, 1, ui, m)
                                        if (idx >= 0) coeffs[idx] -= q[f][w].get(kf, pf)
                                    }
                                }
                            }
                        }
                    }
                }
                constraints.add(LinearConstraint(coeffs, Relationship.EQ, 0.0))
            }
        }

        // THM3: Marginal balance for ni in 1:F(i)-1, i<>f
        for (i in 0 until M) {
            if (i == f) continue
            for (ni in 1 until F[i]) {
                val coeffs = DoubleArray(nVars)
                // LHS: arrivals from j<>i,j<>f with BB(m,j)==0
                for (j in 0 until M) {
                    if (j == i || j == f) continue
                    for (nj in 1..F[j]) {
                        for (kj in 0 until K[j]) {
                            for (hj in 0 until K[j]) {
                                for (ui in 0 until K[i]) {
                                    for (m in 0 until MR) {
                                        if (BB.get(m, j).toInt() == 0) {
                                            val idx = p2Idx(j, nj, kj, i, ni, ui, m)
                                            if (idx >= 0) coeffs[idx] += q[j][i].get(kj, hj)
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                // LHS: arrivals from j==f with MM(m,0)<>i
                for (nj in 1..F[f]) {
                    for (kj in 0 until K[f]) {
                        for (hj in 0 until K[f]) {
                            for (ui in 0 until K[i]) {
                                for (m in 0 until MR) {
                                    if (MM.get(m, 0).toInt() != i) {
                                        val idx = p2Idx(f, nj, kj, i, ni, ui, m)
                                        if (idx >= 0) coeffs[idx] += q[f][i].get(kj, hj)
                                    }
                                }
                            }
                        }
                    }
                }

                // RHS: departures from i at ni+1 to j<>i,j<>f with BB(m,i)==0
                for (j in 0 until M) {
                    if (j == i || j == f) continue
                    for (nj in 0..F[j]) {
                        for (ki in 0 until K[i]) {
                            for (hi in 0 until K[i]) {
                                for (uj in 0 until K[j]) {
                                    for (m in 0 until MR) {
                                        if (BB.get(m, i).toInt() == 0) {
                                            val idx = p2Idx(j, nj, uj, i, ni + 1, ki, m)
                                            if (idx >= 0) coeffs[idx] -= q[i][j].get(ki, hi)
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                // RHS: departures to j==f with BB(m,i)==0 and nj<F(f)
                for (nj in 0 until F[f]) {
                    for (ki in 0 until K[i]) {
                        for (hi in 0 until K[i]) {
                            for (uj in 0 until K[f]) {
                                for (m in 0 until MR) {
                                    if (BB.get(m, i).toInt() == 0) {
                                        val idx = p2Idx(f, nj, uj, i, ni + 1, ki, m)
                                        if (idx >= 0) coeffs[idx] -= q[i][f].get(ki, hi)
                                    }
                                }
                            }
                        }
                    }
                }
                // RHS: unblocking when BB(m,i)==1 and MM(m,0)==i
                for (m in 0 until MR) {
                    if (BB.get(m, i).toInt() == 1 && MM.get(m, 0).toInt() == i) {
                        for (ki in 0 until K[i]) {
                            for (kf in 0 until K[f]) {
                                for (pf in 0 until K[f]) {
                                    for (w in 0 until M) {
                                        if (w != f && w != i) {
                                            val idx = p2Idx(f, F[f], kf, i, ni + 1, ki, m)
                                            if (idx >= 0) coeffs[idx] -= q[f][w].get(kf, pf)
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                constraints.add(LinearConstraint(coeffs, Relationship.EQ, 0.0))
            }
        }

        // THM3f: Marginal balance for i==f
        for (ni in 0 until F[f]) {
            val coeffs = DoubleArray(nVars)
            // LHS: arrivals from j<>f with BB(m,j)==0 and ni<F(f)
            if (ni < F[f]) {
                for (j in 0 until M) {
                    if (j == f) continue
                    for (nj in 1..F[j]) {
                        for (kj in 0 until K[j]) {
                            for (hj in 0 until K[j]) {
                                for (uf in 0 until K[f]) {
                                    for (m in 0 until MR) {
                                        if (BB.get(m, j).toInt() == 0) {
                                            val idx = p2Idx(j, nj, kj, f, ni, uf, m)
                                            if (idx >= 0) coeffs[idx] += q[j][f].get(kj, hj)
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            // RHS: departures from f at ni+1 to j<>f (only m=0, no blocking)
            for (j in 0 until M) {
                if (j == f) continue
                for (nj in 0..F[j]) {
                    for (kf in 0 until K[f]) {
                        for (hf in 0 until K[f]) {
                            for (uj in 0 until K[j]) {
                                if (ni < F[f]) {
                                    val idx = p2Idx(j, nj, uj, f, ni + 1, kf, 0)
                                    if (idx >= 0) coeffs[idx] -= q[f][j].get(kf, hf)
                                }
                            }
                        }
                    }
                }
            }
            constraints.add(LinearConstraint(coeffs, Relationship.EQ, 0.0))
        }

        // THM3I: Blocking depth balance
        for (z in 0 until ZM) {
            val coeffs = DoubleArray(nVars)
            // LHS: arrivals to f at F(f) from j<>f with BB(m,j)==0 and ZZ(m)==z
            for (j in 0 until M) {
                if (j == f) continue
                for (nj in 1..F[j]) {
                    for (kj in 0 until K[j]) {
                        for (hj in 0 until K[j]) {
                            for (uf in 0 until K[f]) {
                                for (m in 0 until MR) {
                                    if (BB.get(m, j).toInt() == 0 && ZZ[m] == z) {
                                        val idx = p2Idx(j, nj, kj, f, F[f], uf, m)
                                        if (idx >= 0) coeffs[idx] += q[j][f].get(kj, hj)
                                    }
                                }
                            }
                        }
                    }
                }
            }
            // RHS: departures from f to j<>f with ZZ(m)==z+1
            for (j in 0 until M) {
                if (j == f) continue
                for (nj in 0..F[j]) {
                    for (kf in 0 until K[f]) {
                        for (hf in 0 until K[f]) {
                            for (uj in 0 until K[j]) {
                                for (m in 0 until MR) {
                                    if (ZZ[m] == z + 1) {
                                        val idx = p2Idx(j, nj, uj, f, F[f], kf, m)
                                        if (idx >= 0) coeffs[idx] -= q[f][j].get(kf, hf)
                                    }
                                }
                            }
                        }
                    }
                }
            }
            constraints.add(LinearConstraint(coeffs, Relationship.EQ, 0.0))
        }

        // THM3L: Maximum blocking depth constraint
        for (m in 0 until MR) {
            if (ZZ[m] != ZM - 1) continue
            val coeffs = DoubleArray(nVars)
            // LHS: arrivals from j<>f with BB(m,j)==0 and MM1(m,j)>0
            for (j in 0 until M) {
                if (j == f || BB.get(m, j).toInt() != 0 || MM1.get(m, j).toInt() <= 0) continue
                for (nj in 1..F[j]) {
                    for (kj in 0 until K[j]) {
                        for (hj in 0 until K[j]) {
                            for (uf in 0 until K[f]) {
                                val idx = p2Idx(j, nj, kj, f, F[f], uf, m)
                                if (idx >= 0) coeffs[idx] += q[j][f].get(kj, hj)
                            }
                        }
                    }
                }
            }
            // RHS: uses MM1(m,j) to index into blocking configuration
            for (j in 0 until M) {
                if (j == f || BB.get(m, j).toInt() != 0 || MM1.get(m, j).toInt() <= 0) continue
                val mp = MM1.get(m, j).toInt() - 1 // Convert to 0-based
                for (kf in 0 until K[f]) {
                    for (uf in 0 until K[f]) {
                        for (w in 0 until M) {
                            if (w != f) {
                                val idx = p2Idx(f, F[f], kf, f, F[f], kf, mp)
                                if (idx >= 0) coeffs[idx] -= q[f][w].get(kf, uf)
                            }
                        }
                    }
                }
            }
            constraints.add(LinearConstraint(coeffs, Relationship.EQ, 0.0))
        }

        // THM4: Queue-length bound inequality
        for (j in 0 until M) {
            for (kj in 0 until K[j]) {
                for (i in 0 until M) {
                    for (m in 0 until MR) {
                        val coeffs = DoubleArray(nVars)
                        // LHS: sum_t sum_ht sum_nj sum_nt nt * p2
                        for (t in 0 until M) {
                            for (ht in 0 until K[t]) {
                                for (njt in 0..F[j]) {
                                    for (nt in 1..F[t]) {
                                        val idx = p2Idx(j, njt, kj, t, nt, ht, m)
                                        if (idx >= 0) coeffs[idx] += nt.toDouble()
                                    }
                                }
                            }
                        }
                        // RHS: -N * sum
                        for (hi in 0 until K[i]) {
                            for (njt in 0..F[j]) {
                                for (ni in 1..F[i]) {
                                    val idx = p2Idx(j, njt, kj, i, ni, hi, m)
                                    if (idx >= 0) coeffs[idx] -= N.toDouble()
                                }
                            }
                        }
                        constraints.add(LinearConstraint(coeffs, Relationship.GEQ, 0.0))
                    }
                }
            }
        }

        // Build objective function
        val objectiveCoeffs = DoubleArray(nVars)
        val targetQueue = objectiveQueue - 1
        for (m in 0 until MR) {
            for (ki in 0 until K[targetQueue]) {
                for (ni in 1..F[targetQueue]) {
                    val idx = p2Idx(targetQueue, ni, ki, targetQueue, ni, ki, m)
                    if (idx >= 0) objectiveCoeffs[idx] = 1.0
                }
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

            // Extract utilizations
            for (i in 0 until M) {
                var totalU = 0.0
                for (m in 0 until MR) {
                    for (ki in 0 until K[i]) {
                        for (ni in 1..F[i]) {
                            val idx = p2Idx(i, ni, ki, i, ni, ki, m)
                            if (idx >= 0) totalU += solution.point[idx]
                        }
                    }
                }
                variables["U_${i + 1}"] = totalU
            }

            // Extract effective utilizations
            for (i in 0 until M) {
                var totalE = 0.0
                for (ki in 0 until K[i]) {
                    val idx = eIdx(i, ki)
                    if (idx >= 0) {
                        totalE += solution.point[idx]
                        variables["e_${i + 1}_${ki + 1}"] = solution.point[idx]
                    }
                }
                variables["Ueff_${i + 1}"] = totalE
                variables["pb_${i + 1}"] = variables["U_${i + 1}"]!! - totalE
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
