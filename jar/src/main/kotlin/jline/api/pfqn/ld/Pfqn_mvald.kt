/**
 * Load-Dependent Mean Value Analysis
 * 
 * Implements MVA for closed networks with load-dependent service rates where service
 * capacity varies with queue length. Handles complex systems like multi-server stations
 * and state-dependent service disciplines with optional numerical stabilization.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.ld

import jline.io.Ret
import jline.util.PopulationLattice
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath

/**
 * Mean Value Analysis (MVA) method for closed networks with load dependent service and stabilization
 */
fun pfqn_mvald(L: Matrix, N: Matrix, Z: Matrix, mu: Matrix, stabilize: Boolean = true): Ret.pfqnMVALD {
    var warn = true
    var isNumStable = true
    val M = L.numRows // number of queues
    val R = L.numCols // number of classes
    var prodN = 1.0
    for (i in 0..<N.numRows) {
        for (j in 0..<N.numCols) {
            prodN *= (1 + N[i, j])
        }
    }
    val Xs = Matrix(R, prodN.toInt()) // throughput for a model with station i less
    val pi = arrayOfNulls<Matrix>(M) // marginal queue-length probabilities pi(k)
    for (ist in 0..<M) {
        pi[ist] = Matrix.ones(N.elementSum().toInt() + 1, prodN.toInt())
    }
    var WN = Matrix(M, R)
    var n = PopulationLattice.pprod(N)
    val lGN: MutableList<Double> = ArrayList()
    lGN.add(0.0)
    while (!(n.numRows == 1 && n.numCols == 1 && n.value() == -1.0)) {
        WN = Matrix(WN.numRows, WN.numCols)
        for (s in 0..<R) {
            if (n[s] > 0) {
                for (ist in 0..<M) {
                    WN[ist, s] = 0
                    var k = 0
                    while (k < n.elementSum()) {
                        WN[ist, s] =
                            WN[ist, s] + ((L[ist, s] / mu[ist, k]) * (k + 1) * pi[ist]!![k, PopulationLattice.hashpop(Matrix.oner(
                                n,
                                ArrayList(listOf(s))), N)])
                        k++
                    }
                }
                Xs[s, PopulationLattice.hashpop(n, N)] = n[s] / (Z[s] + Matrix.extractColumn(WN, s, null).elementSum())
            }
        }
        // Compute pi(k|n)
        var k = 0
        while (k < n.elementSum()) {
            for (ist in 0..<M) {
                pi[ist]!![k + 1, PopulationLattice.hashpop(n, N)] = 0
            }
            for (s in 0..<R) {
                if (n[s] > 0) {
                    for (ist in 0..<M) {
                        pi[ist]!![k + 1, PopulationLattice.hashpop(n, N)] = (pi[ist]!![k + 1, PopulationLattice.hashpop(n,
                            N)] + ((L[ist, s] / mu[ist, k]) * Xs[s, PopulationLattice.hashpop(n,
                            N)] * pi[ist]!![k, PopulationLattice.hashpop(Matrix.oner(n, ArrayList(listOf(s))), N)]))
                    }
                }
            }
            k++
        }
        // Compute pi(0|n)
        for (ist in 0..<M) {
            var sumpi = 0.0
            var k = 0
            while (k < n.elementSum()) {
                sumpi += pi[ist]!![k + 1, PopulationLattice.hashpop(n, N)]
                k++
            }
            val p0 = 1 - sumpi
            if (p0 < 0) {
                if (warn) {
                    System.err.println("pfqn_mvald: MVA-LD is numerically unstable on this model, " + "LINE will force all probabilities to be non-negative.")
                    warn = false
                    isNumStable = false
                }
                if (stabilize) {
                    pi[ist]!![0, PopulationLattice.hashpop(n, N)] = FastMath.ulp(1.0)
                } else {
                    pi[ist]!![0, PopulationLattice.hashpop(n, N)] = p0
                }
            } else {
                pi[ist]!![0, PopulationLattice.hashpop(n, N)] = p0
            }
        }
        val nfind = n.find()
        if (!nfind.isEmpty) {
            var idx = nfind.numRows - 1
            while (idx >= 0 && n[nfind[idx].toInt()] < 0) {
                idx--
            }
            val last_nnz = nfind[idx].toInt()
            var sumn = 0.0
            var sumN = 0.0
            var sumnp = 0.0
            for (i in 0..<R) {
                if (i < last_nnz) {
                    sumn += n[i]
                    sumN += N[i]
                } else if (i > last_nnz) {
                    sumnp += n[i]
                }
            }
            if (sumn == sumN && sumnp == 0.0) {
                val logX = FastMath.log(Xs[last_nnz, PopulationLattice.hashpop(n, N)])
                lGN.add(lGN[lGN.size - 1] - logX)
            }
        }
        n = PopulationLattice.pprod(n, N) // get the next population
    }
    val XN = Matrix.extractColumn(Xs, PopulationLattice.hashpop(N, N), null).transpose()
    val newpi = Matrix(M, N.elementSum().toInt() + 1)
    for (i in 0..<M) {
        for (j in 0..<N.elementSum().toInt() + 1) {
            newpi[i, j] = pi[i]!![j, PopulationLattice.hashpop(N, N)]
        }
    }
    val QN = if (WN.isEmpty) {
        Matrix(0, 0)
    } else {
        WN.elementMult(XN.repmat(M, 1), null)
    }
    val UN = Matrix.ones(newpi.numRows, 1).sub(Matrix.extractColumn(newpi, 0, null))
    val CN = Matrix(N.numRows, N.numCols)
    for (i in 0..<CN.numRows) {
        for (j in 0..<CN.numCols) {
            CN[i, j] = N[i, j] / XN[i, j] - Z[i, j]
        }
    }
    return Ret.pfqnMVALD(XN, QN, UN, CN, lGN, isNumStable, newpi)
}
/**
 * PFQN mvald algorithms
 */
@Suppress("unused")
class PfqnMvaldAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}