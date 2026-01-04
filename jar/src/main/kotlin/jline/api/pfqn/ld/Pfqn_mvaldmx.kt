/**
 * @file Load-dependent MVA for mixed open-closed queueing networks
 * 
 * Implements Mean Value Analysis for mixed queueing networks with both open and closed classes
 * and load-dependent service rates. Handles complex systems with state-dependent service
 * mechanisms and provides exact solutions for mixed network topologies.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.ld

import jline.io.Ret
import jline.util.Maths
import jline.util.PopulationLattice
import jline.util.Utils
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath

/**
 * Mean Value Analysis (MVA) method for mixed queueing networks with load-dependent nodes.
 *
 * @param lambda - arrival rates
 * @param D      - service demand matrix
 * @param N      - population vector
 * @param Z      - think times
 * @param mu     - load dependent rates
 * @param S      - number of servers at each station
 * @return - the performance measures for the given network.
 */

fun pfqn_mvaldmx(lambda: Matrix, D: Matrix, N: Matrix, Z: Matrix, mu: Matrix?, S: Matrix?): Ret.pfqnMVALDMX {/*
         * Parameter S is not used at all
         * */
    var mu = mu
    var S = S
    var NfiniteSum = 0.0
    for (i in 0..<N.numRows) {
        for (j in 0..<N.numCols) {
            val num = N[i, j]
            if (java.lang.Double.isFinite(num)) {
                NfiniteSum += num
            }
        }
    }
    if (mu == null && S == null) {
        mu = Matrix.ones(D.numRows, NfiniteSum.toInt())
        S = Matrix.ones(D.numRows, 1)
    }
    if (mu!!.numCols < NfiniteSum) {
        throw RuntimeException("pfqn_mvaldmx: MVALDMX requires to specify the load-dependent rates with one job more than the maximum closed population.")
    }
    val f = lambda.find()
    for (i in 0..<f.numRows) {
        val num = N[f[i].toInt()]
        if (num > 0 && java.lang.Double.isFinite(num)) {
            throw RuntimeException("pfqn_mvaldmx: Arrival rate cannot be specified on closed classes.")
        }
    }
    val M = D.numRows
    val R = D.numCols
    val openClasses = ArrayList<Int>()
    val closedClasses = ArrayList<Int>()
    for (i in 0..<N.length()) {
        if (Utils.isInf(N[i])) {
            openClasses.add(i)
        } else {
            closedClasses.add(i)
        }
    }
    val XN = Matrix(1, R)
    val UN = Matrix(M, R)
    val CN = Matrix(M, R)
    val QN = Matrix(M, R)
    var lGN = 0.0
    val newMu = Matrix(mu.numRows, mu.numCols + 1)
    for (i in 0..<newMu.numRows) {
        for (j in 0..<newMu.numCols) {
            if (j < mu.numCols) {
                newMu[i, j] = mu[i, j]
            } else {
                newMu[i, j] = mu[i, j - 1]
            }
        }
    }
    mu = newMu // we need up to sum(N)+1, but there is limited load dep
    val ret1 = pfqn_mvaldmx_ec(lambda, D, Matrix(mu))
    val EC = ret1.EC
    val E = ret1.E
    val Eprime = ret1.Eprime
    val C = closedClasses.size // number of closed classes
    val Dc = Matrix(D.numRows, C)
    val Nc = Matrix(1, C)
    val Zc = Matrix(1, C)
    for (i in 0..<C) {
        val c = closedClasses[i]
        for (j in 0..<Dc.numRows) {
            Dc[j, i] = D[j, c]
        }
        Nc[0, i] = N[c]
        Zc[0, i] = Z[c]
    }
    val prods = Matrix(1, C) // needed for fast hashing
    for (r in 0..<C) {
        var prod = 1.0
        for (i in 0..<r) {
            prod *= (Nc[i] + 1)
        }
        prods[0, r] = prod
    }
    // Start at nc=(0,...,0)
    var nvec = PopulationLattice.pprod(Nc)
    // Initialize Pc
    val Pc = arrayOfNulls<Matrix>(M)
    var ncProd = 1.0
    for (i in 0..<C) {
        ncProd *= (1 + Nc[i])
    }
    for (ist in 0..<M) {
        Pc[ist] = Matrix((1 + Nc.elementSum()).toInt(), ncProd.toInt())
    }
    val x = Matrix(C, ncProd.toInt())
    val w = arrayOfNulls<Matrix>(M)
    for (ist in 0..<M) {
        w[ist] = Matrix(C, ncProd.toInt())
    }
    for (ist in 0..<M) {
        Pc[ist]!![0, PopulationLattice.hashpop(nvec, Nc, C, prods)] = 1
    }
    val u = Matrix(M, C)
    // Population recursion
    while (!(nvec.numRows == 1 && nvec.numCols == 1 && nvec.value() == -1.0)) {
        val hnvec = PopulationLattice.hashpop(nvec, Nc, C, prods)
        val nc = nvec.elementSum()
        for (ist in 0..<M) {
            for (c in 0..<C) {
                if (nvec[c] > 0) {
                    val hnvec_c = PopulationLattice.hashpop(Matrix.oner(nvec, ArrayList(listOf(c))), Nc, C, prods)
                    // Compute mean residence times
                    var n = 0
                    while (n < nc) {
                        w[ist]!![c, hnvec] = w[ist]!![c, hnvec] + Dc[ist, c] * (n + 1) * EC[ist, n] * Pc[ist]!![n, hnvec_c]
                        n++
                    }
                }
            }
        }
        // Compute tput
        for (c in 0..<C) {
            var sumw = 0.0
            for (ist in 0..<M) {
                sumw += w[ist]!![c, hnvec]
            }
            x[c, hnvec] = nvec[c] / (Zc[c] + sumw)
        }
        for (ist in 0..<M) {
            var n = 0
            while (n < nc) {
                for (c in 0..<C) {
                    if (nvec[c] > 0) {
                        val hnvec_c = PopulationLattice.hashpop(Matrix.oner(nvec, ArrayList(listOf(c))), Nc, C, prods)
                        Pc[ist]!![1 + n, hnvec] =
                            Pc[ist]!![1 + n, hnvec] + Dc[ist, c] * EC[ist, n] * x[c, hnvec] * Pc[ist]!![n, hnvec_c]
                    }
                }
                n++
            }
            var sumpc = 0.0
            var k = 0
            while (k < nc) {
                sumpc += Pc[ist]!![1 + k, hnvec]
                k++
            }
            Pc[ist]!![0, hnvec] = Maths.max(Math.ulp(1.0), 1 - sumpc)
        }
        val nvecFind = nvec.find()
        if (!nvecFind.isEmpty) {
            var idx = nvecFind.numRows - 1
            while (idx >= 0 && nvec[nvecFind[idx].toInt()] <= 0) {
                idx--
            }

            // now compute the normalizing constant
            val last_nnz = nvecFind[idx].toInt()
            var sumnvec: Double
            var sumnc: Double
            var sumnvecp = 0.0
            sumnc = sumnvecp
            sumnvec = sumnc
            for (i in 0..<C) {
                if (i < last_nnz) {
                    sumnvec += nvec[i]
                    sumnc += Nc[i]
                } else if (i > last_nnz) {
                    sumnvecp += nvec[i]
                }
            }
            if (sumnvec == sumnc && sumnvecp == 0.0) {
                val logX = FastMath.log(XN[last_nnz])
                lGN -= logX
            }
        }
        nvec = PopulationLattice.pprod(nvec, Nc)
    }
    // compute performance indexes at Nc for closed classes
    val hnvec = PopulationLattice.hashpop(Nc, Nc, C, prods)
    for (c in 0..<C) {
        val hnvec_c = PopulationLattice.hashpop(Matrix.oner(Nc, ArrayList(listOf(c))), Nc, C, prods)
        for (ist in 0..<M) {
            u[ist, c] = 0
            val sumNc = Nc.elementSum()
            var n = 0
            while (n < sumNc) {
                u[ist, c] = u[ist, c] + (Dc[ist, c] * x[c, hnvec] * Eprime[ist, n]) / E[ist, n] * Pc[ist]!![n, hnvec_c]
                n++
            }
        }
    }
    for (i in 0..<C) {
        // Throughput
        XN[closedClasses[i]] = x[i, hnvec]
        for (j in 0..<M) {
            // Utilization
            UN[j, closedClasses[i]] = u[j, i]
            // Response time
            CN[j, closedClasses[i]] = w[j]!![i, hnvec]
            // Queue length
            QN[j, closedClasses[i]] = XN[closedClasses[i]] * CN[j, closedClasses[i]]
        }
    }
    // Compute performance indexes at Nc for open classes
    for (r in openClasses) {
        XN[r] = lambda[r]
        for (ist in 0..<M) {
            // Queue-length
            QN[ist, r] = 0
            val sumNc = Nc.elementSum()
            run {
                var n = 0
                while (n <= sumNc) {
                    QN[ist, r] = QN[ist, r] + lambda[r] * D[ist, r] * (n + 1) * EC[ist, n] * Pc[ist]!![n, hnvec]
                    n++
                }
            }
            // Response time
            CN[ist, r] = QN[ist, r] / lambda[r]
            // Utilization - the formula from Bruell-Balbo-Ashfari does not
            // match simulation, this appears to be simply lambda_r*D_{ir}
            UN[ist, r] = 0
            var n = -1
            while (n < sumNc) {
                UN[ist, r] = (UN[ist, r] + lambda[r] * Eprime[ist, n + 2] / E[ist, n + 2] * Pc[ist]!![n + 1, hnvec])
                n++
            }
        }
    }
    val newPc = Matrix(M, (1 + Nc.elementSum()).toInt())
    for (ist in 0..<M) {
        for (j in 0..<newPc.numCols) {
            newPc[ist, j] = Pc[ist]!![j, hnvec]
        }
    }
    return Ret.pfqnMVALDMX(XN, QN, UN, CN, lGN, newPc)
}
/**
 * PFQN mvaldmx algorithms
 */
@Suppress("unused")
class PfqnMvaldmxAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}