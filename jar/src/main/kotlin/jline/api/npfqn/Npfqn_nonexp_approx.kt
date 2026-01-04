/**
 * @file NPFQN Non-Exponential Approximation
 * 
 * Implements approximation methods for non-product-form queueing networks
 * with non-exponential service and inter-arrival time distributions.
 * Provides analysis techniques for NPFQN systems beyond exponential assumptions.
 * 
 * @since LINE 3.0
 */
package jline.api.npfqn

import jline.io.Ret.npfqnNonexpApprox
import jline.lang.NetworkStruct
import jline.GlobalConstants
import jline.lang.constant.SchedStrategy
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath
import kotlin.math.pow

/**
 * APIs for evaluating non-product-form queueing networks.
 */

/**
 * Merges MMAP traffic flows with class switching.
 *
 * @param MMAPs  Map of MMAP traffic flows
 * @param prob   Probability matrix for class switching
 * @param config Configuration for merging ("default", "super")
 * @return Merged and normalized MMAP traffic flow
 */

/**
 * Approximates non-product-form queueing networks using the specified method.
 *
 * @param method   Approximation method ("default", "none", "hmva", "interp")
 * @param sn       Network structure
 * @param ST       Service time matrix
 * @param V        Visit ratios matrix
 * @param SCV      Squared coefficient of variation matrix
 * @param Tin      Initial throughput matrix
 * @param Uin      Initial utilization matrix
 * @param gamma    Gamma matrix for correction
 * @param nservers Number of servers matrix
 * @return Result containing updated matrices and additional parameters
 */
fun npfqn_nonexp_approx(method: String,
                        sn: NetworkStruct,
                        ST: Matrix,
                        V: Matrix?,
                        SCV: Matrix,
                        Tin: Matrix,
                        Uin: Matrix,
                        gamma: Matrix,
                        nservers: Matrix): npfqnNonexpApprox {
    val M = sn.nstations
    val rho = Matrix(M, 1)
    rho.zero()
    val scva = Matrix.ones(M, 1)
    val scvs = Matrix.ones(M, 1)
    val eta = Matrix.ones(M, 1)
    val T = Tin.copy()
    val U = Uin.copy()
    val STLocal = ST.copy()  // Create local copy to avoid modifying the original
    val gammaLocal = gamma.copy()  // Create local copy to avoid modifying the original
    val nserversLocal = nservers.copy()  // Create local copy to avoid modifying the original
    when (method) {
        "default", "none", "hmva" -> {
            // no-op - just return the initialized values as in MATLAB
            // Note: rho stays as zeros, but eta should remain as ones
            return npfqnNonexpApprox(STLocal, gammaLocal, nserversLocal, rho, scva, scvs, eta)
        }
        "interp" -> {
            var ist = 0
            while (ist < M) {
                val nnzClasses = Matrix(1, ST.numCols)
                nnzClasses.zero()
                run {
                    var j = 0
                    while (j < STLocal.numCols) {
                        if (java.lang.Double.isFinite(STLocal[ist, j]) && java.lang.Double.isFinite(SCV[ist, j])) {
                            nnzClasses[0, j] = 1
                        }
                        j++
                    }
                }
                var j = 0
                while (j < nnzClasses.numElements) {
                    if (nnzClasses[j] > 0) {
                        rho[ist] = rho[ist] + U[ist, j]
                    }
                    j++
                }
                if (nnzClasses.elementSum() != 0.0) {
                    when (sn.sched[sn.stations[ist]]) {
                        SchedStrategy.FCFS -> {
                            val STinnz = Matrix(1, nnzClasses.elementSum().toInt())
                            val SCVinnz = STinnz.copy()
                            val Tinnz = STinnz.copy()
                            var tempj = 0
                            var j = 0
                            while (j < nnzClasses.numElements) {
                                if (nnzClasses[j] > 0) {
                                    STinnz[tempj] = STLocal[ist, j]
                                    SCVinnz[tempj] = SCV[ist, j]
                                    Tinnz[tempj] = T[ist, j]
                                    tempj++
                                }
                                j++
                            }
                            if (STinnz.elementMax() - STinnz.elementMin() > 0 || SCVinnz.elementMax() > 1 + GlobalConstants.FineTol || SCVinnz.elementMin() < 1 - GlobalConstants.FineTol) {
                                scva[ist] = 1.0
                                scvs[ist] = SCVinnz.mult(Tinnz.transpose()).toDouble() / Tinnz.elementSum()
                                gammaLocal[ist] = (rho[ist].pow(nserversLocal[ist]) + rho[ist]) / 2
                                if (scvs[ist] > 1 - 1e-6 && scvs[ist] < 1 + 1e-6 && nserversLocal[ist] == 1.0) {
                                    eta[ist] = rho[ist]
                                } else {
                                    eta[ist] = FastMath.exp(-2 * (1 - rho[ist]) / (scvs[ist] + scva[ist] * rho[ist]))
                                }
                                val order = 8
                                val ai = FastMath.pow(rho[ist], order)
                                val bi = FastMath.pow(rho[ist], order)
                                var k = 0
                                while (k < nnzClasses.numElements) {
                                    if (nnzClasses[k] > 0 && sn.rates[ist, k] > 0) {
                                        STLocal[ist, k] =
                                            FastMath.max(0.0, 1 - ai) * STLocal[ist, k] + ai * (bi * eta[ist] + FastMath.max(0.0,
                                                1 - bi) * gammaLocal[ist]) * (nserversLocal[ist] / Tinnz.elementSum())
                                    }
                                    k++
                                }
                                nserversLocal[ist] = 1.0
                            }
                        }
                        else -> {} // no-op
                    }
                }
                ist++
            }
        }

        else -> throw IllegalArgumentException("Unknown approximation method: $method")
    }
    return npfqnNonexpApprox(STLocal, gammaLocal, nserversLocal, rho, scva, scvs, eta)
}