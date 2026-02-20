/**
 * @file Loss Network Analysis via Erlang Fixed Point
 * 
 * Implements fixed-point algorithms for analyzing loss networks using Erlang
 * formulas. Provides efficient computation of blocking probabilities and
 * performance metrics for circuit-switched networks with finite capacity.
 * 
 * @since LINE 3.0
 */
package jline.api.lossn

import jline.GlobalConstants
import jline.io.Ret
import jline.util.Maths
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath
import org.apache.commons.math3.util.MathArrays
import java.util.*

/**
 * This method calculates the Erlang fixed point approximation for loss networks.
 *
 *
 * Calls (jobs) on route (class) r arrive according to a Poisson rate nu_r, r=1,...,R.
 * Call service times on route r have unit mean.
 *
 *
 * The link capacity requirements are defined as:
 * sum<sub>r</sub> A(j, r) n(j, r) &lt; C(j)
 * for all links j=1,...,J, where n(j, r) counts the calls on route r on link j.
 *
 * @param nuVec A Matrix (1xR) representing the arrival rate of route (class) r = 1,...,R.
 * @param Amat  A Matrix (J,R) representing the capacity requirement of link j for route r = 1,...,R.
 * @param cVec  A Matrix (J,1) representing the available capacity of link j.
 * @return lossnErlangFPReturn which contains:
 * - qLen (1xR): mean queue-length for route r = 1,...,R calls.
 * - loss (1xR): loss probability for route r = 1,...,R calls.
 * - eBlock (Jx1): blocking probability for link j = 1,...,J.
 *
 *
 * Note: nu_r may be replaced by a utilization rho_r=nu_r/mu_r, where mu_r is the service rate for route r.
 * Example: lossn_erlangfp(new Matrix("[0.3,0.1]"), new Matrix("[1,1;1,4]"), new Matrix("[1,3]")).qLen.print();
 */
fun lossn_erlangfp(nuVec: Matrix, Amat: Matrix, cVec: Matrix): Ret.lossnErlangFP {
    val nu = nuVec.toArray1D()
    val A = Amat.toArray2D()
    val c = cVec.toArray1D()
    val R = nu.size
    val J = c.size
    val E = DoubleArray(J)
    Arrays.fill(E, 0.5)
    val E_1 = DoubleArray(J)
    Arrays.fill(E_1, 0.0)
    var niter = 0

    while (MathArrays.distance(E, E_1) > 1e-8) {
        niter++
        System.arraycopy(E, 0, E_1, 0, J)

        for (j in 0..<J) {
            var rhoj_1 = 0.0

            for (r in 0..<R) {
                if (A[j][r] > 0) {
                    var termj = nu[r] * A[j][r]

                    for (i in 0..<J) {
                        if (A[i][r] > GlobalConstants.Zero) {
                            termj *= FastMath.pow(1 - E_1[i], A[i][r])
                        }
                    }
                    rhoj_1 += termj
                }
            }
            rhoj_1 /= (1 - E_1[j])
            E[j] = ErlangB(rhoj_1, c[j])
        }
    }

    val QLen = nu.copyOf(nu.size)

    for (r in 0..<R) {
        for (j in 0..<J) {
            QLen[r] *= FastMath.pow(1 - E[j], A[j][r])
        }
    }

    for (i in QLen.indices) {
        QLen[i] = FastMath.max(QLen[i], 0.0) // Ensure no negative values
    }

    val Loss = DoubleArray(R)
    for (i in 0..<R) {
        Loss[i] = QLen[i] / nu[i]
    }

    return Ret.lossnErlangFP(Matrix(QLen), Matrix(Loss), Matrix(E), niter)
}

/**
 * Calculates the Erlang B formula for a given arrival rate and capacity.
 *
 *
 * The Erlang B formula is used to compute the blocking probability in a loss system
 * where calls arrive according to a Poisson process and there are a fixed number of servers.
 *
 * @param nu The traffic intensity or offered load, which is the product of the arrival rate
 * and the average service time.
 * @param C  The total capacity or number of servers.
 * @return The blocking probability, which is the probability that a call is blocked due to all
 * servers being busy.
 */
private fun ErlangB(nu: Double, C: Double): Double {
    var den = 0.0
    var i = 0
    while (i <= C) {
        den += FastMath.exp(i * FastMath.log(nu) - Maths.factln(i))
        i++
    }
    val blockProb = C * FastMath.log(nu) - Maths.factln(C) - FastMath.log(den)
    return FastMath.exp(blockProb)
}