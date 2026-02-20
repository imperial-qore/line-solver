package jline.api.qsys

import jline.io.Ret.qsys
import org.apache.commons.math3.util.FastMath

/**
 * Analyzes a G/G/1 queueing system using Kobayashi's approximation.
 *
 * @param lambda Arrival rate.
 * @param mu     Service rate.
 * @param ca     Coefficient of variation of the arrival process.
 * @param cs     Coefficient of variation of the service time.
 * @return qsysReturn containing average waiting time (W) and modified utilization (rhohat).
 */
fun qsys_gig1_approx_kobayashi(lambda: Double, mu: Double, ca: Double, cs: Double): qsys {
    val rho = lambda / mu
    val rhohat = FastMath.exp(-2 * (1 - rho) / (rho * (ca * ca + cs * cs / rho)))
    val W = rhohat / (1 - rhohat) / lambda
    return qsys(W, rhohat)
}
/**
 * Queueing system gig1 approx kobayashi algorithms
 */
@Suppress("unused")
class QsysGig1ApproxKobayashiAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}