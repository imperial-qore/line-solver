package jline.api.qsys

import jline.io.Ret.qsys
import org.apache.commons.math3.util.FastMath

/**
 * Analyzes a G/G/1 queueing system using the Kramer-Langenbach-Belz (KLB) approximation.
 *
 * @param lambda Arrival rate.
 * @param mu     Service rate.
 * @param ca     Coefficient of variation of the arrival process.
 * @param cs     Coefficient of variation of the service time.
 * @return qsysReturn containing average waiting time (W) and modified utilization (rhohat).
 */
fun qsys_gig1_approx_klb(lambda: Double, mu: Double, ca: Double, cs: Double): qsys {
    // Kramer-Langenbach-Belz formula
    val rho = lambda / mu
    val g = if (ca <= 1) {
        FastMath.exp(-2 * (1 - rho) * FastMath.pow(1 - ca * ca, 2) / (3 * rho * (ca * ca + cs * cs)))
    } else {
        FastMath.exp(-(1 - rho) * (ca * ca - 1) / (ca * ca + 4 * cs * cs))
    }
    val W = 1.0 / mu * ((rho / (1 - rho)) * ((cs * cs + ca * ca) / 2) * g + 1)
    val rhohat = W * lambda / (1 + W * lambda)
    return qsys(W, rhohat)
}
/**
 * Queueing system gig1 approx klb algorithms
 */
@Suppress("unused")
class QsysGig1ApproxKlbAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}