/**
 * @file Multi-class repairman model sampling method for normalizing constants
 * 
 * Implements importance sampling for computing normalizing constants in multi-class
 * repairman models. Uses hybrid uniform and exponential sampling strategies with
 * logarithmic scaling for numerical stability in high-population regimes.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.nc

import jline.io.Ret
import jline.util.Maths
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath
import java.util.*
import kotlin.math.exp

/**
 * APIs for evaluating Product-Form Queueing Networks.
 */


fun pfqn_mmsample2(L: Matrix, N: Matrix, Z: Matrix, samples: Int): Ret.pfqnNc {
    val L_local = L.copy()
    val Z_local = Z.copy()
    val R = N.numElements
    val scaleFactor = 1e-7 + FastMath.min(L_local.elementMin(), Z_local.elementMin())
    L_local.scaleEq(1.0 / scaleFactor)
    Z_local.scaleEq(1.0 / scaleFactor)

    val c = 0.5
    val numSamples1 = FastMath.ceil(c * samples).toInt()
    val numSamples2 = FastMath.ceil(samples * (1 - c)).toInt()
    val v = DoubleArray(numSamples1 + numSamples2)
    val du = DoubleArray(numSamples1 + numSamples2)

    val rand = Random()
    for (i in 0..<numSamples1) {
        v[i] = rand.nextDouble() // we sample more below the mean of the exponential
    }
    val lv = Maths.logSpace(0.0, 5.0, numSamples2)
    if (numSamples2 >= 0) System.arraycopy(lv, 0, v, numSamples1, numSamples2)

    for (i in du.indices) {
        du[i] = if (i == 0) 0.0 else v[i] - v[i - 1]
    }

    // u  = repmat(v',1,R);
    val u = Matrix(v.size, R, v.size * R)
    for (i in v.indices) {
        for (j in 0..<R) {
            u[i, j] = v[i]
        }
    }

    val ZL = Matrix(v.size, R, v.size * R)
    for (r in 0..<R) {
        for (i in v.indices) {
            ZL[i, r] = FastMath.log((Z_local[0, r] + L_local[0, r]) * u[i, r])
        }
    }

    val lGmat = Matrix(v.size, 1)
    for (i in v.indices) {
        var sum = 0.0
        for (j in 0..<R) {
            sum += ZL[i, j] * N[j]
        }
        lGmat[i, 0] = du[i] - v[i] + sum
    }

    var lG = lGmat.elementMax() - N.factln().elementSum()
    lG = lG + N.elementSum() * FastMath.log(scaleFactor)

    return Ret.pfqnNc(exp(lG), lG)
}
/**
 * PFQN mmsample2 algorithms
 */
@Suppress("unused")
class PfqnMmsample2Algo {
    companion object {
        // Class documentation marker for Dokka
    }
}