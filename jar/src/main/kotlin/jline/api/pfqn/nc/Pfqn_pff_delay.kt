/**
 * @file Product-form factor computation for delay stations
 * 
 * Computes the product-form factor for delay stations in closed queueing networks.
 * Calculates the term Z[k]^n[k]/n[k]! for each class k, which represents the contribution
 * of delay stations to the normalizing constant in product-form networks.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.nc

import jline.util.Maths
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath

/**
 * Compute the product-form factor relatively to a Delay station
 *
 * @param Z - think times at the Delay station
 * @param n - number of jobs for each class
 * @return product of terms Z[k]^n[k]/n[k]! for all classes k
 */
fun pfqn_pff_delay(Z: Matrix, n: Matrix): Double {
    val R = n.length()
    if (n.sumRows().sumCols().value() == 0.0) {
        return 1.0
    }

    var f = 0.0
    for (r in 0..<R) {
        if (Z[r] > 0) {
            f += FastMath.log(Z[r]) * n[r]
            f -= Maths.factln(n[r])
        } else if (n[r] > 0) {
            return 0.0
        }
    }
    return FastMath.exp(f)
}
/**
 * PFQN pff delay algorithms
 */
@Suppress("unused")
class PfqnPffDelayAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}