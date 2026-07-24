/**
 * @file Product-form factor computation for delay stations
 *
 * Computes the product-form factor for delay stations in closed queueing networks.
 * Calculates the term Z[k]^n[k]/n[k]! for each class k, which represents the contribution
 * of delay stations to the normalizing constant in product-form networks.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc;

import jline.util.Maths;
import jline.util.matrix.Matrix;
import org.apache.commons.math3.util.FastMath;

public final class Pfqn_pff_delay {
    private Pfqn_pff_delay() {}

    /**
     * Compute the product-form factor relatively to a Delay station.
     *
     * @param Z think times at the Delay station
     * @param n number of jobs for each class
     * @return product of terms Z[k]^n[k]/n[k]! for all classes k
     */
    public static double pfqn_pff_delay(Matrix Z, Matrix n) {
        int R = n.length();
        if (n.sumRows().sumCols().value() == 0.0) {
            return 1.0;
        }

        double f = 0.0;
        for (int r = 0; r < R; r++) {
            if (Z.get(r) > 0) {
                f += FastMath.log(Z.get(r)) * n.get(r);
                f -= Maths.factln((int) n.get(r));
            } else if (n.get(r) > 0) {
                return 0.0;
            }
        }
        return FastMath.exp(f);
    }
}
