/**
 * @file CTMC randomization (uniformization) method
 *
 * Converts a continuous-time Markov chain into an equivalent discrete-time chain
 * using the randomization technique. This transformation enables the use of DTMC
 * algorithms for CTMC analysis and is fundamental to numerical CTMC solution methods.
 *
 * @since LINE 3.0
 */
package jline.api.mc;

import jline.util.Pair;
import jline.util.matrix.Matrix;

public final class Ctmc_randomization {
    private Ctmc_randomization() {}

    /**
     * Convert a CTMC to a DTMC using randomization technique
     *
     * @param Q Infinitesimal generator matrix of the CTMC
     * @param q Optional uniformization rate. If not provided, defaults to 1.05*max(|Q|)
     * @return Pair containing the transition matrix P and the uniformization rate q
     */
    public static Pair<Matrix, Double> ctmc_randomization(Matrix Q, Double q) {
        double uniformizationRate;
        if (q != null) {
            uniformizationRate = q;
        } else {
            // see _kb/03-api-layer.md for rationale
            uniformizationRate = 1.05 * Q.elementMaxAbs();
        }

        int n = Q.getNumRows();
        Matrix I = Matrix.eye(n);
        Matrix P = Q.scale(1.0 / uniformizationRate).add(I);

        Matrix stochasticP = Dtmc_makestochastic.dtmc_makestochastic(P);

        return new Pair<Matrix, Double>(stochasticP, uniformizationRate);
    }

    public static Pair<Matrix, Double> ctmc_randomization(Matrix Q) {
        return ctmc_randomization(Q, null);
    }
}
