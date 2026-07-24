/**
 * @file Kullback-Leibler divergence metric
 *
 * Implements the Kullback-Leibler divergence D_KL(P||Q) = sum p_i log(p_i/q_i) measuring
 * the relative entropy between two probability distributions. A fundamental measure
 * in information theory for distribution comparison and model selection.
 *
 * @since LINE 3.0
 */
package jline.api.measures;

import jline.GlobalConstants;
import jline.util.matrix.Matrix;

public final class Ms_kullbackleibler {
    private Ms_kullbackleibler() {}

    /**
     * Kullback-Leibler divergence between two probability distributions.
     * Part of Shannon's entropy family.
     * Measures information lost when Q is used to approximate P.
     *
     * @param P exact probability distribution
     * @param Q model probability distribution
     * @return KL divergence
     */
    public static double ms_kullbackleibler(Matrix P, Matrix Q) {
        if (P.length() != Q.length()) {
            throw new IllegalArgumentException("Distributions must have the same size");
        }

        double kl = 0.0;
        for (int i = 0; i < P.length(); i++) {
            double pi = P.get(i);
            double qi = Q.get(i);
            if (pi > 0) {
                if (qi == 0.0) {
                    return GlobalConstants.Inf;
                }
                kl += pi * Math.log(pi / qi);
            }
        }

        return kl;
    }
}
