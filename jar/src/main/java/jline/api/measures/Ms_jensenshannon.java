/**
 * @file Jensen-Shannon divergence metric
 *
 * Implements Jensen-Shannon divergence, a symmetric and bounded version of
 * Kullback-Leibler divergence. Measures the similarity between probability
 * distributions with values between 0 and 1, commonly used in phylogenetics and text analysis.
 *
 * @since LINE 3.0
 */
package jline.api.measures;

import jline.util.matrix.Matrix;

public final class Ms_jensenshannon {
    private Ms_jensenshannon() {}

    /**
     * Jensen-Shannon divergence between two probability distributions.
     * Part of Shannon's entropy family. Symmetric version of KL divergence.
     *
     * @param P first probability distribution
     * @param Q second probability distribution
     * @return Jensen-Shannon divergence
     */
    public static double ms_jensenshannon(Matrix P, Matrix Q) {
        if (P.length() != Q.length()) {
            throw new IllegalArgumentException("Distributions must have the same size");
        }

        double sum1 = 0.0;
        double sum2 = 0.0;

        for (int i = 0; i < P.length(); i++) {
            double pi = P.get(i);
            double qi = Q.get(i);
            double m = (pi + qi) / 2.0;

            if (pi > 0 && m > 0) {
                sum1 += pi * Math.log(pi / m);
            }

            if (qi > 0 && m > 0) {
                sum2 += qi * Math.log(qi / m);
            }
        }

        return 0.5 * (sum1 + sum2);
    }
}
