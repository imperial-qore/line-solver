/**
 * @file Sequential convolution of APH distributions
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.Pair;
import jline.util.matrix.Matrix;

import java.util.List;

public final class Aph_convseq {
    private Aph_convseq() {}

    /**
     * Performs sequential convolution of multiple APH distributions.
     *
     * @param aphParams List of (alpha, T) pairs representing APH distributions to convolve
     * @return A Pair of (alpha, T) representing the convolved APH distribution
     * @throws IllegalArgumentException if the input list is empty
     */
    public static Pair<Matrix, Matrix> aph_convseq(List<Pair<Matrix, Matrix>> aphParams) {
        if (aphParams.isEmpty()) {
            throw new IllegalArgumentException("aphParams list cannot be empty");
        }

        if (aphParams.size() == 1) {
            return aphParams.get(0);
        }

        // Start with first two APHs convolved using sequential pattern
        Pair<Matrix, Matrix> first = Aph_simplify.aph_simplify(
                aphParams.get(0).getLeft(), aphParams.get(0).getRight(),
                aphParams.get(1).getLeft(), aphParams.get(1).getRight(),
                1.0, 1.0, 1);
        Matrix alpha = first.getLeft();
        Matrix T = first.getRight();

        // Convolve remaining APHs one by one
        for (int i = 2; i < aphParams.size(); i++) {
            Pair<Matrix, Matrix> result = Aph_simplify.aph_simplify(
                    alpha, T,
                    aphParams.get(i).getLeft(), aphParams.get(i).getRight(),
                    1.0, 1.0, 1);
            alpha = result.getLeft();
            T = result.getRight();
        }

        return new Pair<Matrix, Matrix>(alpha, T);
    }
}
