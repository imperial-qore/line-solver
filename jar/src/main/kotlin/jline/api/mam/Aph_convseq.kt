/**
 * @file Sequential convolution of APH distributions
 *
 * Convolves a sequence of APH distributions by repeatedly applying aph_simplify
 * with sequential pattern.
 *
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix

/**
 * Performs sequential convolution of multiple APH distributions.
 *
 * This function takes a list of APH distributions (each represented by their
 * initial probability vector alpha and generator matrix T) and convolves them
 * sequentially using the aph_simplify function with pattern=1 (sequence).
 *
 * @param aphParams List of (alpha, T) pairs representing APH distributions to convolve
 * @return A Pair of (alpha, T) representing the convolved APH distribution
 * @throws IllegalArgumentException if the input list is empty
 */
fun aph_convseq(aphParams: List<Pair<Matrix, Matrix>>): Pair<Matrix, Matrix> {
    if (aphParams.isEmpty()) {
        throw IllegalArgumentException("aphParams list cannot be empty")
    }

    // If only one APH, return it unchanged
    if (aphParams.size == 1) {
        return aphParams[0]
    }

    // Start with first two APHs convolved using sequential pattern
    var (alpha, T) = aph_simplify(
        aphParams[0].first, aphParams[0].second,
        aphParams[1].first, aphParams[1].second,
        1.0, 1.0, 1  // pattern=1 for sequential
    )

    // Convolve remaining APHs one by one
    for (i in 2 until aphParams.size) {
        val result = aph_simplify(
            alpha, T,
            aphParams[i].first, aphParams[i].second,
            1.0, 1.0, 1  // pattern=1 for sequential
        )
        alpha = result.first
        T = result.second
    }

    return Pair(alpha, T)
}

/**
 * APH convolution sequence algorithms
 */
@Suppress("unused")
class AphConvseqAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}
