package jline.api.mc;

import jline.util.matrix.Matrix;

/**
 * DTMC random transition matrix generation.
 *
 * <p>Generates random stochastic transition matrices for discrete-time Markov chains
 * with controllable properties. Essential for algorithm testing, benchmarking,
 * and creating synthetic DTMC models with specified characteristics.
 *
 * @since LINE 3.0
 */
public final class Dtmc_rand {
    private Dtmc_rand() {}

    /**
     * Form a random infinitesimal generator of a DTMC
     *
     * @param length size of random matrix
     * @return Infinitesimal generator of CTMC
     */
    public static Matrix dtmc_rand(int length) {
        Matrix rand_matrix = new Matrix(length, length);
        rand_matrix.randMatrix(length);
        return Dtmc_makestochastic.dtmc_makestochastic(rand_matrix);
    }
}
