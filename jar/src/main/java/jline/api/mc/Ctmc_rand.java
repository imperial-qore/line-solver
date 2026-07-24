package jline.api.mc;

import jline.util.matrix.Matrix;

/**
 * CTMC random generator matrix construction.
 *
 * <p>Generates random infinitesimal generator matrices for continuous-time Markov chains
 * with specified properties. Used for testing CTMC algorithms and creating benchmark
 * problems with controllable statistical characteristics.
 *
 * @since LINE 3.0
 */
public final class Ctmc_rand {
    private Ctmc_rand() {}

    /**
     * Form a random infinitesimal generator of a CTMC
     *
     * @param length size of random matrix
     * @return Infinitesimal generator of CTMC
     */
    public static Matrix ctmc_rand(int length) {
        Matrix rand_matrix = new Matrix(length, length);
        rand_matrix.randMatrix(length);
        return Ctmc_makeinfgen.ctmc_makeinfgen(rand_matrix);
    }
}
