package jline.api.moment;

/**
 * Joint central moments from a joint survival array.
 *
 * <p>Composes the four edges that separate the two vertices, tail -&gt;
 * binomial -&gt; factorial -&gt; raw -&gt; central, reading the means off the raw
 * array. This is the whole path from a solver that produces survival
 * probabilities (a closed queueing network through its normalizing constants, a
 * CTMC through its stationary distribution, a simulator through a histogram) to
 * the covariances and the higher central moments.
 *
 * <p>Reference:
 * A. Heindl and A. van de Liefvoort. Moment conversions for discrete
 * distributions. PMCCS, 2003, Section 4.
 *
 * @since LINE 3.0
 */
public final class Moment_joint_central_from_tail {
    private Moment_joint_central_from_tail() {}

    /**
     * Converts a joint survival array into joint central moments.
     *
     * @param t flattened joint survival array in row-major order, covering the
     *          support, element 0 being 1
     * @param dims extents of the array, dims[j] = n_j+1, every dims[j] &gt;= 2
     * @return flattened joint central moments, same layout; the entry of
     *         multi-order e_j+e_l is the covariance of N_j and N_l
     */
    public static double[] moment_joint_central_from_tail(double[] t, int[] dims) {
        double[] b = Moment_joint_binomial_from_tail.moment_joint_binomial_from_tail(t, dims);
        double[] f = Moment_joint_factorial_from_binomial
                .moment_joint_factorial_from_binomial(b, dims);
        double[] m = Moment_joint_raw_from_factorial.moment_joint_raw_from_factorial(f, dims);
        return Moment_joint_central_from_raw.moment_joint_central_from_raw(m, dims);
    }
}
