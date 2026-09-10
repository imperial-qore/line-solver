/**
 * @file Markovian Arrival Process kurtosis computation
 *
 * Computes kurtosis of MAP inter-arrival times measuring tail heaviness and distribution shape.
 * Important for characterizing extreme behavior and heavy-tailed properties in arrival processes.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.MatrixCell;

public final class Map_kurt {
    private Map_kurt() {}

    /**
     * Computes the kurtosis of the inter-arrival times in a Markovian Arrival Process (MAP).
     * <p>
     * The kurtosis is computed using the formula:
     * KURT = (m4 - 4*m3*m1 + 6*m2*m1^2 - 3*m1^4) / Var(X)^2
     * <p>
     * where mi is the i-th moment and Var(X) is the variance of the inter-arrival times.
     *
     * @param MAP The Markovian Arrival Process stored in a MatrixCell, containing the
     *            (D0, D1) matrices
     * @return The kurtosis of the inter-arrival times
     */
    public static double map_kurt(MatrixCell MAP) {
        // Compute the first four moments
        double m1 = Map_moment.map_moment(MAP, 1);
        double m2 = Map_moment.map_moment(MAP, 2);
        double m3 = Map_moment.map_moment(MAP, 3);
        double m4 = Map_moment.map_moment(MAP, 4);

        // Compute variance
        double variance = Map_var.map_var(MAP);

        // Compute kurtosis using the formula
        double numerator = m4 - 4.0 * m3 * m1 + 6.0 * m2 * m1 * m1 - 3.0 * m1 * m1 * m1 * m1;
        double denominator = variance * variance;

        return numerator / denominator;
    }
}
