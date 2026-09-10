/**
 * @file Markovian Arrival Process skewness computation
 *
 * Computes skewness of MAP inter-arrival times measuring asymmetry in distributions.
 * Important for statistical characterization and shape analysis of arrival patterns.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import org.apache.commons.math3.util.FastMath;

import java.util.HashMap;
import java.util.Map;

public final class Map_skew {
    private Map_skew() {}

    /**
     * Computes the skewness of the inter-arrival times for a MAP.
     *
     * <p>The skewness measures the asymmetry of the distribution of inter-arrival times. It is
     * calculated using the third central moment normalized by the cube of the standard deviation.
     *
     * @param D0 the hidden transition matrix of the MAP
     * @param D1 the visible transition matrix of the MAP
     * @return the skewness of the inter-arrival times
     */
    public static double map_skew(Matrix D0, Matrix D1) {
        Map<Integer, Double> m = new HashMap<Integer, Double>();
        for (int i = 1; i <= 3; i++) {
            m.put(i, Map_moment.map_moment(D0, D1, i));
        }
        double M3 = m.get(3) - 3 * m.get(2) * m.get(1) + 2 * FastMath.pow(m.get(1), 3);
        return M3 / FastMath.pow(Math.sqrt(Map_scv.map_scv(D0, D1)) * m.get(1), 3);
    }

    /**
     * Computes the skewness of the inter-arrival times for a MAP using a MatrixCell.
     *
     * <p>This method is a convenience overload that extracts the D0 and D1 matrices from the provided
     * MatrixCell and calculates the skewness of the inter-arrival times.
     *
     * @param MAP the MatrixCell representing the MAP, containing the D0 and D1 matrices
     * @return the skewness of the inter-arrival times
     */
    public static double map_skew(MatrixCell MAP) {
        return map_skew(MAP.get(0), MAP.get(1));
    }
}
