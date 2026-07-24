/**
 * @file Markovian Arrival Process cumulative distribution function
 *
 * Computes CDF values for MAP inter-arrival times using CTMC uniformization techniques.
 * Essential for probability analysis and performance evaluation of stochastic arrival processes.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.api.mc.Ctmc_foxglynn;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Map_cdf {
    private Map_cdf() {}

    /**
     * Computes the cumulative distribution function (CDF) of the inter-arrival times of a Markovian Arrival Process (MAP).
     *
     * The MAP is represented by two matrices: D0 and D1. D0 is the hidden transition matrix, representing
     * transitions without an observed event, while D1 is the visible transition matrix, representing transitions
     * with an observed event. The CDF values are calculated for a given set of points.
     *
     * @param D0     the hidden transition matrix of the MAP
     * @param D1     the visible transition matrix of the MAP
     * @param points a matrix containing the points at which to compute the CDF
     * @return a matrix containing the CDF values corresponding to the provided points
     */
    public static Matrix map_cdf(Matrix D0, Matrix D1, Matrix points) {
        Matrix CDFVals = new Matrix(1, points.length());
        Matrix pie = Map_pie.map_pie(D0, D1);
        Matrix e1 = Matrix.ones(D0.getNumRows(), 1);

        // Uniformization writes exp(D0*t) as a Poisson mixture of powers of
        // I + D0/lambda and is valid only when that matrix is nonnegative, i.e.
        // when D0 is a proper sub-generator. A matrix-exponential process (ME,
        // CME, RAP) has negative off-diagonal entries in D0 by construction, and
        // the mixture then loses all cancellation: a CME of order 11 returned a
        // "CDF" of 2.3e5. Those D0 are evaluated with a direct matrix
        // exponential instead, which is what MATLAB map_cdf does for every MAP.
        boolean isSubGenerator = true;
        for (int i = 0; i < D0.getNumRows() && isSubGenerator; i++) {
            for (int j = 0; j < D0.getNumCols(); j++) {
                if (i != j && D0.get(i, j) < 0) {
                    isSubGenerator = false;
                    break;
                }
            }
        }

        double nanVal = 0.0;
        for (int t = 0; t < points.length(); t++) {
            Matrix output = isSubGenerator
                    ? Ctmc_foxglynn.ctmc_foxglynn(pie, D0, points.get(t)).mult(e1)
                    : pie.mult(D0.scale(points.get(t)).expm()).mult(e1);
            double val = 1 - output.value();
            if (Double.isNaN(val)) {
                val = nanVal;
            } else {
                // after it finds the first non-zero, set nanVal to 1.0
                nanVal = 1.0;
            }
            CDFVals.set(0, t, val);
        }
        return CDFVals;
    }

    /**
     * CDF of MAP inter-arrival times when MAP is a MatrixCell.
     */
    public static Matrix map_cdf(MatrixCell MAP, Matrix points) {
        return map_cdf(MAP.get(0), MAP.get(1), points);
    }
}
