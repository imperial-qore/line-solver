/**
 * @file Markovian Arrival Process eigenvalue-based correlation analysis
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.MatrixCell;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.util.FastMath;

public final class Map_gamma2 {
    private Map_gamma2() {}

    /**
     * Returns the largest non-unit eigenvalue (both real and imaginary parts) of the
     * embedded Discrete-Time Markov Chain (DTMC) of a given Markovian Arrival Process (MAP).
     *
     * @param MAP The Markovian Arrival Process stored in a MatrixCell.
     * @return A double array of size 2: [real part, abs(imag part)] of the largest non-unit eigenvalue.
     */
    public static double[] map_gamma2(MatrixCell MAP) {
        org.apache.commons.math3.linear.RealMatrix realMatrix =
                MatrixUtils.createRealMatrix(Map_embedded.map_embedded(MAP).toArray2D());

        EigenDecomposition eigDecomp = new EigenDecomposition(realMatrix);
        double[] realPart = eigDecomp.getRealEigenvalues();
        double[] imgPart = eigDecomp.getImagEigenvalues();

        // Find the index of the maximum absolute value
        int maxIndex = realPart.length - 1;
        int secondMaxIndex = realPart.length - 1;
        double maxAbsValue = -Double.MIN_VALUE;
        for (int i = 0; i < realPart.length; i++) {
            double abs = FastMath.sqrt(realPart[i] * realPart[i] + imgPart[i] * imgPart[i]);
            if (abs > maxAbsValue) {
                maxAbsValue = abs;
                secondMaxIndex = maxIndex;
                maxIndex = i;
            }
        }

        return new double[]{realPart[secondMaxIndex], FastMath.abs(imgPart[secondMaxIndex])};
    }
}
