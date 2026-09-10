/**
 * @file Markovian Arrival Process feasible block matrix construction
 *
 * Constructs feasible MAP representations when exact moment matching fails by adjusting
 * parameters and ensuring mathematical constraints are satisfied. Used for robust MAP fitting.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import static jline.io.InputOutput.line_warning;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Map_feasblock {
    private Map_feasblock() {}

    /**
     * Fits the most similar feasible MAP when exact moment matching fails.
     * Ensures feasibility constraints are met.
     *
     * @param E1 First moment (mean)
     * @param E2 Second moment
     * @param E3 Third moment
     * @param G2 Autocorrelation decay ratio rho(i)/rho(i-1)
     * @param OPT Optional parameter ('scv' means E2 is squared coefficient of variation)
     * @return MAP as Matrix[] where result[0] = D0 and result[1] = D1
     */
    public static Matrix[] map_feasblock(double E1, double E2, double E3, double G2, String OPT) {
        double actualE2 = E2;
        double actualE3 = E3;

        // Handle exponential case
        if (actualE2 == 2 * E1 * E1) {
            Matrix D0 = new Matrix(new double[][] {{-1.0, 0.0}, {0.0, -1.0}});
            Matrix D1 = new Matrix(new double[][] {{0.5, 0.5}, {0.5, 0.5}});
            MatrixCell scaledMAP = Map_scale.map_scale(D0, D1, E1);
            return new Matrix[] {scaledMAP.get(0), scaledMAP.get(1)};
        }

        // Handle OPT parameter
        if (OPT != null && OPT.equalsIgnoreCase("scv")) {
            actualE2 = (1 + E2) * E1 * E1;
        }

        double tolerance = 1e-10;  // kpcfit_tol equivalent

        // Check feasibility constraints and adjust if necessary
        if (actualE2 <= 2 * E1 * E1) {
            line_warning("map_feasblock", "E2 failure (SCV<=1), setting SCV=1.001");
            actualE2 = (2 + tolerance) * E1 * E1;
        }

        double minE3 = (3.0 / 2.0) * actualE2 * actualE2 / E1;
        if (actualE3 <= minE3) {
            line_warning("map_feasblock", "E3 failure, setting E3=(3/2+1e-6)*E2^2/E1");
            actualE3 = (3.0 / 2.0 + tolerance) * actualE2 * actualE2 / E1;
        }

        // Call map_block with adjusted parameters
        return Map_block.map_block(E1, actualE2, actualE3, G2);
    }

    public static Matrix[] map_feasblock(double E1, double E2, double E3, double G2) {
        return map_feasblock(E1, E2, E3, G2, null);
    }
}
