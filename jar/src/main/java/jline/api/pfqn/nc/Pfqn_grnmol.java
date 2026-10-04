/**
 * @file Normalizing constant using Grundmann-Moeller quadrature
 *
 * Implements the Grundmann-Moeller simplex quadrature rule for computing normalizing
 * constants in product-form queueing networks.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.util.FastMath;

import jline.util.Maths;
import jline.util.matrix.Matrix;

public final class Pfqn_grnmol {
    private Pfqn_grnmol() {}

    /**
     * Compute the normalizing constant using Grundmann-Moeller quadrature
     */
    public static double pfqn_grnmol(Matrix L, Matrix N) {
        int M = L.getNumRows();
        int R = L.getNumCols();
        double G = 0.0;
        int S = (int) FastMath.ceil((N.elementSum() - 1) / 2);

        for (int i = 0; i <= S; i++) {
            double cVal = (double) (2 * (S - i) + M);
            double w = FastMath.pow(2.0, (double) (-2 * S))
                    * FastMath.pow(-1.0, (double) i)
                    * FastMath.pow(cVal, (double) (2 * S + 1))
                    / FastMath.exp(Maths.factln(i))
                    / FastMath.exp(Maths.factln(i + (int) cVal));

            List<int[]> compositions = generateCompositions(M, S - i);
            double Hi = 0.0;

            for (int[] comp : compositions) {
                double prodVal = 1.0;
                for (int r = 0; r < R; r++) {
                    double sumLr = 0.0;
                    for (int m = 0; m < M; m++) {
                        sumLr += (2.0 * comp[m] + 1.0) / cVal * L.get(m, r);
                    }
                    prodVal *= FastMath.pow(sumLr, N.get(r));
                }
                Hi += prodVal;
            }
            G += w * Hi;
        }

        double logCoeff = Maths.factln(N.elementSum() + M - 1) - Matrix.factln(N).elementSum();
        G = G * FastMath.exp(logCoeff);

        return G;
    }

    /**
     * Generate all compositions of integer k into m non-negative parts.
     */
    private static List<int[]> generateCompositions(int m, int k) {
        List<int[]> result = new ArrayList<int[]>();
        if (m == 1) {
            result.add(new int[]{k});
            return result;
        }
        if (k == 0) {
            result.add(new int[m]);
            return result;
        }
        for (int i = 0; i <= k; i++) {
            List<int[]> subComps = generateCompositions(m - 1, k - i);
            for (int[] sub : subComps) {
                int[] comp = new int[m];
                comp[0] = i;
                System.arraycopy(sub, 0, comp, 1, m - 1);
                result.add(comp);
            }
        }
        return result;
    }
}
