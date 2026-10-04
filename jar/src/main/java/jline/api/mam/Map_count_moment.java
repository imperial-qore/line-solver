/**
 * @file Markovian Arrival Process counting process moment analysis
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.Maths;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Map_count_moment {
    private Map_count_moment() {}

    /**
     * Computes power moments of counts at resolution t for a Markovian Arrival Process (MAP).
     */
    public static double map_count_moment(MatrixCell MAP, double t, int order) {
        int n = MAP.get(0).getNumRows();
        Matrix theta = Map_prob.map_prob(MAP);
        Matrix e = Matrix.ones(n, 1);

        return computeHigherOrderMoment(MAP, t, theta, e, order);
    }

    /**
     * Computes multiple power moments of counts at resolution t for a MAP.
     */
    public static double[] map_count_moment(MatrixCell MAP, double t, int[] orders) {
        double[] results = new double[orders.length];
        for (int i = 0; i < orders.length; i++) {
            results[i] = map_count_moment(MAP, t, orders[i]);
        }
        return results;
    }

    private static double mgfunc(double z, MatrixCell MAP, double t, Matrix theta, Matrix e) {
        Matrix D0 = MAP.get(0);
        Matrix D1 = MAP.get(1);

        double expZ = Math.exp(z);
        Matrix D1_scaled = D1.scale(expZ);
        Matrix generator = D0.add(1.0, D1_scaled);

        Matrix generatorScaled = generator.scale(t);
        Matrix expMatrix = Maths.matrixExp(generatorScaled);

        Matrix temp = theta.mult(expMatrix);
        return temp.mult(e).get(0, 0);
    }

    private static double computeHigherOrderMoment(MatrixCell MAP, double t, Matrix theta, Matrix e, int order) {
        // Order-th derivative of the count MGF at z=0 via a consistent 5-point central
        // stencil (mirrors MATLAB map_count_moment's numerical branch, derivest orders 1..4).
        // The 5-point stencils are O(h^4) accurate for orders 1-2 and O(h^2) for orders 3-4;
        // the step is chosen per order to balance truncation against floating-point rounding.
        final int points = 5;
        double[] weights = computeFiniteDifferenceWeights(order);
        double h = (order <= 2) ? 1e-3 : 5e-3;
        int offset = points / 2; // 2

        double result = 0.0;
        for (int i = 0; i < points; i++) {
            double z = (i - offset) * h;
            double mgf_val = mgfunc(z, MAP, t, theta, e);
            result += weights[i] * mgf_val;
        }

        return result / Math.pow(h, (double) order);
    }

    /**
     * Central-difference stencil weights on the 5-point grid {-2h,-h,0,h,2h} for the
     * derivative of the given order. Only orders 1-4 are supported numerically; higher
     * orders require symbolic differentiation (not available without a CAS).
     */
    private static double[] computeFiniteDifferenceWeights(int order) {
        switch (order) {
            case 1: return new double[]{1.0 / 12.0, -2.0 / 3.0, 0.0, 2.0 / 3.0, -1.0 / 12.0};
            case 2: return new double[]{-1.0 / 12.0, 4.0 / 3.0, -5.0 / 2.0, 4.0 / 3.0, -1.0 / 12.0};
            case 3: return new double[]{-0.5, 1.0, 0.0, -1.0, 0.5};
            case 4: return new double[]{1.0, -4.0, 6.0, -4.0, 1.0};
            default:
                throw new UnsupportedOperationException(
                        "map_count_moment supports numerical moments of order 1-4; got order " + order);
        }
    }
}
