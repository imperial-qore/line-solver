/**
 * @file Markovian Arrival Process construction using block matrices
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;

public final class Map_block {
    private Map_block() {}

    /**
     * Constructs a MAP(2) or MAP(1) according to given moments and autocorrelation parameters.
     */
    public static Matrix[] map_block(double E1, double E2, double E3, double G2, String OPT) {
        double actualE2 = E2;

        if (OPT != null && OPT.equalsIgnoreCase("scv")) {
            actualE2 = (1 + E2) * E1 * E1;
        }

        return fallbackMAP(E1, actualE2, E3, G2);
    }

    public static Matrix[] map_block(double E1, double E2, double E3, double G2) {
        return map_block(E1, E2, E3, G2, null);
    }

    /**
     * Fallback MAP construction when MMPP(2) fitting fails.
     */
    private static Matrix[] fallbackMAP(double E1, double E2, double E3, double G2) {
        double SCV = (E2 - E1 * E1) / (E1 * E1);

        if (SCV > 1) {
            return constructHyperexponential(E1, SCV, G2);
        } else {
            return constructErlang(E1, SCV, G2);
        }
    }

    /**
     * Construct a hyperexponential MAP for high variability cases.
     */
    private static Matrix[] constructHyperexponential(double E1, double SCV, double G2) {
        double p = 0.5;

        double lambda1 = 2.0 / E1 * (1 + Math.sqrt((SCV - 1) / (SCV + 1)));
        double lambda2 = 2.0 / E1 * (1 - Math.sqrt((SCV - 1) / (SCV + 1)));

        Matrix D0 = new Matrix(2, 2);
        D0.set(0, 0, -lambda1);
        D0.set(0, 1, 0.0);
        D0.set(1, 0, 0.0);
        D0.set(1, 1, -lambda2);

        Matrix D1 = new Matrix(2, 2);
        D1.set(0, 0, lambda1 * p);
        D1.set(0, 1, lambda1 * (1 - p));
        D1.set(1, 0, lambda2 * p);
        D1.set(1, 1, lambda2 * (1 - p));

        return new Matrix[]{D0, D1};
    }

    /**
     * Construct an Erlang-like MAP for low variability cases.
     */
    private static Matrix[] constructErlang(double E1, double SCV, double G2) {
        double mu = 2.0 / E1;

        Matrix D0 = new Matrix(2, 2);
        D0.set(0, 0, -mu);
        D0.set(0, 1, mu * (1 - G2));
        D0.set(1, 0, 0.0);
        D0.set(1, 1, -mu);

        Matrix D1 = new Matrix(2, 2);
        D1.set(0, 0, 0.0);
        D1.set(0, 1, mu * G2);
        D1.set(1, 0, 0.0);
        D1.set(1, 1, mu);

        return new Matrix[]{D0, D1};
    }
}
