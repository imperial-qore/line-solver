package jline.lib.kpctoolbox;

import jline.lib.kpctoolbox.mmpp.MMPP;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * Facade class for MMPP2 fitting functions.
 */
public final class MMPP2 {

    private MMPP2() {}

    public static MatrixCell fit(Matrix moments) {
        double[] m = moments.toArray1D();
        return MMPP.mmpp2_fit3(m[0], m[1], m[2], 0.5);
    }

    public static MatrixCell fit(Matrix moments, double acf1) {
        double[] m = moments.toArray1D();
        double scv = (m[1] - m[0] * m[0]) / (m[0] * m[0]);
        double rho0 = (1 - 1 / scv) / 2;
        double g2 = (rho0 != 0.0) ? acf1 / rho0 : 0.5;
        return MMPP.mmpp2_fit3(m[0], m[1], m[2], g2);
    }

    public static MatrixCell fit1(Matrix moments) {
        double[] m = moments.toArray1D();
        return MMPP.mmpp2_fit1(m[0], m[1], m[2], m[3]);
    }

    public static MatrixCell fit2(Matrix moments, double acf1) {
        double[] m = moments.toArray1D();
        return MMPP.mmpp2_fit4(m[0], m[1], m[2], acf1);
    }

    public static MatrixCell fit3(Matrix moments, double acf1, double acf2) {
        double[] m = moments.toArray1D();
        double g2 = (acf1 != 0.0) ? acf2 / acf1 : 0.5;
        return MMPP.mmpp2_fit2(m[0], m[1], m[2], g2);
    }
}
