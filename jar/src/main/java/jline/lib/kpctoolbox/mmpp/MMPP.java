/**
 * @file Markov Modulated Poisson Process (MMPP) functions
 *
 * Ported from MATLAB: matlab/lib/kpctoolbox/mmpp/
 *
 * @since LINE 3.0
 */
package jline.lib.kpctoolbox.mmpp;

import org.apache.commons.math3.util.FastMath;

import jline.api.mam.*;
import jline.api.mam.Mmpp2_fitc;
import jline.api.mam.Mmpp2_fitc_approx;
import jline.api.mam.Mmpp2_fit;
import jline.api.mam.Mmpp2_fit1;
import jline.api.mam.Mmpp_rand;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class MMPP {
    private MMPP() {}

    public static MatrixCell mmpp2_fit3(double E1, double E2, double E3, double G2) {
        return Mmpp2_fit.mmpp2_fit(E1, E2, E3, G2);
    }

    public static MatrixCell mmpp2_fit1(double mean, double scv, double skew, double idc) {
        return Mmpp2_fit1.mmpp2_fit1(mean, scv, skew, idc);
    }

    public static MatrixCell mmpp2_fit2(double mean, double scv, double skew, double g2) {
        if (scv == 1.0) {
            return Map_exponential.map_exponential(mean);
        }
        double E1 = mean;
        double E2 = (1 + scv) * E1 * E1;
        double E3 = -(2 * E1 * E1 * E1 - 3 * E1 * E2 - skew * FastMath.pow(E2 - E1 * E1, 1.5));
        return mmpp2_fit3(E1, E2, E3, g2);
    }

    public static MatrixCell mmpp2_fit4(double mean, double scv, double skew, double acf1) {
        double E1 = mean;
        double E2 = (1 + scv) * E1 * E1;
        double E3;
        if (skew == -1.0) {
            E3 = -1.0;
        } else {
            E3 = -(2 * E1 * E1 * E1 - 3 * E1 * E2 - skew * FastMath.pow(E2 - E1 * E1, 1.5));
        }
        double rho0 = (1 - 1 / scv) / 2;
        double g2 = acf1 / rho0;
        return mmpp2_fit3(E1, E2, E3, g2);
    }

    public static MatrixCell mmpp2_fitc(double mu, double bt1, double bt2, double binf,
                                        double m3t2, double t1, double t2) {
        Matrix[] result = Mmpp2_fitc.mmpp2_fitc(mu, bt1, bt2, binf, m3t2, t1, t2);
        MatrixCell MAP = new MatrixCell(2);
        MAP.set(0, result[0]);
        MAP.set(1, result[1]);
        return MAP;
    }

    public static MatrixCell mmpp2_fitc_approx(double a, double bt1, double bt2, double binf,
                                               double m3t2, double t1, double t2) {
        Matrix[] result = Mmpp2_fitc_approx.mmpp2_fitc_approx(a, bt1, bt2, binf, m3t2, t1, t2);
        MatrixCell MAP = new MatrixCell(2);
        MAP.set(0, result[0]);
        MAP.set(1, result[1]);
        return MAP;
    }

    public static MatrixCell mmpp2_fitc_theoretical(MatrixCell MAP) {
        return mmpp2_fitc_theoretical(MAP, 1.0, 10.0, 1e8);
    }

    public static MatrixCell mmpp2_fitc_theoretical(MatrixCell MAP, double t1, double t2, double tinf) {
        double a = Map_count_mean.map_count_mean(MAP, t1) / t1;
        double bt1 = Map_count_var.map_count_var(MAP, t1) / (a * t1);
        double bt2 = Map_count_var.map_count_var(MAP, t2) / (a * t2);
        double binf = Map_count_var.map_count_var(MAP, tinf) / (a * tinf);
        double[] mt2 = Map_count_moment.map_count_moment(MAP, t2, new int[]{1, 2, 3});
        double m3t2 = mt2[2] - 3 * mt2[1] * mt2[0] + 2 * mt2[0] * mt2[0] * mt2[0];

        return mmpp2_fitc(a, bt1, bt2, binf, m3t2, t1, t2);
    }

    public static MatrixCell mmpp_rand(int K) {
        return Mmpp_rand.mmpp_rand(K);
    }

    public static MatrixCell mmpp_rand() {
        return mmpp_rand(2);
    }
}
