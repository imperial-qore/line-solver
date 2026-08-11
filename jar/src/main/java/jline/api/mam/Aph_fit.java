/**
 * @file Absorbing Phase-type distribution general fitting algorithms
 *
 * Fits APH distributions to specified moments using optimization and approximation techniques.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.util.FastMath;

import jline.GlobalConstants;
import jline.util.Utils;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Aph_fit {
    private Aph_fit() {}

    /**
     * Fits an acyclic phase-type (APH) distribution to the given moments of a random variable.
     */
    public static MatrixCell aph_fit(double e1, double e2, double e3, int nmax) {
        MatrixCell APH = new MatrixCell();
        if (Utils.isInf(e2) || Utils.isInf(e3)) {
            return Map_exponential.map_exponential(e1);
        }

        // Exponential moment set (scv == 1 and matching third moment) is a
        // degeneracy for the general APH(2) case formulas (division by zero
        // on the n2 boundary); the exponential is itself a valid phase-type,
        // so return it directly (matches MATLAB aph_fit).
        double scv = (e2 - FastMath.pow(e1, 2)) / FastMath.pow(e1, 2);
        if (FastMath.abs(scv - 1.0) < GlobalConstants.FineTol
                && FastMath.abs(e3 - 6.0 * FastMath.pow(e1, 3)) < GlobalConstants.FineTol) {
            return Map_exponential.map_exponential(e1);
        }

        double n2 = e2 / FastMath.pow(e1, 2);
        double n3 = e3 / e1 / e2;

        boolean n2_feas = false;
        boolean n3_ubfeas = false;
        boolean n3_lbfeas = false;
        double n = 1.0;
        double un = 0.0;
        double un_1 = un;
        while ((!n2_feas || !n3_lbfeas || !n3_ubfeas) && n < nmax) {
            n++;
            double pn =
                    ((n + 1) * (n2 - 2) / (3 * n2 * (n - 1)))
                            * (-2 * FastMath.sqrt(n + 1) / FastMath.sqrt(4 * (n + 1) - 3 * n * n2) - 1);
            double an = (n2 - 2) / (pn * (1 - n2)
                    + FastMath.sqrt(FastMath.pow(pn, 2) + pn * n * (n2 - 2) / (n - 1)));
            double ln =
                    ((3 + an) * (n - 1) + 2 * an) / ((n - 1) * (1 + an * pn))
                            - (2 * an * (n + 1)) / (2 * (n - 1) + an * pn * (n * an + 2 * n - 2));

            un =
                    (1 / (FastMath.pow(n, 2) * n2)) * (2 * (n - 2) * (n * n2 - n - 1)
                            * FastMath.sqrt(1 + n * (n2 - 2) / (n - 1))
                            + (n + 2) * (3 * n * n2 - 2 * n - 2));

            if (n2 >= (n + 1) / n && n2 <= (n + 4) / (n + 1)) {
                n2_feas = true;
                if (n3 >= ln) {
                    n3_lbfeas = true;
                }
            } else if (n2 >= (n + 4) / (n + 1)) {
                n2_feas = true;
                if (n3 >= n2 * (n + 1) / n) {
                    n3_lbfeas = true;
                }
            }

            if (n2 >= (n + 1) / n && n2 <= n / (n - 1)) {
                n2_feas = true;
                if (n3 <= un) {
                    n3_ubfeas = true;
                }
            } else if (n2 >= n / (n - 1)) {
                n2_feas = true;
                if (n3 < GlobalConstants.Inf) {
                    n3_ubfeas = true;
                }
            }
        }

        if ((!n2_feas || !n3_lbfeas || !n3_ubfeas) || (n == (double) nmax)) {
            System.out.print("'cannot match moment set exactly'");
            n2 = (n + 1) / n;
            n3 = 2 * n2 - 1;
        }

        if (n2 <= n / (n - 1) || n3 <= 2 * n2 - 1) {
            double b =
                    2 * (4 - n * (3 * n2 - 4)) / (n2 * (4 + n - n * n3)
                            + FastMath.sqrt(n * n2)
                                    * FastMath.sqrt(12 * FastMath.pow(n2, 2) * (n + 1)
                                            + 16 * n3 * (n + 1)
                                            + n2 * (n * (n3 - 15) * (n3 + 1) - 8 * (n3 + 3))));
            double a = (b * n2 - 2) * (n - 1) * b / ((b - 1) * n);
            double p = (b - 1) / a;

            double lambda = 1.0;
            double mu = lambda * (n - 1) / a;
            int nInt = (int) n;
            Matrix alpha = new Matrix(1, nInt, nInt);
            alpha.set(0, 0, p);
            alpha.set(0, nInt - 1, 1 - p);
            Matrix T = new Matrix(nInt, nInt, (int) FastMath.pow(n, 2));
            {
                int i = 0;
                while (i < n - 1) {
                    T.set(i, i, -mu);
                    T.set(i, i + 1, mu);
                    i++;
                }
            }
            T.set(nInt - 1, nInt - 1, -lambda);
            Matrix neg_T = T.copy();
            neg_T.scaleEq(-1.0);
            // MATLAB: -T*ones(n,1)*alpha, i.e. a full COLUMN of ones so that
            // (n x n)(n x 1)(1 x n) yields the rank-one exit matrix
            Matrix one = new Matrix(nInt, 1, nInt);
            int i = 0;
            while (i < nInt) {
                one.set(i, 0, 1);
                i++;
            }
            APH = Map_scale.map_scale(
                    Map_normalize.map_normalize(T, neg_T.mult(one).mult(alpha)).get(0),
                    Map_normalize.map_normalize(T, neg_T.mult(one).mult(alpha)).get(1),
                    e1);
        } else if (n2 > n / (n - 1) && n3 > un_1) {
            double K1 = n - 1;
            double K2 = n - 2;
            double K3 = 3 * n2 - 2 * n3;
            double K4 = n3 - 3;
            double K5 = n - n2;
            double K6 = 1 + n2 - n3;
            double K7 = n + n2 - n * n2;
            double K8 = 3 + 3 * FastMath.pow(n2, 2) + n3 - 3 * n2 * n3;
            // K9..K22 follow MATLAB aph_fit case 2, which relies on complex
            // arithmetic: the inner square root of K9 is negative for many
            // feasible moment sets and the imaginary parts cancel in f. Real
            // doubles would produce NaN, so the chain is evaluated in C.
            Complex K9inner = new Complex(
                    -16 * FastMath.pow(K1, 2) * FastMath.pow(K7, 6)
                            + FastMath.pow(4 * K1 * FastMath.pow(K5, 3)
                                    + FastMath.pow(K1, 2) * K2 * FastMath.pow(K4, 2) * n * FastMath.pow(n2, 2)
                                    + 4 * K2 * n * n2 * (K4 * FastMath.pow(n, 2) - 3 * K6 * n2 + K8 * n), 2), 0.0)
                    .sqrt();
            Complex K9 = K9inner
                    .add(4 * FastMath.pow(K2, 2) * K3 * FastMath.pow(n, 2) * n2
                            + FastMath.pow(K1, 2) * K2 * FastMath.pow(K4, 2) * n * FastMath.pow(n2, 2)
                            + 4 * K1 * K5 * (FastMath.pow(K5, 2) - 3 * K2 * K6 * n * n2))
                    .multiply(108 * FastMath.pow(K1, 2));
            double K10 = FastMath.pow(K4, 2) / (4 * FastMath.pow(K3, 2)) - K5 / (K1 * K3 * n2);
            Complex K9cbrt = K9.pow(1.0 / 3.0);
            Complex K11 = K9cbrt.reciprocal()
                    .multiply(FastMath.pow(2.0, 1.0 / 3) * (3 * FastMath.pow(K5, 2) + K2 * (K3 + 2 * K4) * n * n2)
                            / (K3 * n2));
            Complex K12 = K9cbrt.divide(3 * FastMath.pow(2.0, 7.0 / 3) * FastMath.pow(K1, 2) * K3 * n2);
            Complex K13 = K11.add(K12).add(K10).sqrt();
            Complex K14 = K13.reciprocal()
                    .multiply((6 * K1 * K3 * K4 * K5 + 4 * K2 * FastMath.pow(K3, 2) * n
                            - FastMath.pow(K1, 2) * FastMath.pow(K4, 3) * n2)
                            / (4 * FastMath.pow(K1, 2) * FastMath.pow(K3, 3) * n2));
            double K15 = -K4 / (2 * K3);
            Complex K16 = new Complex(2 * K10, 0.0).subtract(K11).subtract(K12).subtract(K14).sqrt();
            Complex K17 = new Complex(2 * K10, 0.0).subtract(K11).subtract(K12).add(K14).sqrt();
            Complex K18 = new Complex(
                    81 * FastMath.pow(4 * FastMath.pow(K5, 3)
                            + 4 * K2 * K4 * K5 * n * n2
                            + K1 * K2 * FastMath.pow(K4, 2) * n * FastMath.pow(n2, 2), 2)
                            - 48 * FastMath.pow(3 * FastMath.pow(K5, 2) + 2 * K2 * K4 * n * n2, 3), 0.0)
                    .sqrt().negate()
                    .add(36 * FastMath.pow(K5, 3) + 36 * K2 * K4 * K5 * n * n2
                            + 9 * K1 * K2 * FastMath.pow(K4, 2) * n * FastMath.pow(n2, 2));
            Complex K18cbrt = K18.pow(1.0 / 3.0);
            Complex K19 = new Complex(-K5 / (K1 * K4 * n2), 0.0)
                    .subtract(K18cbrt.reciprocal()
                            .multiply(FastMath.pow(2.0, 2.0 / 3) * (3 * FastMath.pow(K5, 2) + 2 * K2 * K4 * n * n2)
                                    / (FastMath.pow(3.0, 1.0 / 3) * K1 * K4 * n2)))
                    .subtract(K18cbrt.divide(FastMath.pow(6.0, 2.0 / 3) * K1 * K4 * n2));
            double K20 =
                    6 * K1 * K3 * K4 * K5 + 4 * K2 * FastMath.pow(K3, 2) * n
                            - FastMath.pow(K1, 2) * FastMath.pow(K4, 3) * n2;
            Complex K21 = K11.add(K12).add(K5 / (2 * n * K1 * K3));
            Complex K22 = K21.multiply(K21).multiply(4.0)
                    .subtract(n * K2 / (n2 * FastMath.pow(K1, 2) * K3))
                    .sqrt()
                    .add(3 * FastMath.pow(K4, 2) / (4 * FastMath.pow(K3, 2)) - 3 * K5 / (K1 * K3 * n2))
                    .sqrt();
            Complex fC = Complex.ZERO;
            if (n3 > un_1 && n3 < 3 * n2 / 2) {
                fC = K13.add(K15).subtract(K17);
            } else if (n3 == 2 * n2 / 2) {
                fC = K19;
            } else if (n3 > 3 * n2 / 2 && K20 > 0) {
                fC = K13.negate().add(K15).add(K16);
            } else if (K20 == 0.0) {
                fC = K22.add(K15);
            } else if (K20 < 0) {
                fC = K13.add(K15).add(K17);
            }
            // Imaginary parts cancel for feasible moment sets (as in MATLAB)
            double f = fC.getReal();
            double a = 2 * (f - 1) * (n - 1) / ((n - 1) * (n2 * FastMath.pow(f, 2) - 2 * f + 2) - n);
            double p = (f - 1) * a;
            double lambda = 1.0;
            double mu = lambda * (n - 1) / a;
            int nInt = (int) n;
            Matrix alpha = new Matrix(1, nInt, nInt);
            alpha.set(0, 0, p);
            alpha.set(0, 1, 1 - p);
            Matrix T = new Matrix(nInt, nInt, (int) FastMath.pow(n, 2));
            {
                int i = 0;
                while (i < n - 1) {
                    T.set(i, i, -mu);
                    T.set(i, i + 1, mu);
                    i++;
                }
            }
            T.set(nInt - 1, nInt - 1, -mu);
            T.set(0, 0, -lambda);
            T.set(0, 1, lambda);
            Matrix neg_T = T.copy();
            neg_T.scaleEq(-1.0);
            // MATLAB: -T*ones(n,1)*alpha (column of ones, see case 1)
            Matrix one = new Matrix(nInt, 1, nInt);
            int i = 0;
            while (i < nInt) {
                one.set(i, 0, 1);
                i++;
            }
            APH = Map_scale.map_scale(
                    Map_normalize.map_normalize(T, neg_T.mult(one).mult(alpha)).get(0),
                    Map_normalize.map_normalize(T, neg_T.mult(one).mult(alpha)).get(1),
                    e1);
        } else {
            System.out.print("moment set cannot be matched with an APH distribution");
        }

        return APH;
    }

    /**
     * Fits APH using the default maximum order of 10.
     */
    public static MatrixCell aph_fit(double e1, double e2, double e3) {
        return aph_fit(e1, e2, e3, 10);
    }
}
