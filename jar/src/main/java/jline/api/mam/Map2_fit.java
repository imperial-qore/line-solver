/**
 * @file Markovian Arrival Process two-phase fitting algorithms
 *
 * Fits MAP(2) processes to match specified moments and autocorrelation decay rates.
 * Fundamental algorithm for constructing realistic arrival process models from statistical data.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import java.util.Random;

import org.apache.commons.math3.util.FastMath;

import jline.io.Ret;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Map2_fit {
    private Map2_fit() {}

    /**
     * Fits a 2-phase Markovian Arrival Process (MAP2) to match the given moments
     * and decay rate of the autocorrelation function.
     */
    public static Ret.mamMAPFitReturn map2_fit(double e1, double e2, double e3, double g2) {
        Ret.mamMAPFitReturn result = new Ret.mamMAPFitReturn();
        double r1 = e1;
        double r2 = e2 / 2;
        double h2 = (r2 - FastMath.pow(r1, 2)) / FastMath.pow(r1, 2);
        if (e3 == -1.0) {
            double scv = (e2 - FastMath.pow(e1, 2)) / FastMath.pow(e1, 2);
            if (1 <= scv && scv < 3) {
                if (g2 < 0) {
                    double h3 = h2 - FastMath.pow(h2, 2);
                    e3 = 12 * FastMath.pow(e1, 3) * h2 + 6 * FastMath.pow(e1, 3) * h3
                            + 6 * FastMath.pow(e1, 3) * (1 + FastMath.pow(h2, 2));
                } else {
                    e3 = 1.501 * FastMath.pow(e2, 2) / e1;
                }
            } else if (3 <= scv) {
                e3 = 1.501 * FastMath.pow(e2, 2) / e1;
            } else if (0 < scv && scv < 1) {
                e3 = (1 + 1e-10) * (12 * FastMath.pow(e1, 3) * h2
                        + 6 * FastMath.pow(e1, 3) * (h2 * (1 - h2 - 2 * FastMath.sqrt(-h2)))
                        + 6 * FastMath.pow(e1, 3) * (1 + FastMath.pow(h2, 2)));
            }
        }

        if (e3 == -2.0) {
            double scv = (e2 - FastMath.pow(e1, 2)) / FastMath.pow(e1, 2);
            if (scv >= 1) {
                e3 = (3.0 / 2 + 1e-6) * FastMath.pow(e2, 2) / e1;
            } else if (0 < scv && scv < 1) {
                double h3 = h2 * (1 - h2 - 2 * FastMath.sqrt(-h2));
                e3 = 6 * FastMath.pow(e1, 3) * (FastMath.pow(h2, 2) + h3);
            }
        }

        if (e3 == -3.0) {
            double scv = (e2 - FastMath.pow(e1, 2)) / FastMath.pow(e1, 2);
            if (scv >= 1) {
                e3 = FastMath.pow(10.0, 6);
            } else if (0 < scv && scv < 1) {
                double h3 = FastMath.pow(-h2, 2);
                e3 = 6 * FastMath.pow(e1, 3) * (FastMath.pow(h2, 2) + h3);
            }
        }

        if (e3 == -4.0) {
            double scv = (e2 - FastMath.pow(e1, 2)) / FastMath.pow(e1, 2);
            double r = new Random().nextDouble();
            if (scv >= 1) {
                e3 = r * (3.0 / 2 + 1e-6) * FastMath.pow(e2, 2) / e1 + (1 - r) * FastMath.pow(10.0, 6);
            } else if (0 < scv && scv < 1) {
                double h3 = r * FastMath.pow(-h2, 2) + (1 - r) * h2 * (1 - h2 - 2 * FastMath.sqrt(-h2));
                e3 = 6 * FastMath.pow(e1, 3) * (FastMath.pow(h2, 2) + h3);
            }
        }

        if (e3 > -1 && e3 < 0) {
            double scv = (e2 - FastMath.pow(e1, 2)) / FastMath.pow(e1, 2);
            double r = FastMath.abs(e3);
            if (scv >= 1) {
                e3 = r * (3.0 / 2 + 1e-6) * FastMath.pow(e2, 2) / e1 + (1 - r) * FastMath.pow(10.0, 6);
            } else if (0 < scv && scv < 1) {
                double h3 = r * h2 * (1 - h2 - 2 * FastMath.sqrt(-h2)) + (1 - r) * FastMath.pow(-h2, 2);
                e3 = 6 * FastMath.pow(e1, 3) * (FastMath.pow(h2, 2) + h3);
            }
        }

        double r3 = e3 / 6;
        double h3 = (r3 * r1 - FastMath.pow(r2, 2)) / FastMath.pow(r1, 4);
        // MATLAB: b = h3 + h2^2 - h2; c = sqrt(b^2 + 4*h2^3) - both in h2, not r1
        double b = h3 + FastMath.pow(h2, 2) - h2;
        double c = FastMath.sqrt(FastMath.pow(b, 2) + 4 * FastMath.pow(h2, 3));

        if (r1 <= 0) {
            result.error = 10.0;
            return result;
        }

        if (h2 == 0.0) {
            if (h3 == 0.0 && g2 == 0.0) {
                result.MAP = Map_exponential.map_exponential(e1);
                result.error = 0.0;
                return result;
            } else {
                result.error = 20.0;
                return result;
            }
        }

        if (h2 > 0 && h3 > 0) {
            if (b >= 0) {
                if ((b - c) / (b + c) <= g2) {
                    MatrixCell MAP = new MatrixCell();
                    Matrix D0 = new Matrix(2, 2, 4);
                    D0.set(0, 0, -(2 * h2 + b - c));
                    D0.set(0, 1, 0);
                    D0.set(1, 0, 0);
                    D0.set(1, 1, -(2 * h2 + b + c));
                    D0.scaleEq(1.0 / (2 * r1 * h3));
                    MAP.set(0, D0);

                    Matrix D1 = new Matrix(2, 2, 4);
                    D1.set(0, 0, (2 * h2 + b - c) * (1 - b / c + g2 * (1 + b / c)));
                    D1.set(0, 1, (2 * h2 + b - c) * (1 + b / c) * (1 - g2));
                    D1.set(1, 0, (2 * h2 + b + c) * (1 - b / c) * (1 - g2));
                    D1.set(1, 1, (2 * h2 + b + c) * (1 + b / c + g2 * (1 - b / c)));
                    D1.scaleEq(1.0 / (4 * r1 * h3));
                    MAP.set(1, D1);

                    result.MAP = MAP;
                } else {
                    result.error = 51.0;
                }
            } else if (b < 0) {
                if (0 <= g2 && g2 < 1) {
                    MatrixCell MAP = new MatrixCell();
                    Matrix D0 = new Matrix(2, 2, 4);
                    D0.set(0, 0, -(2 * h2 + b - c));
                    D0.set(0, 1, 0);
                    D0.set(1, 0, 0);
                    D0.set(1, 1, -(2 * h2 + b + c));
                    D0.scaleEq(1.0 / (2 * r1 * h3));
                    MAP.set(0, D0);

                    Matrix D1 = new Matrix(2, 2, 4);
                    D1.set(0, 0, (2 * h2 + b - c) * (1 - b / c + g2 * (1 + b / c)));
                    D1.set(0, 1, (2 * h2 + b - c) * (1 + b / c) * (1 - g2));
                    D1.set(1, 0, (2 * h2 + b + c) * (1 - b / c) * (1 - g2));
                    D1.set(1, 1, (2 * h2 + b + c) * (1 + b / c + g2 * (1 - b / c)));
                    D1.scaleEq(1.0 / (4 * r1 * h3));
                    MAP.set(1, D1);

                    result.MAP = MAP;
                } else if (-(h3 + FastMath.pow(h2, 2)) / h2 <= g2 && g2 < 0) {
                    double a = (h3 + FastMath.pow(h2, 2)) / h2;
                    double d1 = ((1 - a) * (2 * h2 * g2 + b - c) + g2 * (b + c) - (b - c))
                            / ((1 - a) * (2 * h2 + b - c) + 2 * c);
                    double d2 = ((g2 - 1) * (b - c)) / ((1 - a) * (2 * h2 + b - c) + 2 * c);
                    MatrixCell MAP = new MatrixCell();
                    Matrix D0 = new Matrix(2, 2, 4);
                    D0.set(0, 0, -(2 * h2 + b - c));
                    D0.set(0, 1, (2 * h2 + b - c) * (1 - a));
                    D0.set(1, 0, 0);
                    D0.set(1, 1, -(2 * h2 + b + c));
                    D0.scaleEq(1.0 / (2 * r1 * h3));
                    MAP.set(0, D0);

                    Matrix D1 = new Matrix(2, 2, 4);
                    D1.set(0, 0, (2 * h2 + b - c) * d1);
                    D1.set(0, 1, (2 * h2 + b - c) * (a - d1));
                    D1.set(1, 0, (2 * h2 + b + c) * d2);
                    D1.set(1, 1, (2 * h2 + b + c) * (1 - d2));
                    D1.scaleEq(1.0 / (2 * r1 * h3));
                    MAP.set(1, D1);

                    result.MAP = MAP;
                } else {
                    result.error = 52.0;
                }
            }
            if (result.MAP.size() > 0 && !Map_isfeasible.map_isfeasible(result.MAP)) {
                result.error = -1.0;
            }
            return result;
        } else if (-0.25 <= h2 && h2 < 0
                && h2 * (1 - h2 - 2 * FastMath.sqrt(-h2)) <= h3
                && h3 <= -FastMath.pow(h2, 2)) {
            if (g2 >= 0) {
                if (g2 <= -FastMath.pow(h2 + FastMath.sqrt(-h3), 2) / h2) {
                    double a = (2 * h2 + b - c) * (h2 + FastMath.sqrt(-h3)) / (2 * h2 * FastMath.sqrt(-h3));
                    c = -c;
                    double d1 = ((1 - a) * (2 * h2 * g2 + b - c) + g2 * (b + c) - (b - c))
                            / ((1 - a) * (2 * h2 + b - c) + 2 * c);
                    double d2 = ((g2 - 1) * (b - c)) / ((1 - a) * (2 * h2 + b - c) + 2 * c);
                    MatrixCell MAP = new MatrixCell();
                    Matrix D0 = new Matrix(2, 2, 4);
                    D0.set(0, 0, -(2 * h2 + b - c));
                    D0.set(0, 1, (2 * h2 + b - c) * (1 - a));
                    D0.set(1, 0, 0);
                    D0.set(1, 1, -(2 * h2 + b + c));
                    D0.scaleEq(1.0 / (2 * r1 * h3));
                    MAP.set(0, D0);

                    Matrix D1 = new Matrix(2, 2, 4);
                    D1.set(0, 0, (2 * h2 + b - c) * d1);
                    D1.set(0, 1, (2 * h2 + b - c) * (a - d1));
                    D1.set(1, 0, (2 * h2 + b + c) * d2);
                    D1.set(1, 1, (2 * h2 + b + c) * (1 - d2));
                    D1.scaleEq(1.0 / (2 * r1 * h3));
                    MAP.set(1, D1);

                    result.MAP = MAP;
                } else {
                    result.error = 53.0;
                }
            } else if (g2 < 0) {
                if (g2 >= -(h3 + FastMath.pow(h2, 2)) / h2) {
                    double a = (h3 + FastMath.pow(h2, 2)) / h2;
                    c = -c;
                    double d1 = ((1 - a) * (2 * h2 * g2 + b - c) + g2 * (b + c) - (b - c))
                            / ((1 - a) * (2 * h2 + b - c) + 2 * c);
                    double d2 = ((g2 - 1) * (b - c)) / ((1 - a) * (2 * h2 + b - c) + 2 * c);
                    MatrixCell MAP = new MatrixCell();
                    Matrix D0 = new Matrix(2, 2, 4);
                    D0.set(0, 0, -(2 * h2 + b - c));
                    D0.set(0, 1, (2 * h2 + b - c) * (1 - a));
                    D0.set(1, 0, 0);
                    D0.set(1, 1, -(2 * h2 + b + c));
                    D0.scaleEq(1.0 / (2 * r1 * h3));
                    MAP.set(0, D0);

                    Matrix D1 = new Matrix(2, 2, 4);
                    D1.set(0, 0, (2 * h2 + b - c) * d1);
                    D1.set(0, 1, (2 * h2 + b - c) * (a - d1));
                    D1.set(1, 0, (2 * h2 + b + c) * d2);
                    D1.set(1, 1, (2 * h2 + b + c) * (1 - d2));
                    D1.scaleEq(1.0 / (2 * r1 * h3));
                    MAP.set(1, D1);

                    result.MAP = MAP;
                } else {
                    result.error = 54.0;
                }
            }
            if (result.MAP.size() > 0 && !Map_isfeasible.map_isfeasible(result.MAP)) {
                result.error = -1.0;
            }
            return result;
        } else {
            if (!(-0.25 <= h2 && h2 < 0
                    && h2 * (1 - h2 - 2 * FastMath.sqrt(-h2)) <= h3
                    && h3 <= -FastMath.pow(h2, 2))) {
                result.error = 30.0;
            } else if ((h2 > 0 && h3 < 0) || h2 * (1 - h2 - 2 * FastMath.sqrt(-h2)) > h3
                    || h3 <= FastMath.pow(h2, 2)) {
                result.error = 40.0;
            } else {
                result.error = 60.0; // Infeasible moment set
            }
        }

        return result;
    }

    /**
     * Fits a 2-phase Markovian Arrival Process (MAP2) using the first three moments.
     */
    public static Ret.mamMAPFitReturn map2_fit(double e1, double e2, double e3) {
        return map2_fit(e1, e2, -1.0, e3);
    }
}
