package jline.api.mam;

import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.SimpleBounds;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.MultivariateOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.BOBYQAOptimizer;
import org.apache.commons.math3.analysis.MultivariateFunction;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * Fits a second-order Marked MMPP using optimization.
 */
public final class Mmpp2_fit_count_approx {
    private Mmpp2_fit_count_approx() {}

    public static Matrix[] mmpp2_fit_count_approx(final double a, final double bt1, final double bt2,
                                                   final double binf, final double m3t2, final double t1, final double t2) {
        ObjectiveFunction objectiveFunction = new ObjectiveFunction(new MultivariateFunction() {
            @Override
            public double value(double[] params) {
                double l1 = params[0];
                double l2 = params[1];
                double r1 = params[2];
                double r2 = params[3];

                double xa = (l1 * r2 + l2 * r1) / (r1 + r2);
                double factor = a / xa;

                double xbt1 = (r1 * (2 * Math.pow(l1, 2) * Math.pow(r2, 2) * t1 * factor - 2 * Math.pow(l2, 2) * r2 - 2 * Math.pow(l1, 2) * r2 +
                        2 * Math.pow(l2, 2) * Math.pow(r2, 2) * t1 * factor + 4 * l1 * l2 * r2 +
                        2 * Math.pow(l1, 2) * r2 * Math.exp(-r1 * t1 * factor - r2 * t1 * factor) +
                        2 * Math.pow(l2, 2) * r2 * Math.exp(-r1 * t1 * factor - r2 * t1 * factor) -
                        4 * l1 * l2 * Math.pow(r2, 2) * t1 * factor -
                        4 * l1 * l2 * r2 * Math.exp(-r1 * t1 * factor - r2 * t1 * factor)) +
                        Math.pow(r1, 2) * (2 * r2 * t1 * factor * Math.pow(l1, 2) - 4 * r2 * t1 * factor * l1 * l2 +
                                2 * r2 * t1 * factor * Math.pow(l2, 2))) /
                        (t1 * factor * Math.pow(r1 + r2, 3) * (l1 * r2 + l2 * r1)) + 1;

                double xbt2;
                if (Math.abs(t1 - t2) > 1e-10) {
                    xbt2 = (r1 * (2 * Math.pow(l1, 2) * Math.pow(r2, 2) * t2 * factor - 2 * Math.pow(l2, 2) * r2 - 2 * Math.pow(l1, 2) * r2 +
                            2 * Math.pow(l2, 2) * Math.pow(r2, 2) * t2 * factor + 4 * l1 * l2 * r2 +
                            2 * Math.pow(l1, 2) * r2 * Math.exp(-r1 * t2 * factor - r2 * t2 * factor) +
                            2 * Math.pow(l2, 2) * r2 * Math.exp(-r1 * t2 * factor - r2 * t2 * factor) -
                            4 * l1 * l2 * Math.pow(r2, 2) * t2 * factor -
                            4 * l1 * l2 * r2 * Math.exp(-r1 * t2 * factor - r2 * t2 * factor)) +
                            Math.pow(r1, 2) * (2 * r2 * t2 * factor * Math.pow(l1, 2) - 4 * r2 * t2 * factor * l1 * l2 +
                                    2 * r2 * t2 * factor * Math.pow(l2, 2))) /
                            (t2 * factor * Math.pow(r1 + r2, 3) * (l1 * r2 + l2 * r1)) + 1;
                } else {
                    xbt2 = xbt1;
                }

                double xbinf = ((2 * r2 * Math.pow(l1, 2) - 4 * r2 * l1 * l2 + 2 * r2 * Math.pow(l2, 2)) * Math.pow(r1, 2) +
                        (2 * Math.pow(l1, 2) * Math.pow(r2, 2) - 4 * l1 * l2 * Math.pow(r2, 2) + 2 * Math.pow(l2, 2) * Math.pow(r2, 2)) * r1) /
                        (Math.pow(r1 + r2, 3) * (l1 * r2 + l2 * r1)) + 1;

                double t = t2 * factor;
                double d = r1 + r2;
                double p = (l1 - l2) * (r1 - r2);
                double xg3t = Math.pow(xa, 3) * Math.pow(t, 3) +
                        3 * Math.pow(xa, 2) * (xbinf - 1) * Math.pow(t, 2) +
                        3 * xa * (xbinf - 1) / d * (p / d - xa) * t +
                        3 * xa / Math.pow(d, 2) * (xbinf - 1) * (p + xa * d) * t * Math.exp(-t * d) -
                        6 * xa / Math.pow(d, 3) * (xbinf - 1) * p * (1 - Math.exp(-t * d));
                double xm3t2 = xg3t - 3 * xa * t * (xa * t - 1) * xbt2 - xa * t * (xa * t - 1) * (xa * t - 2);

                double obj = 0.0;
                obj += Math.pow(xa / a - 1, 2);
                obj += Math.pow(xbt1 / bt1 - 1, 2);
                if (Math.abs(t1 - t2) > 1e-10) {
                    obj += Math.pow(xbt2 / bt2 - 1, 2);
                }
                obj += Math.pow(xbinf / binf - 1, 2);
                obj += Math.pow(xm3t2 / m3t2 - 1, 2);

                return obj;
            }
        });

        double[] initialGuess = new double[]{
                a * 0.75,
                a * 1.5,
                1.0 / 3.0,
                2.0 / 3.0
        };

        double[] lowerBounds = new double[]{1e-6, 1e-6, 0.0, 0.0};
        double[] upperBounds = new double[]{Double.MAX_VALUE, Double.MAX_VALUE, Double.MAX_VALUE, Double.MAX_VALUE};

        MultivariateOptimizer optimizer = new BOBYQAOptimizer(9);

        PointValuePair optimum = optimizer.optimize(
                objectiveFunction,
                GoalType.MINIMIZE,
                new InitialGuess(initialGuess),
                new SimpleBounds(lowerBounds, upperBounds),
                new MaxEval(10000)
        );

        double[] xopt = optimum.getPoint();
        double l1 = xopt[0];
        double l2 = xopt[1];
        double r1 = xopt[2];
        double r2 = xopt[3];

        Matrix D0 = new Matrix(2, 2);
        D0.set(0, 0, -(l1 + r1));
        D0.set(0, 1, r1);
        D0.set(1, 0, r2);
        D0.set(1, 1, -(l2 + r2));

        Matrix D1 = new Matrix(2, 2);
        D1.set(0, 0, l1);
        D1.set(0, 1, 0.0);
        D1.set(1, 0, 0.0);
        D1.set(1, 1, l2);

        MatrixCell fit = Map_scale.map_scale(D0, D1, 1.0 / a);
        return new Matrix[]{fit.get(0), fit.get(1)};
    }
}
