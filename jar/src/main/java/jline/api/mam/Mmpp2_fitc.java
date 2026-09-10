/**
 * @file Markov Modulated Poisson Process counting-based fitting
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.solvers.BisectionSolver;

import jline.util.matrix.Matrix;

public final class Mmpp2_fitc {
    private Mmpp2_fitc() {}

    /**
     * Fits a MMPP(2) according to [Heffes and Lucantoni, 1986].
     */
    public static Matrix[] mmpp2_fitc(final double mu, double bt1, double bt2,
                                      final double binf, double m3t2, final double t1, double t2) {
        // Degenerate case
        if (Math.abs(binf - 1) < 1e-8 && Math.abs(binf - bt1) < 1e-8) {
            return new Matrix[]{
                    new Matrix(new double[]{-mu}),
                    new Matrix(new double[]{mu})
            };
        }

        if (!(binf > bt1 && bt1 > 1)) {
            return new Matrix[]{
                    new Matrix(new double[]{-mu}),
                    new Matrix(new double[]{mu})
            };
        }

        final double c = (binf - 1) / (binf - bt1);
        final double z = -c * Math.exp(-c);

        BisectionSolver solver = new BisectionSolver(1e-12);
        UnivariateFunction function = new UnivariateFunction() {
            @Override
            public double value(double w) {
                return z - w * Math.exp(w);
            }
        };

        double w;
        try {
            w = solver.solve(10000, function, -1.0, 10.0, 1.0);
        } catch (Exception e) {
            if (z > -1 / Math.E) {
                w = -1 + Math.sqrt(2 * (1 + Math.E * z));
            } else {
                w = 1.0;
            }
        }

        double x = (w + c) / t1;

        double k1 = Math.pow(mu, 3) * Math.pow(t2, 3);
        double k2 = 3 * Math.pow(mu, 2) * (binf - 1) * Math.pow(t2, 2);
        double k3 = 3 * mu * (binf - 1) / x * t2;
        double k4 = 3 * mu / Math.pow(x, 2) * (binf - 1) * t2 * Math.exp(-x * t2);
        double k5 = 6 * mu / Math.pow(x, 3) * (binf - 1) * (1 - Math.exp(-x * t2));
        double g1t2 = m3t2 + 3 * mu * t2 * (mu * t2 - 1) * bt2 + mu * t2 * (mu * t2 - 1) * (mu * t2 - 2);
        double h = (g1t2 - k1 - k2 - k3 * (-mu) - k4 * mu * x) / ((k3 / x) + k4 - k5);

        double r1;
        double r2;
        double l1;
        double l2;

        if (Math.abs(h) < 1e-4) {
            r1 = x / 2;
            r2 = x / 2;
            l2 = mu - 0.5 * Math.sqrt(2 * (binf - 1) * mu * x);
            l1 = mu + 0.5 * Math.sqrt(2 * (binf - 1) * mu * x);
        } else {
            double y = (binf - 1) * mu * Math.pow(x, 3) / (2 * Math.pow(h, 2));
            double tempR1 = x / 2 * (1 + 1 / Math.sqrt(4 * y + 1));
            double tempR2 = x - tempR1;

            if (tempR1 < tempR2) {
                double tmp = tempR1;
                tempR1 = tempR2;
                tempR2 = tmp;
            }
            r1 = tempR1;
            r2 = tempR2;

            double w_val = h / (r1 - r2);
            double w_min = -mu / r1 * (r1 + r2);
            double w_max = mu / r2 * (r1 + r2);

            if (w_val < w_min || w_val > w_max) {
                // ignore third moment to achieve feasibility; reassign r1/r2 so the
                // returned D0 uses the feasible rates (matches MATLAB mmpp2_fitc.m)
                double z_val = (binf - 1) * Math.pow(x, 3) * mu;
                double u = x * z_val / (2 * Math.pow(mu, 2) * Math.pow(x, 2) + z_val);
                r1 = u + (x - u) / 2;
                r2 = x - r1;
                double delta = Math.sqrt(z_val / (2 * r1 * r2));
                l2 = mu - r2 / x * delta;
                l1 = l2 + delta;
            } else {
                l2 = mu - h / (r1 - r2) * (r2 / (r1 + r2));
                l1 = h / (r1 - r2) + l2;
            }
        }

        Matrix D0 = new Matrix(2, 2);
        D0.set(0, 0, -(r1 + l1));
        D0.set(0, 1, r1);
        D0.set(1, 0, r2);
        D0.set(1, 1, -(r2 + l2));

        Matrix D1 = new Matrix(2, 2);
        D1.set(0, 0, l1);
        D1.set(0, 1, 0.0);
        D1.set(1, 0, 0.0);
        D1.set(1, 1, l2);

        return new Matrix[]{D0, D1};
    }
}
