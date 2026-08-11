/**
 * @file M3PP(2,2) multiclass fitting with covariance constraint optimization
 *
 * Implements constrained optimization for fitting M3PP(2,2) parameters given an underlying
 * MMPP(2), using covariance constraints between arrival classes.
 *
 * @since LINE 3.0
 */
package jline.api.mam.m3pp;

import jline.GlobalConstants;
import jline.util.matrix.Matrix;

public class M3pp22_fitc_approx_cov_multiclass {
    private M3pp22_fitc_approx_cov_multiclass() {}

    /**
     * Fits a M3PP(2,2) given the underlying MMPP(2).
     */
    public static Matrix[] m3pp22_fitc_approx_cov_multiclass(Matrix[] mmpp, double[] ai, double st3, double t3) {
        int m = ai.length;

        if (m > 2) {
            throw new IllegalArgumentException("No more than two classes supported");
        }

        Matrix[] fit = mmpp;

        // Degenerate case: Poisson process
        if (fit[0].getNumRows() == 1) {
            Matrix D0 = fit[0];
            Matrix D1 = fit[1];
            Matrix[] result = new Matrix[2 + m];
            result[0] = D0;
            result[1] = D1;
            double a = 0.0;
            for (double v : ai) a += v;
            double[] pi = new double[m];
            for (int i = 0; i < m; i++) pi[i] = ai[i] / a;
            for (int i = 0; i < m; i++) {
                result[2 + i] = D1.scale(pi[i]);
            }
            return result;
        }

        // Just a single class
        if (m == 1) {
            return new Matrix[]{fit[0], fit[1], fit[1]};
        }

        // Extract parameters
        double l1 = fit[1].get(0, 0);
        double l2 = fit[1].get(1, 1);
        double r1 = fit[0].get(0, 1);
        double r2 = fit[0].get(1, 0);
        double t = t3;
        double a1 = ai[0];

        // Calculate coefficients
        double w0 = (2 * r1 * (1 - Math.exp(-(r1 + r2) * t) - (r1 + r2) * t) *
                (Math.pow(a1, 2) * (r1 + r2) - a1 * r2 * (l1 - l2))) /
                (r2 * Math.pow(r1 + r2, 3));

        double w1 = -(2 * r1 * (1 - Math.exp(-(r1 + r2) * t) - (r1 + r2) * t) *
                (2 * a1 * l2 * (r1 + r2) - l2 * r2 * (l1 - l2))) /
                (r2 * Math.pow(r1 + r2, 3));

        double w2 = ((2 * Math.pow(l2, 2) * r2 * t) * (r1 + r2) +
                2 * Math.pow(l2, 2) * r1 * (1 - Math.exp(-(r1 + r2) * t))) /
                (r2 * Math.pow(r1 + r2, 2)) - (2 * Math.pow(l2, 2) * t) / r2;

        double w3 = (r1 + r2) / (l2 * r1);
        double w4 = (l1 * r2) / (l2 * r1);

        // Bounds for the first and second root
        double L1 = GlobalConstants.NegInf;
        double L2 = GlobalConstants.NegInf;
        double U1 = GlobalConstants.Inf;
        double U2 = GlobalConstants.Inf;

        // If set to true, the first (second) root is never feasible
        boolean infeasible1 = false;
        boolean infeasible2 = false;

        // Defined for convenience
        double z = w0 - Math.pow(w1, 2) / (4 * w2);

        // Impose square root argument is >= 0
        if (w2 > 0) {
            L1 = Math.max(L1, z);
            L2 = Math.max(L2, z);
        } else if (w2 < 0) {
            U1 = Math.min(U1, z);
            U2 = Math.min(U2, z);
        }

        // Impose q2 >= 0
        // First root
        if (w1 >= 0) {
            L1 = Math.max(L1, w0);
        } else if (w2 < 0) {
            infeasible1 = true;
        }
        // Second root
        if (w1 <= 0) {
            U2 = Math.min(U2, w0);
        } else if (w2 > 0) {
            infeasible2 = true;
        }

        // Impose q1 >= 0
        double tmp = 2 * a1 * w3 * w2 + w1;
        // First root
        if (tmp >= 0) {
            U1 = Math.min(U1, z + Math.pow(tmp, 2) / (4 * w2));
        } else if (w2 > 0) {
            infeasible1 = true;
        }
        // Second root
        if (tmp <= 0) {
            L2 = Math.max(L2, z + Math.pow(tmp, 2) / (4 * w2));
        } else if (w2 < 0) {
            infeasible2 = true;
        }

        // Impose q2 <= 1
        tmp = 2 * w2 + w1;
        // First root
        if (tmp >= 0) {
            U1 = Math.min(U1, z + Math.pow(tmp, 2) / (4 * w2));
        } else if (w2 > 0) {
            infeasible1 = true;
        }
        // Second root
        if (tmp <= 0) {
            L2 = Math.max(L2, z + Math.pow(tmp, 2) / (4 * w2));
        } else if (w2 < 0) {
            infeasible2 = true;
        }

        // Impose q1 <= 1
        tmp = 2 * a1 * w2 * w3 - 2 * w2 * w4 + w1;
        // First root
        if (tmp >= 0) {
            L1 = Math.max(L1, z + Math.pow(tmp, 2) / (4 * w2));
        } else if (w2 < 0) {
            infeasible1 = true;
        }
        // Second root
        if (tmp <= 0) {
            U2 = Math.min(U2, z + Math.pow(tmp, 2) / (4 * w2));
        } else if (w2 > 0) {
            infeasible2 = true;
        }

        if (infeasible1 && infeasible2) {
            throw new IllegalStateException("Empty feasibility region. This should not happen.");
        }

        // Compute feasible covariance
        double sigma;
        int root;

        if (infeasible2) {
            sigma = Math.max(Math.min(st3, U1), L1);
            root = 1;
        } else if (infeasible1) {
            sigma = Math.max(Math.min(st3, U2), L2);
            root = 2;
        } else {
            double sigma1 = Math.max(Math.min(st3, U1), L1);
            double sigma2 = Math.max(Math.min(st3, U2), L2);
            if (Math.abs(sigma1 - st3) < Math.abs(sigma2 - st3)) {
                sigma = sigma1;
                root = 1;
            } else {
                sigma = sigma2;
                root = 2;
            }
        }

        // Compute parameters
        double q2;
        if (root == 1) {
            q2 = (-w1 + Math.sqrt(Math.pow(w1, 2) - 4 * w2 * (w0 - sigma))) / (2 * w2);
        } else {
            q2 = (-w1 - Math.sqrt(Math.pow(w1, 2) - 4 * w2 * (w0 - sigma))) / (2 * w2);
        }

        double q1 = (a1 * (r1 + r2) - l2 * q2 * r1) / (l1 * r2);

        // Check feasibility
        double tol = 1e-8;
        double q1Final;
        double q2Final;

        if (q1 >= -tol && q1 <= 1 + tol && q2 >= -tol && q2 <= 1 + tol) {
            q1Final = Math.min(Math.max(q1, 0.0), 1.0);
            q2Final = Math.min(Math.max(q2, 0.0), 1.0);
        } else {
            throw new IllegalStateException("Parameters are infeasible. This should not happen.");
        }

        // Assemble M3PP[2]
        Matrix D0 = fit[0];
        Matrix D1 = fit[1];

        Matrix D1_class1 = new Matrix(2, 2);
        D1_class1.set(0, 0, D1.get(0, 0) * q1Final);
        D1_class1.set(1, 1, D1.get(1, 1) * q2Final);

        Matrix D1_class2 = new Matrix(2, 2);
        D1_class2.set(0, 0, D1.get(0, 0) * (1 - q1Final));
        D1_class2.set(1, 1, D1.get(1, 1) * (1 - q2Final));

        return new Matrix[]{D0, D1, D1_class1, D1_class2};
    }
}
