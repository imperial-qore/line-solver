/**
 * @file M/G/1 queueing system analysis with FB/LAS scheduling
 *
 * @since LINE 3.0
 */
package jline.api.qsys;

import jline.io.Ret;
import jline.util.matrix.Matrix;

public final class Qsys_mg1_fb {
    private Qsys_mg1_fb() {}

    /**
     * Analyzes an M/G/1 queueing system with FB (Feedback/LAS) scheduling.
     */
    public static Ret.qsys_prio qsys_mg1_fb(Matrix lambda, Matrix mu, Matrix cs) {
        double[] lambdaArr = lambda.toArray1D();
        double[] muArr = mu.toArray1D();
        double[] csArr = cs.toArray1D();

        if (lambdaArr.length != muArr.length || lambdaArr.length != csArr.length) {
            throw new IllegalArgumentException("lambda, mu, and cs must have the same length");
        }

        int K = lambdaArr.length;

        for (int i = 0; i < K; i++) {
            if (lambdaArr[i] <= 0.0) throw new IllegalArgumentException("lambda[" + i + "] must be positive");
            if (muArr[i] <= 0.0) throw new IllegalArgumentException("mu[" + i + "] must be positive");
            if (csArr[i] < 0.0) throw new IllegalArgumentException("cs[" + i + "] must be non-negative");
        }

        double[] rhoI = new double[K];
        double rhoTotal = 0.0;
        for (int i = 0; i < K; i++) {
            rhoI[i] = lambdaArr[i] / muArr[i];
            rhoTotal += rhoI[i];
        }
        if (rhoTotal >= 1.0) {
            throw new IllegalStateException("System is unstable: utilization rho = " + rhoTotal + " >= 1");
        }

        boolean allExp = true;
        for (int i = 0; i < K; i++) {
            if (Math.abs(csArr[i] - 1.0) >= 1e-6) {
                allExp = false;
                break;
            }
        }

        double[] W = new double[K];
        if (allExp) {
            // see _kb/03-api-layer.md for rationale
            for (int k = 0; k < K; k++) {
                double muK = muArr[k];
                double xMax = 20.0 / muK;
                int nIntervals = 4000; // even
                double h = xMax / nIntervals;
                double acc = 0.0;
                for (int s = 0; s <= nIntervals; s++) {
                    double x = s * h;
                    double t = fbResponseExp(x, lambdaArr, muArr);
                    double fx = muK * Math.exp(-muK * x);
                    double wgt = (s == 0 || s == nIntervals) ? 1.0 : (s % 2 == 1 ? 4.0 : 2.0);
                    acc += wgt * t * fx;
                }
                W[k] = acc * h / 3.0;
            }
        } else {
            // General-service heuristic (MATLAB qsys_mg1_fb_general):
            // evaluate the conditional response at x = 1/mu(k)
            for (int k = 0; k < K; k++) {
                double x = 1.0 / muArr[k];
                double rhoX = 0.0;
                for (int i = 0; i < K; i++) {
                    double integralFbar;
                    if (Math.abs(csArr[i] - 1.0) < 1e-10) {
                        integralFbar = (1.0 - Math.exp(-muArr[i] * x)) / muArr[i];
                    } else {
                        integralFbar = Math.min(x, 1.0 / muArr[i]);
                    }
                    rhoX += lambdaArr[i] * integralFbar;
                }

                double numerator = 0.0;
                for (int i = 0; i < K; i++) {
                    double integralTFbar;
                    if (Math.abs(csArr[i] - 1.0) < 1e-10) {
                        double muI = muArr[i];
                        integralTFbar = (1.0 - Math.exp(-muI * x) * (1.0 + muI * x)) / (muI * muI);
                    } else {
                        integralTFbar = Math.min(x * x / 2.0, 1.0 / (muArr[i] * muArr[i]));
                    }
                    numerator += lambdaArr[i] * integralTFbar;
                }

                if (rhoX >= 1.0) {
                    W[k] = Double.POSITIVE_INFINITY;
                } else {
                    double waitingTerm = numerator / ((1.0 - rhoX) * (1.0 - rhoX));
                    double serviceTerm = x / (1.0 - rhoX);
                    W[k] = waitingTerm + serviceTerm;
                }
            }
        }

        double Q = 0.0;
        for (int i = 0; i < K; i++) {
            Q += lambdaArr[i] * W[i];
        }
        double rhohat = Q / (1.0 + Q);

        return new Ret.qsys_prio(new Matrix(W), rhohat);
    }

    /**
     * FB conditional mean response time T(x) for exponential job sizes,
     * matching MATLAB compute_fb_response: T(x) = m(x)/(1-rho(x))^2 + x/(1-rho(x))
     * with rho(x), m(x) the truncated load and truncated second-moment term.
     */
    private static double fbResponseExp(double x, double[] lambdaArr, double[] muArr) {
        int K = lambdaArr.length;
        double rhoX = 0.0;
        double numerator = 0.0;
        for (int i = 0; i < K; i++) {
            double muI = muArr[i];
            double intFbar = (1.0 - Math.exp(-muI * x)) / muI;
            rhoX += lambdaArr[i] * intFbar;
            double intTFbar = (1.0 - Math.exp(-muI * x) * (1.0 + muI * x)) / (muI * muI);
            numerator += lambdaArr[i] * intTFbar;
        }
        if (rhoX >= 1.0) {
            return Double.POSITIVE_INFINITY;
        }
        return numerator / ((1.0 - rhoX) * (1.0 - rhoX)) + x / (1.0 - rhoX);
    }
}
