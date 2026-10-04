/**
 * @file Markovian Arrival Process autocorrelation decay rate analysis
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import org.apache.commons.math3.fitting.leastsquares.LeastSquaresBuilder;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.Pair;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Map_gamma {
    private Map_gamma() {}

    /**
     * Estimates the autocorrelation decay rate of a MAP.
     *
     * For MAPs of order higher than 2 the ACF is not geometric, so the decay rate
     * is obtained by fitting rho_k = RHO0*gamma^k in the least-squares sense, with
     * RHO0 = (1 - 1/SCV)/2 held fixed, mirroring MATLAB map_gamma.
     *
     * Distinct from {@link Map_gamma2}, which returns the second largest eigenvalue
     * of the embedded DTMC; the two agree for order 2 only.
     *
     * @param D0 the hidden transition matrix
     * @param D1 the visible transition matrix
     * @return the autocorrelation decay rate
     */
    public static double map_gamma(Matrix D0, Matrix D1) {
        return map_gamma(D0, D1, 1000);
    }

    /**
     * Estimates the autocorrelation decay rate of a MAP.
     *
     * @param D0 the hidden transition matrix
     * @param D1 the visible transition matrix
     * @param limit the maximum lag considered when fitting the ACF curve
     * @return the autocorrelation decay rate
     */
    public static double map_gamma(Matrix D0, Matrix D1, int limit) {
        int n = D0.getNumRows();

        if (n == 1) {
            // Poisson process: no correlation
            return 0.0;
        }

        if (n == 2) {
            // second-order MAP: geometric ACF
            double acf1 = Map_acf.map_acf(D0, D1, 1).value();
            if (FastMath.abs(acf1) < 1e-8) {
                // phase-type
                return 0.0;
            }
            return Map_acf.map_acf(D0, D1, 2).value() / acf1;
        }

        // higher-order MAP: the ACF is not geometric, so fit rho_k = RHO0*gamma^k
        // over the lags 1:(limit/10):limit, as MATLAB does
        int step = FastMath.max(1, limit / 10);
        int count = (limit - 1) / step + 1;
        final double[] lag = new double[count];
        Matrix lags = new Matrix(1, count);
        for (int i = 0; i < count; i++) {
            lag[i] = 1 + (double) i * step;
            lags.set(0, i, lag[i]);
        }

        double m1 = Map_mean.map_mean(D0, D1);
        double m2 = Map_moment.map_moment(D0, D1, 2);
        double scv = (m2 - m1 * m1) / (m1 * m1);
        final double rho0 = 0.5 * (1.0 - 1.0 / scv);

        Matrix acf = Map_acf.map_acf(D0, D1, lags);
        double[] rho = new double[count];
        for (int i = 0; i < count; i++) {
            rho[i] = acf.get(i);
        }

        return fitGeometricAcf(lag, rho, rho0, 0.99);
    }

    /**
     * Fits rho_k = rho0*gamma^k for gamma by robust nonlinear least squares.
     *
     * Mirrors MATLAB nlinfit with RobustWgtFun='fair': an ordinary least-squares
     * fit, then iteratively reweighted fits with the fair weight w=1/(1+|r|),
     * tuning constant 1.4, residuals adjusted by the leverage of the
     * least-squares Jacobian and scaled by a MAD estimate of sigma.
     *
     * @param lag the lags at which rho was evaluated
     * @param rho the autocorrelation at those lags
     * @param rho0 the fixed lag-0 coefficient (1 - 1/SCV)/2
     * @param start the initial guess for gamma
     * @return the fitted decay rate gamma
     */
    private static double fitGeometricAcf(double[] lag, double[] rho, double rho0, double start) {
        double gamma = leastSquares(lag, rho, rho0, start, null);

        // Leverage of the least-squares Jacobian, as advised by DuMouchel &
        // O'Brien. Held fixed across the reweighting, as nlinfit does. For a
        // single parameter the QR of an n-by-1 Jacobian reduces to normalising
        // it, so the leverage is just its squared unit entries.
        double norm2 = 0.0;
        for (int i = 0; i < lag.length; i++) {
            double j = rho0 * lag[i] * FastMath.pow(gamma, lag[i] - 1.0);
            norm2 += j * j;
        }
        double[] adjust = new double[lag.length];
        for (int i = 0; i < lag.length; i++) {
            double j = rho0 * lag[i] * FastMath.pow(gamma, lag[i] - 1.0);
            double h = norm2 > 0 ? FastMath.min(0.9999, j * j / norm2) : 0.0;
            adjust[i] = 1.0 / FastMath.sqrt(1.0 - h);
        }

        // A near-perfect fit would drive the MAD estimate of sigma to zero and
        // make every point an outlier, so floor it against the spread of the
        // response
        double mean = 0.0;
        for (int i = 0; i < rho.length; i++) mean += rho[i];
        mean /= rho.length;
        double var = 0.0;
        for (int i = 0; i < rho.length; i++) var += (rho[i] - mean) * (rho[i] - mean);
        var /= (rho.length - 1);
        double tinyS = 1e-6 * FastMath.sqrt(var);
        if (tinyS == 0.0) tinyS = 1.0;

        final double TUNE = 1.4; // fair
        double delta = FastMath.sqrt(Math.ulp(1.0));
        double[] weights = new double[lag.length];
        for (int iter = 0; iter < 200; iter++) {
            double previous = gamma;
            double[] radj = new double[lag.length];
            double[] absRadj = new double[lag.length];
            for (int i = 0; i < lag.length; i++) {
                radj[i] = (rho[i] - rho0 * FastMath.pow(gamma, lag[i])) * adjust[i];
                absRadj[i] = FastMath.abs(radj[i]);
            }
            // one parameter, so no residual is dropped from the MAD
            double sigma = median(absRadj) / 0.6745;
            double scale = FastMath.max(sigma, tinyS) * TUNE;
            for (int i = 0; i < lag.length; i++) {
                weights[i] = 1.0 / (1.0 + FastMath.abs(radj[i] / scale));
            }
            gamma = leastSquares(lag, rho, rho0, previous, weights);
            if (FastMath.abs(gamma - previous)
                    < delta * FastMath.max(FastMath.abs(gamma), FastMath.abs(previous))) {
                break;
            }
        }

        return gamma;
    }

    /**
     * Solves one (optionally weighted) least-squares fit of rho_k = rho0*gamma^k.
     *
     * @param lag the lags at which rho was evaluated
     * @param rho the autocorrelation at those lags
     * @param rho0 the fixed lag-0 coefficient
     * @param start the initial guess for gamma
     * @param weights the per-point weights, or null for an unweighted fit
     * @return the fitted decay rate gamma
     */
    private static double leastSquares(final double[] lag, final double[] rho, final double rho0,
                                       double start, final double[] weights) {
        final double[] sqrtW = new double[lag.length];
        for (int i = 0; i < lag.length; i++) {
            sqrtW[i] = weights == null ? 1.0 : FastMath.sqrt(weights[i]);
        }
        // The weights scale both the model and the target, so the residual the
        // optimizer minimises is sqrt(w)*(model - rho)
        final double[] target = new double[lag.length];
        for (int i = 0; i < lag.length; i++) {
            target[i] = sqrtW[i] * rho[i];
        }

        MultivariateJacobianFunction model = new MultivariateJacobianFunction() {
            @Override
            public Pair<RealVector, RealMatrix> value(RealVector point) {
                double gamma = point.getEntry(0);
                double[] value = new double[lag.length];
                double[][] jacobian = new double[lag.length][1];
                for (int i = 0; i < lag.length; i++) {
                    value[i] = sqrtW[i] * rho0 * FastMath.pow(gamma, lag[i]);
                    jacobian[i][0] = sqrtW[i] * rho0 * lag[i] * FastMath.pow(gamma, lag[i] - 1.0);
                }
                return new Pair<RealVector, RealMatrix>(new ArrayRealVector(value, false),
                        new Array2DRowRealMatrix(jacobian, false));
            }
        };

        LeastSquaresProblem problem = new LeastSquaresBuilder()
                .start(new double[]{start})
                .model(model)
                .target(target)
                .lazyEvaluation(false)
                .maxEvaluations(100000)
                .maxIterations(100000)
                .build();

        return new LevenbergMarquardtOptimizer().optimize(problem).getPoint().getEntry(0);
    }

    /**
     * Returns the median of the given values. The array is sorted in place.
     *
     * @param values the values
     * @return their median
     */
    private static double median(double[] values) {
        double[] sorted = values.clone();
        java.util.Arrays.sort(sorted);
        int n = sorted.length;
        if (n % 2 == 1) {
            return sorted[n / 2];
        }
        return 0.5 * (sorted[n / 2 - 1] + sorted[n / 2]);
    }

    /**
     * Estimates the autocorrelation decay rate of a MAP.
     *
     * @param MAP the MAP, as {D0,D1}
     * @return the autocorrelation decay rate, or 0 for a null MAP
     */
    public static double map_gamma(MatrixCell MAP) {
        if (MAP != null) {
            return map_gamma(MAP.get(0), MAP.get(1));
        } else {
            return 0.0;
        }
    }
}
