package jline.api.map;

import jline.util.matrix.Matrix;

/**
 * MAP/M/1-PS Sojourn Time Distribution
 *
 * Computes the complementary distribution function of sojourn time in a
 * MAP/M/1 processor-sharing queue using the algorithm from:
 *
 * Masuyama, H., and Takine, T. (2003). Sojourn time distribution in a
 * MAP/M/1 processor-sharing queue. Operations Research Letters, 31(6), 406-412.
 */
public final class MAPM1PSCdfRespT {
    private MAPM1PSCdfRespT() {}

    public static double[] computeCdf(Matrix C, Matrix D, double mu, double[] x) {
        return computeCdf(C, D, mu, x, 1e-11, 1e-10);
    }

    public static double[] computeCdf(Matrix C, Matrix D, double mu, double[] x, double epsilon) {
        return computeCdf(C, D, mu, x, epsilon, 1e-10);
    }

    public static double[] computeCdf(Matrix C, Matrix D, double mu, double[] x,
                                      double epsilon, double epsilonPrime) {
        int M = C.getNumRows();
        if (M != C.getNumCols()) {
            throw new IllegalArgumentException("C must be square");
        }
        if (M != D.getNumRows() || M != D.getNumCols()) {
            throw new IllegalArgumentException("C and D must have same dimensions");
        }
        if (!(mu > 0)) {
            throw new IllegalArgumentException("Service rate mu must be positive");
        }

        Matrix pi = computeStationaryDistribution(C, D);

        Matrix e = Matrix.ones(M, 1);
        double lambda = pi.mult(D).mult(e).get(0, 0);

        double rho = lambda / mu;
        if (!(rho < 1.0)) {
            throw new IllegalArgumentException("System is unstable (rho = " + rho + " >= 1)");
        }

        Matrix R = computeRMatrix(C, D, mu);

        Matrix I = Matrix.eye(M);
        Matrix pi0 = pi.mult(I.sub(1.0, R));

        int Nepsilon = determineNEpsilon(pi0, R, D, lambda, epsilon, e);

        double theta = 0.0;
        for (int i = 0; i < M; i++) {
            double v = Math.abs(C.get(i, i));
            if (v > theta) {
                theta = v;
            }
        }

        double[] result = new double[x.length];

        for (int idx = 0; idx < x.length; idx++) {
            double xVal = x[idx];

            double thetaPlusMu = theta + mu;
            double meanVal = thetaPlusMu * xVal;

            int L;
            int Kmax;

            if (meanVal > 0) {
                L = Math.max(0, (int) Math.floor(meanVal - 10 * Math.sqrt(meanVal)));
                Kmax = (int) Math.ceil(meanVal + 10 * Math.sqrt(meanVal));
            } else {
                L = 0;
                Kmax = 0;
            }

            Matrix[][] h = computeHRecursive(C, D, mu, Nepsilon, Kmax, theta);

            double WBar = 0.0;

            for (int n = 0; n <= Nepsilon; n++) {
                Matrix weight;
                if (n == 0) {
                    weight = pi0.mult(D);
                } else {
                    weight = pi0.mult(Matrix.pow(R, n)).mult(D);
                }

                Matrix sumK = Matrix.zeros(M, 1);
                for (int k = L; k <= Kmax; k++) {
                    double poissonTerm = poissonPmf(k, meanVal);
                    sumK = sumK.add(1.0, h[n][k].scale(poissonTerm));
                }

                WBar += (1.0 / lambda) * weight.mult(sumK).get(0, 0);
            }

            result[idx] = WBar;
        }

        return result;
    }

    private static Matrix computeStationaryDistribution(Matrix C, Matrix D) {
        int M = C.getNumRows();
        Matrix Q = C.add(1.0, D);

        Matrix A = Q.transpose();

        for (int j = 0; j < M; j++) {
            A.set(M - 1, j, 1.0);
        }

        Matrix b = Matrix.zeros(M, 1);
        b.set(M - 1, 0, 1.0);

        Matrix piT = Matrix.zeros(M, 1);
        Matrix.solve(A, b, piT);
        return piT.transpose();
    }

    private static Matrix computeRMatrix(Matrix C, Matrix D, double mu) {
        int M = C.getNumRows();

        if (M == 1) {
            double a = mu;
            double b = C.get(0, 0) - mu;
            double c = D.get(0, 0);

            double discriminant = b * b - 4 * a * c;
            if (discriminant < 0) {
                throw new IllegalArgumentException("No real solution for R matrix");
            }

            double R1 = (-b - Math.sqrt(discriminant)) / (2 * a);
            double R2 = (-b + Math.sqrt(discriminant)) / (2 * a);

            double Rval;
            if (R1 >= 0 && R1 < 1) {
                Rval = R1;
            } else if (R2 >= 0 && R2 < 1) {
                Rval = R2;
            } else {
                throw new IllegalStateException("No valid R in [0,1)");
            }

            Matrix result = new Matrix(1, 1);
            result.set(0, 0, Rval);
            return result;
        } else {
            Matrix I = Matrix.eye(M);
            Matrix CminusMuI = C.sub(mu, I);

            Matrix R = D.scale(-1.0).mult(CminusMuI.inv());

            int maxIter = 1000;
            double tol = 1e-12;

            for (int iter = 0; iter < maxIter; iter++) {
                Matrix Rold = R.copy();

                Matrix F = D.add(1.0, R.mult(CminusMuI)).add(mu, R.mult(R));
                Matrix Fprime = CminusMuI.add(2 * mu, R);

                R = R.sub(1.0, F.mult(Fprime.inv()));

                if (R.sub(1.0, Rold).elementMaxAbs() < tol) {
                    break;
                }
            }

            for (int i = 0; i < M; i++) {
                for (int j = 0; j < M; j++) {
                    if (R.get(i, j) < 0) {
                        R.set(i, j, 0.0);
                    }
                }
            }

            return R;
        }
    }

    private static int determineNEpsilon(Matrix pi0, Matrix R, Matrix D,
                                          double lambda, double epsilon, Matrix e) {
        int M = R.getNumRows();
        double cumsumProb = 0.0;
        Matrix Rpower = Matrix.eye(M);

        for (int n = 0; n <= 1000; n++) {
            cumsumProb += (1.0 / lambda) * pi0.mult(Rpower).mult(D).mult(e).get(0, 0);
            if (cumsumProb > 1 - epsilon) {
                return n;
            }
            Rpower = Rpower.mult(R);
        }

        return 100;
    }

    private static Matrix[][] computeHRecursive(Matrix C, Matrix D, double mu,
                                                  int N, int K, double theta) {
        int M = C.getNumRows();
        Matrix I = Matrix.eye(M);
        Matrix e = Matrix.ones(M, 1);

        double thetaPlusMu = theta + mu;
        Matrix thetaIPlusC = I.scale(theta).add(1.0, C);

        Matrix[][] h = new Matrix[N + 1][K + 1];
        for (int n = 0; n <= N; n++) {
            for (int k = 0; k <= K; k++) {
                h[n][k] = Matrix.zeros(M, 1);
            }
        }

        for (int n = 0; n <= N; n++) {
            h[n][0] = e.copy();
        }

        for (int k = 0; k < K; k++) {
            for (int n = 0; n <= N; n++) {
                Matrix term1 = Matrix.zeros(M, 1);
                Matrix term2 = thetaIPlusC.mult(h[n][k]);
                Matrix term3 = Matrix.zeros(M, 1);

                if (n > 0) {
                    term1 = h[n - 1][k].scale(((double) (n * mu)) / (n + 1));
                }

                if (n < N) {
                    term3 = D.mult(h[n + 1][k]);
                }

                h[n][k + 1] = term1.add(1.0, term2).add(1.0, term3).scale(1.0 / thetaPlusMu);
            }
        }

        return h;
    }

    private static double poissonPmf(int k, double lambda) {
        if (lambda == 0.0) {
            return (k == 0) ? 1.0 : 0.0;
        }

        double logProb = k * Math.log(lambda) - lambda;

        for (int i = 1; i <= k; i++) {
            logProb -= Math.log((double) i);
        }

        return Math.exp(logProb);
    }
}
