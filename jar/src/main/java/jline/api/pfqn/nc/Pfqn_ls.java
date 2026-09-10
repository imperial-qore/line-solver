/**
 * @file Logistic sampling method for normalizing constant computation
 *
 * Implements the logistic sampling approach for computing normalizing constants in
 * closed product-form queueing networks. Uses importance sampling with multivariate
 * normal distributions centered at the Logistic expansion maximum.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc;

import org.apache.commons.math3.distribution.MultivariateNormalDistribution;
import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.util.FastMath;

import jline.GlobalConstants;
import jline.io.Ret;
import jline.util.Maths;
import jline.util.matrix.Matrix;

public final class Pfqn_ls {
    private Pfqn_ls() {}

    /**
     * Logistic sampling method to compute the normalizing constant.
     */
    public static Ret.pfqnNc pfqn_ls(Matrix L, Matrix N, Matrix Z, long I, long seed) {
        int M = L.getNumRows();
        int R = L.getNumCols();
        Matrix Lsum = new Matrix(M, 1);
        for (int i = 0; i < M; i++) {
            Lsum.set(i, L.sumRows(i));
        }
        Matrix L_new = new Matrix(0, R);
        for (int i = 0; i < M; i++) {
            Matrix L_row_i = new Matrix(1, R);
            Matrix.extract(L, i, i + 1, 0, R, L_row_i, 0, 0);
            if (Lsum.get(i) > GlobalConstants.CoarseTol) {
                L_new = Matrix.concatRows(L_new, L_row_i, null);
            }
        }
        M = L_new.getNumRows();
        R = L_new.getNumCols();
        double[][] sample = null;
        double lGn;

        if (L_new.isEmpty() || N.isEmpty() || N.elementSum() == 0.0
                || L_new.elementSum() < GlobalConstants.CoarseTol) {
            // Z can be empty here, and an empty class contributes 0, not 0*log(0).
            lGn = -Matrix.factln(N).elementSum();
            for (int r = 0; r < N.length(); r++) {
                if (N.get(r) > 0) {
                    lGn += N.get(r) * FastMath.log(Z.isEmpty() ? 0.0 : Z.sumCols(r));
                }
            }
        } else if (Z.isEmpty()) {
            Ret.pfqnLeFpi ret = Pfqn_le_fpi.pfqn_le_fpi(L_new, N);
            Matrix umax = ret.u;
            Matrix A = Pfqn_le_hessian.pfqn_le_hessian(L_new, N, umax.transpose());
            Matrix A_t = A.transpose();
            for (int i = 0; i < A.getNumRows(); i++) {
                for (int j = 0; j < A.getNumCols(); j++) {
                    A.set(i, j, (A.get(i, j) + A_t.get(i, j)) / 2.0);
                }
            }
            Matrix iA = A.inv();
            Matrix x0 = new Matrix(1, M - 1);
            for (int i = 0; i < M - 1; i++) {
                x0.set(i, FastMath.log(umax.get(i) / umax.get(M - 1)));
            }

            double[] x0_array = new double[M - 1];
            for (int i = 0; i < M - 1; i++) {
                x0_array[i] = x0.get(i);
            }
            double[][] iA_array = new double[iA.getNumRows()][];
            for (int i = 0; i < iA.getNumRows(); i++) {
                double[] tmp_row = new double[iA.getNumCols()];
                for (int j = 0; j < iA.getNumCols(); j++) {
                    tmp_row[j] = iA.get(i, j);
                }
                iA_array[i] = tmp_row;
            }
            MersenneTwister mt = new MersenneTwister((int) seed);
            MultivariateNormalDistribution mvd = new MultivariateNormalDistribution(mt, x0_array, iA_array);

            if (sample == null) {
                sample = mvd.sample((int) I);
            }

            // log-domain evaluation: exponentiating the integrand or the normal
            // density overflows/underflows once lG is large (hundreds of nats).
            double[] lr = new double[(int) I];
            double lrmax = Double.NEGATIVE_INFINITY;
            for (long i = 0; i < I; i++) {
                double[] sample_i = sample[(int) i];
                double lT = Maths.simplex_logfun(sample_i, L_new, N);
                double ldpdf = logmvnpdf_prec(sample_i, x0_array, A);
                lr[(int) i] = lT - ldpdf;
                if (lr[(int) i] > lrmax) {
                    lrmax = lr[(int) i];
                }
            }
            double meanexp = 0.0;
            for (int i = 0; i < (int) I; i++) {
                meanexp += FastMath.exp(lr[i] - lrmax);
            }
            meanexp /= I;

            Matrix tmp = new Matrix(1, N.length() + 1);
            Matrix.extract(N, 0, 1, 0, N.length(), tmp, 0, 0);
            tmp.set(N.length(), (double) (M - 1));
            lGn = Maths.multinomialln(tmp) + Maths.factln(M - 1) + lrmax + FastMath.log(meanexp);
        } else {
            Ret.pfqnLeFpiZ ret = Pfqn_le_fpiZ.pfqn_le_fpiZ(L_new, N, Z);
            Matrix umax = ret.u;
            double vmax = ret.v;
            Matrix A = Pfqn_le_hessianZ.pfqn_le_hessianZ(L_new, N, Z, umax.transpose(), vmax);
            Matrix A_t = A.transpose();
            for (int i = 0; i < A.getNumRows(); i++) {
                for (int j = 0; j < A.getNumCols(); j++) {
                    A.set(i, j, (A.get(i, j) + A_t.get(i, j)) / 2.0);
                }
            }
            Matrix iA = A.inv();
            Matrix x0 = new Matrix(1, M);
            for (int i = 0; i < M - 1; i++) {
                x0.set(i, FastMath.log(umax.get(i) / umax.get(M - 1)));
            }
            x0.set(M - 1, FastMath.log(vmax));

            double[] x0_array = new double[M];
            for (int i = 0; i < M; i++) {
                x0_array[i] = x0.get(i);
            }
            double[][] iA_array = new double[iA.getNumRows()][];
            for (int i = 0; i < iA.getNumRows(); i++) {
                double[] tmp_row = new double[iA.getNumCols()];
                for (int j = 0; j < iA.getNumCols(); j++) {
                    tmp_row[j] = iA.get(i, j);
                }
                iA_array[i] = tmp_row;
            }
            MersenneTwister mt = new MersenneTwister((int) seed);
            MultivariateNormalDistribution mvd = new MultivariateNormalDistribution(mt, x0_array, iA_array);

            if (sample == null) {
                sample = mvd.sample((int) I);
            }

            double epsilon = 1e-10;
            double eN = epsilon * N.elementSum();
            double eta = N.elementSum() + M * (1 + eN);
            int K = M;

            // see _kb/03-api-layer.md for rationale
            double[] lr = new double[(int) I];
            double lrmax = Double.NEGATIVE_INFINITY;
            for (long i = 0; i < I; i++) {
                double[] sample_i = sample[(int) i];
                double lT = pfqn_ls_helper(sample_i, K, M, eta, eN, L_new, N, Z);
                double ldpdf = logmvnpdf_prec(sample_i, x0_array, A);
                lr[(int) i] = lT - ldpdf;
                if (lr[(int) i] > lrmax) {
                    lrmax = lr[(int) i];
                }
            }
            double meanexp = 0.0;
            for (int i = 0; i < (int) I; i++) {
                meanexp += FastMath.exp(lr[i] - lrmax);
            }
            meanexp /= I;
            lGn = -Matrix.factln(N).elementSum() + lrmax + FastMath.log(meanexp);
        }
        double Gn = FastMath.exp(lGn);
        return new Ret.pfqnNc(Gn, lGn);
    }

    public static Ret.pfqnNc pfqn_ls(Matrix L, Matrix N, Matrix Z, long I) {
        return pfqn_ls(L, N, Z, I, 23000L);
    }

    public static Ret.pfqnNc pfqn_ls(Matrix L, Matrix N, long I) {
        return pfqn_ls(L, N, new Matrix(1, L.getNumCols()), I, 23000L);
    }

    public static Ret.pfqnNc pfqn_ls(Matrix L, Matrix N, long I, long seed) {
        return pfqn_ls(L, N, new Matrix(1, L.getNumCols()), I, seed);
    }

    /**
     * Auxiliary function used in the logistic sampling method.
     */
    static double pfqn_ls_helper(double[] x, int K, int M, double eta, double eN,
                                 Matrix L, Matrix N, Matrix Z) {
        double res = -Math.exp(x[K - 1]) + K * (1 + eN) * x[M - 1];
        double tmp1 = 0.0;
        double tmp2 = 0.0;
        for (int i = 0; i < K - 1; i++) {
            tmp1 += x[i];
            tmp2 += FastMath.exp(x[i]);
        }
        tmp2 = FastMath.log(tmp2 + 1);
        tmp2 *= -eta;
        res += (tmp1 + tmp2);

        Matrix L_row_K = new Matrix(1, L.getNumCols());
        Matrix L_first_K_minus_one_rows = new Matrix(K - 1, L.getNumCols());
        Matrix x_first_K_minus_one_elements = new Matrix(1, K - 1);
        Matrix.extract(L, K - 1, K, 0, L.getNumCols(), L_row_K, 0, 0);
        Matrix.extract(L, 0, K - 1, 0, L.getNumCols(), L_first_K_minus_one_rows, 0, 0);
        for (int i = 0; i < K - 1; i++) {
            x_first_K_minus_one_elements.set(i, FastMath.exp(x[i]));
        }
        for (int i = 0; i < K - 1; i++) {
            for (int j = 0; j < L.getNumCols(); j++) {
                L_first_K_minus_one_rows.set(i, j, L_first_K_minus_one_rows.get(i, j) * FastMath.exp(x[K - 1]) + Z.get(j));
            }
        }
        x_first_K_minus_one_elements = x_first_K_minus_one_elements.mult(L_first_K_minus_one_rows);
        for (int i = 0; i < L_row_K.length(); i++) {
            L_row_K.set(i, FastMath.log(L_row_K.get(i) * FastMath.exp(x[K - 1]) + Z.get(i)
                    + x_first_K_minus_one_elements.get(i)));
        }
        res += N.mult(L_row_K.transpose()).elementSum();

        // Return the log of the integrand (previously wrapped in exp(), which
        // overflows for large lG); the caller now works in the log domain.
        return res;
    }

    /**
     * Log-density of N(x0, inv(A)) evaluated at x, computed from the precision
     * matrix A to avoid multivariate-normal-density underflow at large lG.
     */
    static double logmvnpdf_prec(double[] x, double[] x0, Matrix A) {
        int d = x.length;
        double logdetA = FastMath.log(FastMath.abs(A.det()));
        double[] dvec = new double[d];
        for (int i = 0; i < d; i++) {
            dvec[i] = x[i] - x0[i];
        }
        double q = 0.0;
        for (int i = 0; i < d; i++) {
            double Adi = 0.0;
            for (int j = 0; j < d; j++) {
                Adi += A.get(i, j) * dvec[j];
            }
            q += dvec[i] * Adi;
        }
        return 0.5 * logdetA - (d / 2.0) * FastMath.log(2.0 * FastMath.PI) - 0.5 * q;
    }
}
