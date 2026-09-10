/**
 * @file Parameter sanitization for product-form queueing network models
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc;

import java.util.Arrays;
import java.util.Comparator;

import org.apache.commons.math3.util.FastMath;

import jline.GlobalConstants;
import jline.io.Ret;
import jline.util.Maths;
import jline.util.matrix.Matrix;

public final class Pfqn_nc_sanitize {
    private Pfqn_nc_sanitize() {}

    /**
     * Sanitizes product-form model parameters to avoid degeneracies.
     */
    public static Ret.pfqnNcSanitize pfqn_nc_sanitize(final Matrix lambda, Matrix L, Matrix N, Matrix Z, double atol) {
        Matrix L_new = L.copy();
        Matrix Z_new = Z.copy();
        L_new.removeNaN();
        Z_new.removeNaN();
        Matrix L_tmp = new Matrix(L_new.getNumRows(), 0);
        Matrix N_tmp = new Matrix(N.getNumRows(), 0);
        Matrix Z_tmp = new Matrix(Z_new.getNumRows(), 0);
        Matrix lambda_tmp = new Matrix(lambda.getNumRows(), 0);
        for (int i = 0; i < N.length(); i++) {
            if (FastMath.abs(N.get(i)) >= GlobalConstants.FineTol
                    && !(L_new.sumCols(i) + Z_new.sumCols(i) < atol)) {
                Matrix L_col = Matrix.extractColumn(L_new, i, null);
                Matrix N_col = Matrix.extractColumn(N, i, null);
                Matrix Z_col = Matrix.extractColumn(Z_new, i, null);
                Matrix lambda_col = Matrix.extractColumn(lambda, i, null);
                L_tmp = Matrix.concatColumns(L_tmp, L_col, null);
                N_tmp = Matrix.concatColumns(N_tmp, N_col, null);
                Z_tmp = Matrix.concatColumns(Z_tmp, Z_col, null);
                lambda_tmp = Matrix.concatColumns(lambda_tmp, lambda_col, null);
            }
        }
        L_new = L_tmp;
        Z_new = Z_tmp;
        Matrix N_new = N_tmp;
        Matrix lambda_new = lambda_tmp;

        double lGremaind = 0.0;

        // see _kb/03-api-layer.md for rationale
        Matrix Z_zeroDemand = new Matrix(Z_new.getNumRows(), 0);
        Matrix N_zeroDemand = new Matrix(N_new.getNumRows(), 0);
        L_tmp = new Matrix(L_new.getNumRows(), 0);
        N_tmp = new Matrix(N_new.getNumRows(), 0);
        Z_tmp = new Matrix(Z_new.getNumRows(), 0);
        lambda_tmp = new Matrix(lambda_new.getNumRows(), 0);

        for (int i = 0; i < L_new.getNumCols(); i++) {
            Matrix L_col = Matrix.extractColumn(L_new, i, null);
            Matrix N_col = Matrix.extractColumn(N_new, i, null);
            Matrix Z_col = Matrix.extractColumn(Z_new, i, null);
            double LmaxCol = 0.0;
            for (int j = 0; j < L_new.getNumRows(); j++) {
                LmaxCol = FastMath.max(LmaxCol, L_new.get(j, i));
            }
            if (LmaxCol < atol) {
                Z_zeroDemand = Matrix.concatColumns(Z_zeroDemand, Z_col, null);
                N_zeroDemand = Matrix.concatColumns(N_zeroDemand, N_col, null);
            } else {
                L_tmp = Matrix.concatColumns(L_tmp, L_col, null);
                N_tmp = Matrix.concatColumns(N_tmp, N_col, null);
                Z_tmp = Matrix.concatColumns(Z_tmp, Z_col, null);
                if (lambda_new.getNumCols() > 0) {
                    lambda_tmp = Matrix.concatColumns(lambda_tmp,
                            Matrix.extractColumn(lambda_new, i, null), null);
                }
            }
        }

        for (int j = 0; j < N_zeroDemand.getNumCols(); j++) {
            double Nr = N_zeroDemand.get(0, j);
            double Zr = 0.0;
            for (int i = 0; i < Z_zeroDemand.getNumRows(); i++) {
                Zr += Z_zeroDemand.get(i, j);
            }
            lGremaind += Nr * FastMath.log(Zr) - Maths.factln(Nr);
        }
        L_new = L_tmp;
        Z_new = Z_tmp;
        N_new = N_tmp;
        lambda_new = lambda_tmp;

        if (!L_new.isEmpty()) {
            Matrix Lmax = new Matrix(1, L_new.getNumCols());
            for (int i = 0; i < Lmax.length(); i++) {
                Matrix L_col_i = Matrix.extractColumn(L_new, i, null);
                Lmax.set(i, L_col_i.elementMax());
            }
            if (Lmax.isEmpty()) {
                Lmax = new Matrix(1, Z_new.getNumCols());
                Lmax.ones();
            }
            Matrix repmat_Lmax_L = Lmax.repmat(L_new.getNumRows(), 1);
            Matrix repmat_Lmax_Z = Lmax.repmat(Z_new.getNumRows(), 1);
            for (int i = 0; i < L_new.getNumRows(); i++) {
                for (int j = 0; j < L_new.getNumCols(); j++) {
                    L_new.set(i, j, L_new.get(i, j) / repmat_Lmax_L.get(i, j));
                }
            }

            for (int i = 0; i < Z_new.getNumRows(); i++) {
                for (int j = 0; j < Z_new.getNumCols(); j++) {
                    Z_new.set(i, j, Z_new.get(i, j) / repmat_Lmax_Z.get(i, j));
                }
            }

            Matrix Lmax_log = Lmax.copy().transpose();
            for (int i = 0; i < Lmax_log.length(); i++) {
                Lmax_log.set(i, FastMath.log(Lmax_log.get(i)));
            }
            lGremaind += N_new.mult(Lmax_log).get(0);

            if (!Z_new.isEmpty()) {
                Integer[] index = new Integer[Z_new.getNumCols()];
                for (int i = 0; i < index.length; i++) {
                    index[i] = i;
                }

                final Matrix Z_sum_col = new Matrix(1, Z_new.getNumCols());
                for (int i = 0; i < Z_sum_col.length(); i++) {
                    Z_sum_col.set(i, Z_new.sumCols(i));
                }

                Arrays.sort(index, new Comparator<Integer>() {
                    @Override
                    public int compare(Integer i1, Integer i2) {
                        return Double.compare(Z_sum_col.get(i1), Z_sum_col.get(i2));
                    }
                });

                // see _kb/03-api-layer.md for rationale
                int[] order = new int[index.length];
                for (int i = 0; i < index.length; i++) {
                    order[i] = index[i];
                }
                L_new = permuteColumns(L_new, order);
                Z_new = permuteColumns(Z_new, order);
                N_new = permuteColumns(N_new, order);
                lambda_new = permuteColumns(lambda_new, order);
            }
        }

        // Ensure zero think time classes are anyway first. The test is on the
        // COLUMN SUM of Z, so it stays correct for multi-row delay matrices.
        if (!Z_new.isEmpty()) {
            int ncols = L_new.getNumCols();
            int[] order = new int[ncols];
            int pos = 0;
            for (int i = 0; i < ncols; i++) {
                if (Z_new.sumCols(i) < atol) {
                    order[pos++] = i;
                }
            }
            for (int i = 0; i < ncols; i++) {
                if (!(Z_new.sumCols(i) < atol)) {
                    order[pos++] = i;
                }
            }
            L_new = permuteColumns(L_new, order);
            N_new = permuteColumns(N_new, order);
            Z_new = permuteColumns(Z_new, order);
            lambda_new = permuteColumns(lambda_new, order);
        }

        return new Ret.pfqnNcSanitize(lambda_new, L_new, N_new, Z_new, lGremaind);
    }

    /**
     * Reorders the columns of A according to the given permutation. An argument
     * with no columns (an absent lambda, or an absent L) is returned unchanged,
     * which mirrors the MATLAB guards on empty arguments.
     */
    private static Matrix permuteColumns(Matrix A, int[] order) {
        if (A.getNumCols() == 0) {
            return A;
        }
        Matrix out = new Matrix(A.getNumRows(), 0);
        for (int i = 0; i < order.length; i++) {
            out = Matrix.concatColumns(out, Matrix.extractColumn(A, order[i], null), null);
        }
        return out;
    }
}
