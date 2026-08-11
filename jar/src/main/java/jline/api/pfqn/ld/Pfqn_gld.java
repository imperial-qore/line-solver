/**
 * @file Generalized load-dependent normalizing constant computation
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.ld;

import org.apache.commons.math3.util.FastMath;

import jline.GlobalConstants;
import jline.api.pfqn.nc.Pfqn_nc;
import jline.io.Ret;
import jline.solvers.SolverOptions;
import jline.solvers.nc.SolverNC;
import jline.util.Maths;
import jline.util.matrix.Matrix;

public final class Pfqn_gld {
    private Pfqn_gld() {}

    /**
     * Compute the normalizing constant of a single-class load-dependent closed queueing network model
     */
    public static Ret.pfqnNc pfqn_gld(Matrix L, Matrix N, Matrix mu, SolverOptions options) {
        int M = L.getNumRows();
        int R = L.getNumCols();
        Matrix lambda = new Matrix(1, R);
        Double G;
        double lG;

        if (M == 1) {
            Matrix N_tmp = new Matrix(1, 0);
            Matrix L_tmp = new Matrix(1, 0);
            for (int i = 0; i < R; i++) {
                if (L.get(i) > GlobalConstants.FineTol) {
                    Matrix N_tmp2 = new Matrix(1, 1);
                    N_tmp2.fill(N.get(i));
                    Matrix L_tmp2 = new Matrix(1, 1);
                    L_tmp2.fill(FastMath.log(L.get(0, i)));
                    N_tmp = Matrix.concatColumns(N_tmp, N_tmp2, null);
                    L_tmp = Matrix.concatColumns(L_tmp, L_tmp2, null);
                }
            }
            Matrix mu_new;
            if ((int) N.elementSum() >= mu.getNumCols()) {
                mu_new = Matrix.extractRows(mu, 0, 1, null);
            } else {
                mu_new = new Matrix(1, 0);
                int i = 0;
                while (i < N.elementSum()) {
                    Matrix mu_col_i = new Matrix(1, 1);
                    Matrix.extract(mu, 0, 1, i, i + 1, mu_col_i, 0, 0);
                    mu_new = Matrix.concatColumns(mu_new, mu_col_i, null);
                    i++;
                }
            }

            for (int i = 0; i < mu_new.length(); i++) {
                mu_new.set(i, FastMath.log(mu_new.get(i)));
            }

            lG = (Maths.factln(N.elementSum()) - Matrix.factln(N).elementSum()
                    + N_tmp.mult(L_tmp.transpose()).get(0) - mu_new.elementSum());
            G = FastMath.exp(lG);
            return new Ret.pfqnNc(G, lG);
        }

        if (R == 1) {
            Ret.pfqnNc ret = Pfqn_gldsingle.pfqn_gldsingle(L, N, mu, null);
            lG = ret.lG;
            G = ret.G;
            return new Ret.pfqnNc(G, lG);
        }

        if (L.isEmpty()) {
            G = 0.0;
            lG = Double.NEGATIVE_INFINITY;
            return new Ret.pfqnNc(G, lG);
        }

        Matrix mu_new;
        if (mu == null) {
            mu_new = new Matrix(M, (int) N.elementSum());
            mu_new.fill(1.0);
        } else {
            mu_new = mu.copy();
        }
        SolverOptions options_new = (options != null) ? options : SolverNC.defaultOptions();

        boolean isLoadDep = false;
        boolean[] isInfServer = new boolean[M];
        for (int i = 0; i < M; i++) {
            Matrix mu_row_i = new Matrix(1, (int) N.elementSum());
            Matrix.extract(mu_new, i, i + 1, 0, (int) N.elementSum(), mu_row_i, 0, 0);
            boolean flag = true;
            for (int j = 0; j < mu_row_i.getNumCols(); j++) {
                if (FastMath.abs(mu_row_i.get(j) - (j + 1)) > GlobalConstants.FineTol) {
                    flag = false;
                    break;
                }
            }

            if (FastMath.abs(mu_row_i.elementMin() - 1) < GlobalConstants.FineTol
                    && FastMath.abs(mu_row_i.elementMax() - 1) < GlobalConstants.FineTol) {
                isInfServer[i] = false;
            } else if (flag) {
                isInfServer[i] = true;
            } else {
                isInfServer[i] = false;
                isLoadDep = true;
            }
        }

        if (!isLoadDep) {
            Matrix Lli = new Matrix(0, L.getNumCols());
            Matrix Zli = new Matrix(0, L.getNumCols());
            for (int i = 0; i < M; i++) {
                Matrix L_row_i = Matrix.extractRows(L, i, i + 1, null);
                if (isInfServer[i]) {
                    Zli = Matrix.concatRows(Zli, L_row_i, null);
                } else {
                    Lli = Matrix.concatRows(Lli, L_row_i, null);
                }
            }
            if (Lli.isEmpty()) {
                Lli = N.copy();
                Lli.fill(0.0);
            }
            if (Zli.isEmpty()) {
                Zli = N.copy();
                Zli.fill(0.0);
            }
            options_new.method = "exact";
            lG = Pfqn_nc.pfqn_nc(lambda, Lli, N, Zli.sumCols(), options_new).lG;
            G = FastMath.exp(lG);
            return new Ret.pfqnNc(G, lG);
        }

        G = 0.0;
        if (M == 0) {
            lG = FastMath.log(G);
            return new Ret.pfqnNc(G, lG);
        }

        if (FastMath.abs(N.elementMax()) < GlobalConstants.FineTol
                && FastMath.abs(N.elementMin()) < GlobalConstants.FineTol) {
            G = 1.0;
            lG = FastMath.log(G);
            return new Ret.pfqnNc(G, lG);
        }

        if (R == 1) {
            G = Pfqn_gldsingle.pfqn_gldsingle(L, N, mu_new, null).G;
            lG = FastMath.log(G);
            return new Ret.pfqnNc(G, lG);
        }

        G = Pfqn_gld.pfqn_gld(Matrix.extractRows(L, 0, M - 1, null), N,
                Matrix.extractRows(mu_new, 0, M - 1, null), options_new).G;

        for (int r = 0; r < R; r++) {
            if (N.get(r) > GlobalConstants.FineTol) {
                Matrix N_1 = N.copy();
                if (R > 1) {
                    N_1.set(r, N_1.get(r) - 1);
                } else {
                    for (int i = 0; i < N_1.length(); i++) {
                        N_1.set(i, N_1.get(i) - 1);
                    }
                }
                G += L.get(M - 1, r) / mu_new.get(M - 1, 0)
                        * Pfqn_gld.pfqn_gld(L, N_1, Pfqn_mushift.pfqn_mushift(mu, M - 1), options_new).G;
            }
        }
        lG = FastMath.log(G);
        return new Ret.pfqnNc(G, lG);
    }
}
