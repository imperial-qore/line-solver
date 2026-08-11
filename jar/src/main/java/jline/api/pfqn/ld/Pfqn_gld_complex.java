/**
 * @file Complex-valued generalized load-dependent normalizing constant computation
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.ld;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.util.FastMath;

import jline.GlobalConstants;
import jline.io.Ret;
import jline.solvers.SolverOptions;
import jline.solvers.nc.SolverNC;
import jline.util.Maths;
import jline.util.matrix.ComplexMatrix;
import jline.util.matrix.Matrix;

public final class Pfqn_gld_complex {
    private Pfqn_gld_complex() {}

    /**
     * Compute the normalizing constant of a single-class load-dependent closed
     * queueing network model with complex demands.
     */
    public static Ret.pfqnNcComplex pfqn_gld_complex(ComplexMatrix L, Matrix N, Matrix mu, SolverOptions options) {
        int M = L.getNumRows();
        int R = L.getNumCols();
        Complex G;
        Complex lG;

        if (M == 1) {
            Matrix N_tmp = new Matrix(1, 0);
            ComplexMatrix L_tmp = new ComplexMatrix(1, 0);
            for (int i = 0; i < R; i++) {
                if (L.real.get(i) > GlobalConstants.FineTol || L.im.get(i) > GlobalConstants.FineTol) {
                    Matrix N_tmp2 = new Matrix(1, 1);
                    N_tmp2.fill(N.get(i));
                    ComplexMatrix L_tmp2 = new ComplexMatrix(1, 1);
                    L_tmp2.real.fill(0.5 * FastMath.log(Math.pow(L.real.get(0, i), 2.0) + FastMath.pow(L.im.get(0, i), 2)));
                    L_tmp2.im.fill(Math.log(Math.atan(L.im.get(0, i) / L.real.get(0, i))));
                    N_tmp = Matrix.concatColumns(N_tmp, N_tmp2, null);
                    L_tmp.real = Matrix.concatColumns(L_tmp.real, L_tmp2.real, null);
                    L_tmp.im = Matrix.concatColumns(L_tmp.im, L_tmp2.im, null);
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
            lG = new Complex(N_tmp.mult(L_tmp.real.transpose()).get(0),
                    N_tmp.mult(L_tmp.im.transpose()).get(0))
                    .add(Maths.factln(N.elementSum()) - Matrix.factln(N).elementSum())
                    .subtract(mu_new.elementSum());
            G = lG.exp();
            return new Ret.pfqnNcComplex(G, lG);
        }

        if (R == 1) {
            Ret.pfqnNcComplex ret = Pfqn_gldsingle_complex.pfqn_gldsingle_complex(L, N, mu, null);
            lG = ret.lG;
            G = ret.G;
            return new Ret.pfqnNcComplex(G, lG);
        }

        if (L.isEmpty()) {
            G = new Complex(0.0);
            lG = new Complex(GlobalConstants.NegInf);
            return new Ret.pfqnNcComplex(G, lG);
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
            throw new RuntimeException("pfqn_gld_complex is only implemented for load dependent models.");
        }

        G = new Complex(0.0);
        if (M == 0) {
            lG = G.log();
            return new Ret.pfqnNcComplex(G, lG);
        }

        if (FastMath.abs(N.elementMax()) < GlobalConstants.FineTol
                && FastMath.abs(N.elementMin()) < GlobalConstants.FineTol) {
            G = new Complex(1.0);
            lG = G.log();
            return new Ret.pfqnNcComplex(G, lG);
        }

        if (R == 1) {
            return Pfqn_gldsingle_complex.pfqn_gldsingle_complex(L, N, mu_new, null);
        }

        G = pfqn_gld_complex(ComplexMatrix.extractRows(L, 0, M - 1, null),
                N,
                Matrix.extractRows(mu_new, 0, M - 1, null),
                options_new).G;

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
                G = G.add(L.get(M - 1, r).divide(mu_new.get(M - 1, 0))
                        .multiply(pfqn_gld_complex(L, N_1, Pfqn_mushift.pfqn_mushift(mu, M - 1), options_new).G));
            }
        }
        lG = G.log();
        return new Ret.pfqnNcComplex(G, lG);
    }
}
