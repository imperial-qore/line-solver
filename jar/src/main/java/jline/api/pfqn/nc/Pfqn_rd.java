/**
 * @file Reduction Heuristic (RD) method for load-dependent normalizing constants
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.apache.commons.math3.util.FastMath;

import jline.GlobalConstants;
import jline.api.pfqn.ld.Pfqn_gldsingle;
import jline.api.pfqn.mva.Pfqn_mva;
import jline.io.Ret;
import jline.solvers.SolverOptions;
import jline.solvers.nc.SolverNC;
import jline.util.Maths;
import jline.util.Utils;
import jline.util.matrix.Matrix;

public final class Pfqn_rd {
    private Pfqn_rd() {}

    public static Ret.pfqnRd pfqn_rd(Matrix L, Matrix N, Matrix Z, Matrix mu, SolverOptions options) {
        if (options == null) {
            options = SolverNC.defaultOptions();
        }
        String method = options.method;
        int M = L.getNumRows();
        int R = L.getNumCols();
        Matrix lambda = new Matrix(1, R);
        lambda.zero();
        if (N.elementSum() < 0) {
            return new Ret.pfqnRd(GlobalConstants.NegInf);
        }
        boolean isLi = true;
        for (int i = 0; i < M; i++) {
            for (int j = 1; j < R; j++) {
                if (mu.get(i, j) != mu.get(i, 0)) {
                    isLi = false;
                    break;
                }
            }
            if (isLi) {
                for (int j = 1; j < R; j++) {
                    L.set(i, j, L.get(i, j) / mu.get(i, 0));
                }
            }
        }
        if (N.elementSum() == 0.0) {
            return new Ret.pfqnRd(0.0);
        }

        Matrix gamma = new Matrix(M, (int) FastMath.ceil(N.elementSum()));
        gamma.ones();

        Matrix mu_new;
        if ((int) N.elementSum() >= mu.getNumCols()) {
            mu_new = mu.copy();
        } else {
            mu_new = new Matrix(mu.getNumRows(), 0);
            int i = 0;
            while (i < N.elementSum()) {
                Matrix mu_col_i = new Matrix(mu.getNumRows(), 1);
                Matrix.extract(mu, 0, mu.getNumRows(), i, i + 1, mu_col_i, 0, 0);
                mu_new = Matrix.concatColumns(mu_new, mu_col_i, null);
                i++;
            }
        }

        for (int i = 0; i < mu_new.getNumRows(); i++) {
            for (int j = 0; j < mu_new.getNumCols(); j++) {
                if (Double.isNaN(mu_new.get(i, j))) {
                    mu_new.set(i, j, GlobalConstants.Inf);
                }
            }
        }
        Matrix s = new Matrix(M, 1);
        for (int i = 0; i < M; i++) {
            if (mu.get(i, mu.getNumCols() - 1) != GlobalConstants.Inf) {
                Matrix mu_i = new Matrix(1, mu.getNumCols());
                for (int j = 0; j < mu.getNumCols(); j++) {
                    double mu_ij = FastMath.abs(mu.get(i, j) - mu.get(i, mu.getNumCols() - 1));
                    mu_i.set(j, ((mu_ij < options.tol) ? 1 : 0));
                }
                s.set(i, mu_i.find().elementMin());
            } else {
                s.set(i, N.elementSum() - 1);
            }
        }
        Matrix isDelay = new Matrix(new ArrayList<Double>(Collections.nCopies(M, 0.0)));
        new Matrix(new ArrayList<Double>(Collections.nCopies(M, 0.0)));
        Matrix y = L.copy();

        for (int i = 0; i < M; i++) {
            if (Utils.isInf(mu.get(i, (int) s.get(i)))) {
                for (int j = mu.getNumCols() - 1; j >= 0; j--) {
                    if (Double.isFinite(mu.get(i, j))) {
                        s.set(i, (double) j);
                        break;
                    }
                }
            }
            for (int j = 0; j < y.getNumCols(); j++) {
                y.set(i, j, y.get(i, j) / mu.get(i, (int) s.get(i)));
            }
        }
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < gamma.getNumCols(); j++) {
                gamma.set(i, j, mu.get(i, j) / mu.get(i, (int) s.get(i)));
            }
            double max = Double.NaN;
            double Ncum = 0.0;
            for (int j = 0; j < mu.getNumCols(); j++) {
                Ncum += 1.0;
                double x = FastMath.abs(mu.get(i, j) - Ncum);
                if (Double.isNaN(max) || x > max) {
                    max = x;
                }
            }
            if (max < options.tol) {
                isDelay.set(i, 1.0);
            }
        }
        Matrix beta = Matrix.ones(M, (int) FastMath.ceil(N.elementSum()));
        for (int i = 0; i < M; i++) {
            double beta_ij = gamma.get(i, 0) / (1 - gamma.get(i, 0));
            beta.set(i, 0, Double.isNaN(beta_ij) ? GlobalConstants.Inf : beta_ij);
            int j = 1;
            while (j < FastMath.ceil(N.elementSum())) {
                beta_ij = (1 - gamma.get(i, j - 1)) * (gamma.get(i, j) / (1 - gamma.get(i, j)));
                beta.set(i, j, Double.isNaN(beta_ij) ? GlobalConstants.Inf : beta_ij);
                j++;
            }
        }
        boolean isInf = true;
        int idx = 0;
        while (isInf && idx < beta.getNumElements()) {
            if (Double.isFinite(beta.get(idx))) {
                isInf = false;
            }
            idx++;
        }
        if (isInf) {
            options.method = "default";
            double lG = Pfqn_nc.pfqn_nc(lambda, L, N, Z, options).lG;
            options.method = method;
            return new Ret.pfqnRd(lG);
        }
        double Cgamma = 0.0;
        List<Double> sld = new ArrayList<Double>();
        for (int i = 0; i < s.getNumElements(); i++) {
            if (s.get(i) > 0) {
                sld.add(s.get(i));
            }
        }
        Matrix sldM = new Matrix(sld);
        int vmax = (int) FastMath.min(sldM.elementSum(), FastMath.ceil(N.elementSum()));
        Matrix Y = Pfqn_mva.pfqn_mva(y, N, Matrix.scaleMult(N, 0.0), Matrix.ones(1, M)).X;
        Matrix rhoN = y.mult(Y.transpose());
        Matrix lEN = new Matrix(1, vmax + 1);
        lEN.zero();
        for (int vtot = 0; vtot < vmax; vtot++) {
            lEN.set(vtot + 1, FastMath.log(Math.abs(
                    Pfqn_gldsingle.pfqn_gldsingle(rhoN, Matrix.singleton((double) (vtot + 1)), beta, options).G)));
        }
        double EN;
        for (int vtot = 0; vtot < lEN.getNumElements(); vtot++) {
            EN = FastMath.exp(lEN.get(vtot));
            Cgamma += ((N.elementSum() - Maths.max(0.0, (double) (vtot - 1))) / N.elementSum()) * EN;
        }
        options.method = "default";
        double lGN = Pfqn_nc.pfqn_nc(lambda, y, N, Z, options).lG;
        options.method = method;
        lGN += FastMath.log(Cgamma);
        return new Ret.pfqnRd(lGN, Cgamma);
    }
}
