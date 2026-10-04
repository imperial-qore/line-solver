/**
 * @file PANACEA approximation
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc;

import org.apache.commons.math3.util.FastMath;

import jline.GlobalConstants;
import jline.io.InputOutput;
import jline.io.Ret;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;

public final class Pfqn_panacea {
    private Pfqn_panacea() {}

    public static Ret.pfqnNc pfqn_panacea(Matrix L, Matrix N, Matrix Z) {
        return pfqn_panacea(L, N, Z, new SolverOptions());
    }

    /**
     * Compute the PANACEA approximation
     */
    public static Ret.pfqnNc pfqn_panacea(Matrix L, Matrix N, Matrix Z, SolverOptions options) {
        Matrix Zlocal = Z;
        String method = options.method;
        int M = L.getNumRows();
        int R = L.getNumCols();
        double lG = Double.NaN;
        double G = Double.NaN;
        if (Zlocal.isEmpty() || Zlocal.elementSum() < options.tol) {
            Zlocal = Zlocal.copy();
            Zlocal = Zlocal.add(GlobalConstants.FineTol, Zlocal);
        }

        if (L.isEmpty() || L.elementSum() < options.tol) {
            Matrix tmp1 = Zlocal.sumCols();
            for (int i = 0; i < tmp1.length(); i++) {
                tmp1.set(i, N.get(i) * FastMath.log(tmp1.get(i)));
            }
            lG = -Matrix.factln(N).elementSum() + tmp1.elementSum();
            G = FastMath.exp(lG);
            return new Ret.pfqnNc(G, lG, method);
        }

        Matrix r = new Matrix(M, R);
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < R; j++) {
                r.set(i, j, L.get(i, j) / Zlocal.get(0, j));
            }
        }

        // see _kb/03-api-layer.md for rationale
        double Nt = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < R; j++) {
                if (r.get(i, j) > 0) {
                    Nt = FastMath.max(Nt, 1.0 / r.get(i, j));
                }
            }
        }

        Matrix beta = new Matrix(1, R);
        for (int j = 0; j < R; j++) {
            beta.set(0, j, N.get(0, j) / Nt);
        }

        Matrix gamma = new Matrix(M, R);
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < R; j++) {
                gamma.set(i, j, r.get(i, j) * Nt);
            }
        }

        Matrix alpha = new Matrix(1, M);
        for (int i = 0; i < M; i++) {
            double sum = 0.0;
            for (int j = 0; j < R; j++) {
                sum += N.get(j) * r.get(i, j);
            }
            alpha.set(0, i, 1 - sum);
        }

        Matrix gammatilde = new Matrix(M, R);
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < R; j++) {
                gammatilde.set(i, j, gamma.get(i, j % R) / alpha.get(0, i));
            }
        }

        double minAlpha = Double.MAX_VALUE;
        for (int i = 0; i < M; i++) {
            minAlpha = FastMath.min(minAlpha, alpha.get(0, i));
        }
        if (minAlpha < 0) {
            InputOutput.line_warning("pfqn_panacea", "Model is not in normal usage");
            return new Ret.pfqnNc(Double.NaN, Double.NaN, method);
        }

        double A0 = 1.0;
        double A1 = 0.0;
        for (int j = 0; j < R; j++) {
            Matrix m = new Matrix(1, R);
            m.set(0, j, 2.0);
            A1 -= beta.get(j) * Pfqn_ca.pfqn_ca(gammatilde, m).G;
        }

        double A2 = 0.0;
        for (int j = 0; j < R; j++) {
            Matrix m = new Matrix(1, R);
            m.set(0, j, 3.0);
            A2 += 2 * beta.get(j) * Pfqn_ca.pfqn_ca(gammatilde, m).G;

            m.set(0, j, 4.0);
            A2 += 3 * FastMath.pow(beta.get(j), 2) * Pfqn_ca.pfqn_ca(gammatilde, m).G;

            for (int k = 0; k < R; k++) {
                if (k == j) continue;
                m.zero();
                m.set(0, j, 2.0);
                m.set(0, k, 2.0);
                A2 = A2 + 0.5 * beta.get(j) * beta.get(k) * Pfqn_ca.pfqn_ca(gammatilde, m).G;
            }
        }

        double[] I = new double[]{A0, A1 / Nt, A2 / (Nt * Nt)};
        double sumI = 0.0;
        for (double v : I) sumI += v;

        Matrix tmp1 = Zlocal.sumCols();
        for (int i = 0; i < tmp1.length(); i++) {
            tmp1.set(i, N.get(i) * FastMath.log(tmp1.get(i)));
        }
        lG = -Matrix.factln(N).elementSum() + tmp1.elementSum();
        lG += FastMath.log(sumI);
        for (int s = 0; s < M; s++) {
            lG -= FastMath.log(alpha.get(0, s));
        }
        G = FastMath.exp(lG);

        if (!Double.isFinite(lG)) {
            G = Double.NaN;
            lG = Double.NaN;
        }
        return new Ret.pfqnNc(G, lG, method);
    }
}
