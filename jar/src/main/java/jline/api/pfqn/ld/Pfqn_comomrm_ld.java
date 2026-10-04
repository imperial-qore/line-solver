/**
 * @file COMOM normalizing constant method for load-dependent repairman models
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.ld;

import java.util.HashSet;
import java.util.Set;

import org.apache.commons.math3.util.FastMath;

import jline.GlobalConstants;
import jline.api.pfqn.nc.Pfqn_ca;
import jline.api.pfqn.nc.Pfqn_nc_sanitize;
import jline.io.Ret;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;

public final class Pfqn_comomrm_ld {
    private Pfqn_comomrm_ld() {}

    /**
     * Run the COMOM normalizing constant solution method on a repairman model.
     */
    public static Ret.pfqnComomrmLd pfqn_comomrm_ld(Matrix Lin, Matrix Nin, Matrix Zin, Matrix muIn, SolverOptions options) {
        Matrix L = Lin.copy();
        Matrix N = Nin.copy();
        Matrix Z = Zin.copy();
        Matrix mu = muIn.copy();
        double atol = options.tol;
        N.ceilEq();
        int M = L.getNumRows();
        int R = L.getNumCols();
        int Nt = (int) N.elementSum();
        Z = Z.sumCols();

        if (Z.elementSum() < GlobalConstants.Zero) {
            boolean[] zset = new boolean[M];
            Set<Integer> ignoreIndex = new HashSet<Integer>();
            boolean[] colset = new boolean[R];
            Matrix OneToNt = Matrix.createLike(mu.getRow(0));
            for (int nt = 1; nt <= OneToNt.length(); nt++) {
                OneToNt.set(nt - 1, (double) nt);
            }

            for (int j = 0; j < R; j++) {
                colset[j] = true;
            }

            for (int i = 0; i < M; i++) {
                if (mu.getRow(i).sub(OneToNt).norm() < atol) {
                    zset[i] = true;
                    ignoreIndex.add(Integer.valueOf(i));
                } else {
                    zset[i] = false;
                }
            }
            Z = L.getSlice(zset, colset);
            mu.removeRows(ignoreIndex);
            L.removeRows(ignoreIndex);
        }

        if (L.elementSum() < GlobalConstants.Zero) {
            Ret.pfqnNc ca = Pfqn_ca.pfqn_ca(L, N, Z);
            double G = ca.G;
            double lG = ca.lG;
            Matrix prob = new Matrix(Nt + 1, 1);
            prob.set(Nt, 1.0);
            return new Ret.pfqnComomrmLd(G, lG, prob);
        }

        double lG0;

        Ret.pfqnNcSanitize sanitizedResult = Pfqn_nc_sanitize.pfqn_nc_sanitize(new Matrix(1, R), L, N, Z, atol);
        L = sanitizedResult.L;
        N = sanitizedResult.N;
        Z = sanitizedResult.Z;
        lG0 = sanitizedResult.lGremaind;
        M = L.getNumRows();
        R = L.getNumCols();

        if (Z.isEmpty()) {
            if (L.isEmpty()) {
                double G = FastMath.exp(lG0);
                double lG = lG0;
                Matrix prob = new Matrix(Nt + 1, 1);
                prob.set(0, 1.0);
                return new Ret.pfqnComomrmLd(G, lG, prob);
            }
            Z = new Matrix(1, R);
        } else if (L.isEmpty()) {
            L = new Matrix(1, R);
        }

        if (M == 0) {
            double G = FastMath.exp(lG0);
            double lG = lG0;
            Matrix prob = new Matrix(Nt + 1, 1);
            prob.set(Nt, 1.0);
            return new Ret.pfqnComomrmLd(G, lG, prob);
        }

        if (M != 1) {
            throw new IllegalArgumentException("The solver accepts at most a single queueing station.");
        }

        Matrix h = new Matrix(Nt + 1, 1);
        h.set(Nt, 1.0);
        Matrix scale = new Matrix(Nt, 1);
        int nt = 0;

        for (int r = 0; r < R; r++) {
            Matrix Tr = Matrix.eye(Nt + 1).scale(Z.get(r));
            for (int i = 0; i < Tr.getNumCols() - 1; i++) {
                Tr.set(i, i + 1, L.get(r) * (Nt - i) / mu.get(Nt - i - 1));
            }
            Matrix hT;
            int nr = 0;
            while (nr < N.get(r)) {
                hT = Tr.copy().scale(1.0 / (1.0 + nr));
                h = hT.mult(h);
                scale.set(nt, FastMath.abs(h.sort().elementSum()));
                h.absEq();
                h.scaleEq(1.0 / scale.get(nt));
                nt++;
                nr++;
            }
        }

        double lG = lG0 + scale.log().elementSum();
        double G = FastMath.exp(lG);
        Matrix prob = h.reverse().scale(1 / G);
        prob.divideEq(prob.elementSum());

        return new Ret.pfqnComomrmLd(G, lG, prob);
    }
}
