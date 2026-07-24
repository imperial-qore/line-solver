/**
 * @file Normalizing constant for mixed open-closed networks with limited load dependence.
 *
 * The closed-conditional normalizing constant of a mixed limited load-dependent (LLD)
 * network equals a purely closed load-dependent normalizing constant in which every
 * queueing station i carries the Bruell-Balbo-Afshari effective capacity rate
 * mu_i^eff(n) = 1/EC_i(n), where EC is returned by pfqn_ldmx_ec and folds the open
 * classes into the closed subnetwork. The open classes contribute the separable
 * prefactor lGopen = sum_i log E_i(0), which reduces to -sum_i log(1-rho_i) in the
 * load-independent limit. Mean closed metrics follow from the standard normalizing
 * constant ratios, e.g. X_r = G(N-e_r)/G(N), matching the exact pfqn_mvaldmx solver.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.ld;

import java.util.ArrayList;

import org.apache.commons.math3.util.FastMath;

import jline.io.Ret;
import jline.solvers.SolverOptions;
import jline.util.Utils;
import jline.util.matrix.Matrix;

public final class Pfqn_ncldmx {
    private Pfqn_ncldmx() {}

    /**
     * Normalizing constant for mixed open/closed networks with limited load dependence.
     *
     * @param lambda arrival rate vector (0 on closed classes)
     * @param D      service demand matrix (M x R)
     * @param N      population vector (Inf on open classes)
     * @param Z      think time vector (closed classes)
     * @param mu     load-dependent rate matrix (M x >= sum(N_closed))
     * @param S      number of servers per station (kept for signature parity)
     * @param options solver options forwarded to pfqn_ncld
     * @return closed-conditional constant, open prefactor and method
     */
    public static Ret.pfqnNcldmx pfqn_ncldmx(Matrix lambda, Matrix D, Matrix N, Matrix Z,
                                             Matrix mu, Matrix S, SolverOptions options) {
        int M = D.getNumRows();
        int R = D.getNumCols();
        if (Z == null) {
            Z = new Matrix(1, R);
        }
        ArrayList<Integer> openClasses = new ArrayList<Integer>();
        ArrayList<Integer> closedClasses = new ArrayList<Integer>();
        for (int r = 0; r < R; r++) {
            if (Utils.isInf(N.get(r))) {
                openClasses.add(r);
            } else {
                closedClasses.add(r);
            }
        }
        for (int r : closedClasses) {
            if (lambda.get(r) != 0 && N.get(r) > 0) {
                throw new RuntimeException("pfqn_ncldmx: Arrival rate cannot be specified on closed classes.");
            }
        }
        int C = closedClasses.size();
        int Kc = 0;
        for (int r : closedClasses) {
            Kc += (int) N.get(r);
        }

        // pad mu to at least max(1,Kc) columns, then append one extra column as in
        // pfqn_mvaldmx so that pfqn_ldmx_ec returns EC with Nt = max(1,Kc)+1 columns.
        int minCols = Math.max(1, Kc);
        int padCols = Math.max(mu.getNumCols(), minCols) + 1;
        Matrix mup = new Matrix(M, padCols);
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < padCols; j++) {
                int src = Math.min(j, mu.getNumCols() - 1);
                mup.set(i, j, mu.get(i, src));
            }
        }
        Ret.pfqnLDMXEC ec = Pfqn_ldmx_ec.pfqn_ldmx_ec(lambda, D, mup);
        Matrix EC = ec.EC;
        Matrix E = ec.E;
        double lGopen = 0.0;
        for (int i = 0; i < M; i++) {
            lGopen += FastMath.log(E.get(i, 0));
        }

        if (Kc == 0) {
            return new Ret.pfqnNcldmx(1.0, 0.0, lGopen, "exact");
        }

        Matrix Dc = new Matrix(M, C);
        Matrix Nc = new Matrix(1, C);
        Matrix Zc = new Matrix(1, C);
        Matrix muEff = new Matrix(M, Kc);
        for (int ci = 0; ci < C; ci++) {
            int r = closedClasses.get(ci);
            for (int i = 0; i < M; i++) {
                Dc.set(i, ci, D.get(i, r));
            }
            Nc.set(0, ci, N.get(r));
            Zc.set(0, ci, Z.get(r));
        }
        for (int i = 0; i < M; i++) {
            for (int k = 0; k < Kc; k++) {
                muEff.set(i, k, 1.0 / EC.get(i, k));
            }
        }
        Ret.pfqnNc cc = Pfqn_ncld.pfqn_ncld(Dc, Nc, Zc, muEff, options);
        return new Ret.pfqnNcldmx(cc.G, cc.lG, lGopen, cc.method);
    }
}
