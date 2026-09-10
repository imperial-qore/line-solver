/**
 * @file Normalizing constant and mean measures for mixed open-closed networks with limited load dependence.
 *
 * The closed-conditional normalizing constant of a mixed limited load-dependent (LLD)
 * network equals a purely closed load-dependent normalizing constant in which every
 * queueing station i carries the Bruell-Balbo-Afshari effective capacity rate
 * mu_i^eff(n) = 1/EC_i(n), where EC is returned by pfqn_ldmx_ec and folds the open
 * classes into the closed subnetwork. The open classes contribute the separable
 * prefactor lGopen = sum_i log E_i(0), which reduces to -sum_i log(1-rho_i) in the
 * load-independent limit.
 *
 * Mean measures follow from that identification without ever enumerating the closed
 * population lattice, which is what makes this the normalizing-constant counterpart of
 * pfqn_mvaldmx rather than a rename of it:
 *
 *   - closed throughputs are the ratios X_r = G(N-e_r)/G(N);
 *   - closed queue lengths are the conditional normalizing-constant recursion of the
 *     load-dependent closed network (pfqn_mushift / pfqn_fnc), applied to the
 *     effective-capacity rates;
 *   - open queue lengths are the Bruell-Balbo-Afshari sum
 *     Q_ir = lambda_r D_ir sum_n (n+1) EC_i(n+1) P_i(n) with its saturated tail folded
 *     onto the closed mean. EC_i(n) is constant for n &gt;= b_i, the level where the rate
 *     row stops growing, so writing EC_i(n) = EC_i^inf + delta_i(n) with delta_i(n) = 0
 *     for n &gt;= b_i leaves
 *
 *       Q_ir = lambda_r D_ir [ EC_i^inf (Q_i^closed + 1)
 *                              + sum_{n=0}^{b_i-2} (n+1) delta_i(n+1) P_i(n) ]
 *
 *     using sum_n P_i(n) = 1 and sum_n n P_i(n) = Q_i^closed. Only the first b_i-1
 *     marginal probabilities survive, and b_i is the number of servers, not the
 *     population: a single-server station needs none at all, and the formula collapses
 *     to the classical lambda_r D_ir (1+Q_i^closed)/(1-rho_i).
 *
 * The marginals that remain are themselves normalizing-constant ratios,
 * P_i(n) = sum_{|k|=n} F_i(k) G_{-i}(N-k) / G(N).
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.ld;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.util.FastMath;

import jline.io.Ret;
import jline.solvers.SolverOptions;
import jline.util.Maths;
import jline.util.Utils;
import jline.util.matrix.Matrix;

public final class Pfqn_ncldmx {
    private Pfqn_ncldmx() {}

    /**
     * Normalizing constant and mean measures for mixed open/closed networks with
     * limited load dependence.
     *
     * @param lambda arrival rate vector (0 on closed classes)
     * @param D      service demand matrix (M x R)
     * @param N      population vector (Inf on open classes)
     * @param Z      think time vector (closed classes)
     * @param mu     load-dependent rate matrix (M x &gt;= sum(N_closed))
     * @param S      number of servers per station (kept for signature parity)
     * @param options solver options forwarded to pfqn_ncld
     * @return closed-conditional constant, open prefactor, method, throughputs and queue lengths
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
        Matrix lambdao = new Matrix(1, R);
        lambdao.zero();
        for (int r : openClasses) {
            lambdao.set(r, lambda.get(r));
        }
        Ret.pfqnLDMXEC ec = Pfqn_ldmx_ec.pfqn_ldmx_ec(lambdao, D, mup);
        Matrix EC = ec.EC;
        Matrix E = ec.E;
        double lGopen = 0.0;
        for (int i = 0; i < M; i++) {
            lGopen += FastMath.log(E.get(i, 0));
        }

        int ncol = Math.max(1, Kc);
        Matrix Dc = new Matrix(M, C);
        Matrix Nc = new Matrix(1, C);
        Matrix Zc = new Matrix(1, C);
        Matrix muEff = new Matrix(M, ncol);
        for (int ci = 0; ci < C; ci++) {
            int r = closedClasses.get(ci);
            for (int i = 0; i < M; i++) {
                Dc.set(i, ci, D.get(i, r));
            }
            Nc.set(0, ci, N.get(r));
            Zc.set(0, ci, Z.get(r));
        }
        for (int i = 0; i < M; i++) {
            for (int k = 0; k < ncol; k++) {
                muEff.set(i, k, 1.0 / EC.get(i, k));
            }
        }

        double lG;
        double G;
        String method;
        if (Kc == 0) {
            lG = 0.0;
            G = 1.0;
            method = "exact";
        } else {
            Ret.pfqnNc cc = Pfqn_ncld.pfqn_ncld(Dc, Nc, Zc, muEff, options);
            lG = cc.lG;
            G = cc.G;
            method = cc.method;
        }

        // ---- mean measures ----
        Matrix XN = new Matrix(1, R);
        XN.zero();
        Matrix QN = new Matrix(M, R);
        QN.zero();
        for (int r : openClasses) {
            XN.set(r, lambda.get(r));
        }

        double[] lGr = new double[Math.max(C, 1)];
        if (Kc > 0) {
            for (int rc = 0; rc < C; rc++) {
                if (Nc.get(rc) <= 0) {
                    continue;
                }
                Matrix Ncr = Matrix.oner(Nc, new ArrayList<Integer>(Arrays.asList(rc)));
                lGr[rc] = ncldLocal(Dc, Ncr, Zc, muEff, options);
                XN.set(closedClasses.get(rc), FastMath.exp(lGr[rc] - lG));
            }
            // closed queue lengths: conditional normalizing-constant recursion of the
            // load-dependent closed network, on the effective-capacity rates
            for (int ist = 0; ist < M; ist++) {
                boolean anyDemand = false;
                for (int rc = 0; rc < C; rc++) {
                    if (Dc.get(ist, rc) > 0) {
                        anyDemand = true;
                        break;
                    }
                }
                if (!anyDemand) {
                    continue;
                }
                Matrix muhat = Pfqn_mushift.pfqn_mushift(muEff, ist);
                Ret.pfqnFnc fnc = Pfqn_fnc.pfqn_fnc(Matrix.extractRows(muhat, ist, ist + 1, null));
                Matrix muhatF = fnc.mu;
                double cshift = fnc.c.get(0);
                Matrix Dminus = removeRow(Dc, ist);
                Matrix muminus = removeRow(muEff, ist);
                Matrix DcPlus = Matrix.concatRows(Dc, Matrix.extractRows(Dc, ist, ist + 1, null), null);
                Matrix muhatPlus = Matrix.concatRows(muhat, muhatF, null);
                for (int rc = 0; rc < C; rc++) {
                    if (Nc.get(rc) <= 0 || Dc.get(ist, rc) <= 0) {
                        continue;
                    }
                    Matrix Ncr = Matrix.oner(Nc, new ArrayList<Integer>(Arrays.asList(rc)));
                    double lGhat = ncldLocal(Dc, Ncr, Zc, muhat, options);
                    double lGhatf = ncldLocal(DcPlus, Ncr, Zc, muhatPlus, options);
                    double lGminus = ncldLocal(Dminus, Ncr, Zc, muminus, options);
                    double CQ = (FastMath.exp(lGhatf - lGhat) - 1) + cshift * (FastMath.exp(lGminus - lGhat) - 1);
                    double ldDemand = FastMath.log(Dc.get(ist, rc)) + lGhat
                            - FastMath.log(muEff.get(ist, 0)) - lGr[rc];
                    QN.set(ist, closedClasses.get(rc),
                            FastMath.exp(ldDemand) * XN.get(closedClasses.get(rc)) * (1 + CQ));
                }
            }
        }

        // open queue lengths, with the saturated tail of EC folded onto the closed mean
        if (!openClasses.isEmpty()) {
            for (int ist = 0; ist < M; ist++) {
                double Qtot = 0.0;
                for (int rc = 0; rc < C; rc++) {
                    Qtot += QN.get(ist, closedClasses.get(rc));
                }
                int b = lldLevel(mup, ist);
                double ECinf = EC.get(ist, Math.min(b, EC.getNumCols()) - 1);
                double acc = ECinf * (Qtot + 1.0);
                if (b >= 2) {
                    Matrix Dminus = null;
                    Matrix muminus = null;
                    Matrix Drow = null;
                    Matrix murow = null;
                    if (Kc > 0) {
                        Dminus = removeRow(Dc, ist);
                        muminus = removeRow(muEff, ist);
                        Drow = Matrix.extractRows(Dc, ist, ist + 1, null);
                        murow = Matrix.extractRows(muEff, ist, ist + 1, null);
                    }
                    for (int n = 0; n <= b - 2; n++) {
                        double delta = EC.get(ist, n) - ECinf; // EC_i(n+1), one-based
                        if (delta == 0.0) {
                            continue;
                        }
                        // WITH NO CLOSED POPULATION THE MARGINAL IS DEGENERATE, not absent:
                        // P_i(0) = 1 and P_i(n) = 0 above it, so only the n = 0 term survives
                        // and acc collapses to EC_i(1), the exact open load-dependent mean.
                        // Skipping the loop instead left acc at EC_i^inf, i.e. read a c-server
                        // station as if every arrival found it saturated -- an M/M/3 at
                        // lambda = 1.5 came back with mean 1 against the exact 1.7368.
                        double Pn = Kc == 0
                                ? (n == 0 ? 1.0 : 0.0)
                                : marginal(n, Drow, murow, Nc, Zc, Dminus, muminus, lG, options);
                        acc += (n + 1) * delta * Pn;
                    }
                }
                for (int r : openClasses) {
                    QN.set(ist, r, lambda.get(r) * D.get(ist, r) * acc);
                }
            }
        }

        return new Ret.pfqnNcldmx(G, lG, lGopen, method, XN, QN);
    }

    /**
     * pfqn_ncld, with the empty-station residual network handled explicitly: a network
     * reduced to its think times alone has G(N) = prod_r Z_r^N_r / N_r!, and no G at all
     * when a class has jobs but neither demand nor think time.
     */
    private static double ncldLocal(Matrix L, Matrix N, Matrix Z, Matrix mu, SolverOptions options) {
        if (L.getNumRows() == 0) {
            boolean allZero = true;
            for (int r = 0; r < N.getNumElements(); r++) {
                if (N.get(r) > 0) {
                    allZero = false;
                    break;
                }
            }
            if (allZero) {
                return 0.0;
            }
            double acc = 0.0;
            for (int r = 0; r < N.getNumElements(); r++) {
                if (N.get(r) <= 0) {
                    continue;
                }
                double zr = 0.0;
                for (int i = 0; i < Z.getNumRows(); i++) {
                    zr += Z.get(i, r);
                }
                if (zr <= 0) {
                    return Double.NEGATIVE_INFINITY;
                }
                acc += N.get(r) * FastMath.log(zr) - Maths.factln(N.get(r));
            }
            return acc;
        }
        return Pfqn_ncld.pfqn_ncld(L, N, Z, mu, options).lG;
    }

    /** P_i(n) = sum_{|k|=n, k&lt;=Nc} F_i(k) G_{-i}(Nc-k) / G(Nc). */
    private static double marginal(int n, Matrix Drow, Matrix murow, Matrix Nc, Matrix Zc,
                                   Matrix Dminus, Matrix muminus, double lG, SolverOptions options) {
        int C = Nc.getNumElements();
        int[] cap = new int[C];
        for (int r = 0; r < C; r++) {
            cap[r] = (int) Nc.get(r);
        }
        List<int[]> ks = compositions(n, cap);
        Matrix Zzero = new Matrix(1, C);
        Zzero.zero();
        double P = 0.0;
        for (int[] k : ks) {
            Matrix kv = new Matrix(1, C);
            Matrix rest = new Matrix(1, C);
            for (int r = 0; r < C; r++) {
                kv.set(r, k[r]);
                rest.set(r, Nc.get(r) - k[r]);
            }
            double lF = (n == 0) ? 0.0 : ncldLocal(Drow, kv, Zzero, murow, options);
            double lGbar = ncldLocal(Dminus, rest, Zc, muminus, options);
            P += FastMath.exp(lF + lGbar - lG);
        }
        return P;
    }

    /** Non-negative integer vectors k with sum(k) == n and k &lt;= cap. */
    private static List<int[]> compositions(int n, int[] cap) {
        List<int[]> out = new ArrayList<int[]>();
        int C = cap.length;
        if (C == 0) {
            if (n == 0) {
                out.add(new int[0]);
            }
            return out;
        }
        int[] cur = new int[C];
        compositionsRec(n, cap, 0, cur, out);
        return out;
    }

    private static void compositionsRec(int rem, int[] cap, int idx, int[] cur, List<int[]> out) {
        if (idx == cap.length - 1) {
            if (rem <= cap[idx]) {
                int[] k = cur.clone();
                k[idx] = rem;
                out.add(k);
            }
            return;
        }
        int hi = Math.min(rem, cap[idx]);
        for (int v = 0; v <= hi; v++) {
            cur[idx] = v;
            compositionsRec(rem - v, cap, idx + 1, cur, out);
        }
        cur[idx] = 0;
    }

    /**
     * First column of the trailing constant run of a limited load-dependence row, i.e.
     * the level b with mu(n) = mu(b) for every n &gt;= b. This is the level pfqn_ldmx_ec
     * infers, and hence the one past which EC is constant. One-based.
     */
    private static int lldLevel(Matrix mu, int row) {
        int b = mu.getNumCols();
        if (b == 0) {
            return 1;
        }
        while (b > 1 && mu.get(row, b - 2) == mu.get(row, b - 1)) {
            b--;
        }
        return b;
    }

    /** A copy of A with row {@code skip} removed. */
    private static Matrix removeRow(Matrix A, int skip) {
        int rows = A.getNumRows();
        int cols = A.getNumCols();
        Matrix out = new Matrix(Math.max(rows - 1, 0), cols);
        out.zero();
        int w = 0;
        for (int i = 0; i < rows; i++) {
            if (i == skip) {
                continue;
            }
            for (int j = 0; j < cols; j++) {
                out.set(w, j, A.get(i, j));
            }
            w++;
        }
        return out;
    }
}
