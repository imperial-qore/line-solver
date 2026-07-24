/**
 * Exact analytic performance sensitivities for closed product-form queueing
 * networks. Dispatches to a CoMoM-backed kernel for the single-station
 * repairman model and to differentiated Mean Value Analysis otherwise.
 *
 * <p>Derivatives of the mean measures {X, Q, U, R} with respect to the service
 * demands L(i,r) and think times Z(r) are analytic (exact to machine
 * precision), not finite differences. Queue-length variances and covariances
 * are returned as an exact by-product. Mirrors the MATLAB reference
 * {@code pfqn_sens.m} and the native-Python {@code pfqn_sens}.</p>
 *
 * <p>Reference: Z. Liu and P. Nain, INRIA RR-1144, 1989; X.-R. Cao and D.-J.
 * Ma, Performance Evaluation 26:181-199, 1996; G. Casale, "CoMoM: Efficient
 * Class-Oriented Evaluation of Multiclass Performance Models", IEEE TSE 2011.</p>
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.sens;

import java.util.HashMap;
import java.util.Map;

import org.apache.commons.math3.util.FastMath;

import jline.GlobalConstants;
import jline.api.pfqn.nc.Pfqn_comomrm;
import jline.io.Ret;
import jline.util.matrix.Matrix;

public final class Pfqn_sens {
    private Pfqn_sens() {}

    /**
     * Exact derivatives of {X, Q, U, R} w.r.t. demands L and think times Z for
     * a closed single-server (or residence-multiplicity mi) product-form
     * network with an infinite-server delay Z.
     *
     * @param L  service demand matrix (M x R), L(i,r) = visits_ir / rate_ir
     * @param N  population vector (1 x R)
     * @param Z  think time vector (1 x R), may be empty
     * @param mi station residence multiplicity (1 x M), null for ones
     * @return sensitivities (base measures, Jacobians, queue-length moments)
     */
    public static Ret.pfqnSens pfqn_sens(Matrix L, Matrix N, Matrix Z, Matrix mi) {
        int M = L.getNumRows();
        int R = L.getNumCols();

        N = N.copy().ceil();
        if (N.getNumRows() > 1) {
            N = N.transpose();
        }
        if (Z == null || Z.isEmpty()) {
            Z = new Matrix(1, R);
        } else if (Z.getNumRows() > 1) {
            Z = Z.transpose();
        }
        if (mi == null) {
            mi = Matrix.ones(1, M);
        } else if (mi.getNumRows() > 1) {
            mi = mi.transpose();
        }

        // see _kb/03-api-layer.md for rationale
        double tol = GlobalConstants.FineTol;
        boolean anyPop = false;
        boolean useComom = (M == 1);
        for (int i = 0; i < M && useComom; i++) {
            if (mi.get(0, i) != 1.0) {
                useComom = false;
            }
        }
        for (int r = 0; r < R; r++) {
            if (N.get(0, r) > 0) {
                anyPop = true;
                if (useComom && (!(Z.get(0, r) > tol) || !(L.get(0, r) > tol))) {
                    useComom = false;
                }
            }
        }
        useComom = useComom && anyPop;

        Ret.pfqnSens sens;
        if (useComom) {
            sens = sensComom(L, N, Z);
        } else {
            sens = sensMva(L, N, Z, mi);
        }
        attachSecondMoments(sens, L, N, Z, mi);
        return sens;
    }

    public static Ret.pfqnSens pfqn_sens(Matrix L, Matrix N, Matrix Z) {
        return pfqn_sens(L, N, Z, null);
    }

    // =====================================================================
    // CoMoM-backed kernel (M=1 repairman model)
    // =====================================================================
    private static Ret.pfqnSens sensComom(Matrix L, Matrix N, Matrix Z) {
        int R = L.getNumCols();
        double[] D = new double[R];
        double[] Zarr = new double[R];
        int[] Nn = new int[R];
        for (int r = 0; r < R; r++) {
            D[r] = L.get(0, r);
            Zarr[r] = Z.get(0, r);
            Nn[r] = (int) Math.round(N.get(0, r));
        }

        // parameter list: L(0,r) first, then Z(r) (mirrors sensMva ordering)
        int P = 2 * R;
        int[] paramType = new int[P];
        int[] paramStation = new int[P];
        int[] paramClass = new int[P];
        int[] pL = new int[R];
        int[] pZ = new int[R];
        for (int r = 0; r < R; r++) {
            pL[r] = r;
            paramType[pL[r]] = 0; paramStation[pL[r]] = 0; paramClass[pL[r]] = r;
            pZ[r] = R + r;
            paramType[pZ[r]] = 1; paramStation[pZ[r]] = -1; paramClass[pZ[r]] = r;
        }

        Matrix X = new Matrix(1, R);
        Matrix Q = new Matrix(1, R);
        Matrix U = new Matrix(1, R);
        Matrix C = new Matrix(1, R);
        Matrix dX = new Matrix(R, P);
        Matrix[] dQ = zeros(P, 1, R);
        Matrix[] dU = zeros(P, 1, R);
        Matrix[] dC = zeros(P, 1, R);

        boolean anyPop = false;
        for (int r = 0; r < R; r++) {
            if (Nn[r] > 0) { anyPop = true; break; }
        }
        if (!anyPop) {
            return new Ret.pfqnSens(X, Q, U, C, dX, dQ, dU, dC,
                    paramType, paramStation, paramClass);
        }

        ReplNC nc = new ReplNC(D, Zarr, R);

        // base measures
        for (int r = 0; r < R; r++) {
            if (Nn[r] < 1) {
                continue;
            }
            X.set(0, r, nc.xput(1, Nn, r));
        }
        for (int s = 0; s < R; s++) {
            Q.set(0, s, nc.qmean(Nn, s));
        }
        for (int r = 0; r < R; r++) {
            U.set(0, r, X.get(0, r) * D[r]);
            if (X.get(0, r) > 0) {
                C.set(0, r, Q.get(0, r) / X.get(0, r));
            }
        }

        // Jacobian
        for (int r = 0; r < R; r++) {
            if (Nn[r] < 1) {
                continue;   // empty class: X=Q=0, all derivatives 0
            }
            int[] Nr = Nn.clone();
            Nr[r] -= 1;
            double Xr = X.get(0, r);
            double Qr = Q.get(0, r);
            for (int s = 0; s < R; s++) {
                double drs = (r == s) ? 1.0 : 0.0;

                // L(0,s) parameter
                int p = pL[s];
                double Vrs = Qr * (drs + 2.0 * nc.qplus(Nr, s) - Q.get(0, s));  // Cov[n_r,n_s]
                double dQ_L = Vrs / D[s];
                double dX_L = Xr * (nc.qmean(Nr, s) - Q.get(0, s)) / D[s];
                double dU_L = dX_L * D[r] + Xr * drs;
                dQ[p].set(0, r, dQ_L);
                dX.set(r, p, dX_L);
                dU[p].set(0, r, dU_L);
                if (Xr > 0) {
                    dC[p].set(0, r, (dQ_L * Xr - Qr * dX_L) / (Xr * Xr));
                }

                // Z(s) parameter: d log G_m(n)/dZ_s = G_m(n-1_s)/G_m(n)
                p = pZ[s];
                double dX_Z = Xr * (nc.xput(1, Nr, s) - X.get(0, s));
                double dQ_Z = Qr * (nc.xput(2, Nr, s) - X.get(0, s));
                double dU_Z = dX_Z * D[r];
                dQ[p].set(0, r, dQ_Z);
                dX.set(r, p, dX_Z);
                dU[p].set(0, r, dU_Z);
                if (Xr > 0) {
                    dC[p].set(0, r, (dQ_Z * Xr - Qr * dX_Z) / (Xr * Xr));
                }
            }
        }

        return new Ret.pfqnSens(X, Q, U, C, dX, dQ, dU, dC,
                paramType, paramStation, paramClass);
    }

    /** Replicated-model normalizing-constant moments via memoized CoMoM. */
    private static final class ReplNC {
        private final double[] D;
        private final double[] Z;
        private final int R;
        private final Map<String, Double> cache = new HashMap<String, Double>();

        ReplNC(double[] D, double[] Z, int R) {
            this.D = D;
            this.Z = Z;
            this.R = R;
        }

        /** log normalizing constant of the m-replica model at population n. */
        double lgm(int m, int[] n) {
            int k = 0;
            for (int r = 0; r < R; r++) {
                if (n[r] > 0) {
                    k++;
                }
            }
            if (k == 0) {
                return 0.0;
            }
            StringBuilder sb = new StringBuilder();
            sb.append(m).append('|');
            for (int r = 0; r < R; r++) {
                sb.append(n[r]).append(',');
            }
            String key = sb.toString();
            Double cached = cache.get(key);
            if (cached != null) {
                return cached;
            }
            // strip zero-population classes (they leave the NC unchanged)
            Matrix Ls = new Matrix(1, k);
            Matrix Ns = new Matrix(1, k);
            Matrix Zs = new Matrix(1, k);
            int j = 0;
            for (int r = 0; r < R; r++) {
                if (n[r] > 0) {
                    Ls.set(0, j, D[r]);
                    Ns.set(0, j, n[r]);
                    Zs.set(0, j, Z[r]);
                    j++;
                }
            }
            double lg = Pfqn_comomrm.pfqn_comomrm(Ls, Ns, Zs, m, GlobalConstants.FineTol).lG;
            cache.put(key, lg);
            return lg;
        }

        /** mean class-s queue at population n (single station, m=1 model). */
        double qmean(int[] n, int s) {
            if (n[s] < 1) {
                return 0.0;
            }
            int[] nm = n.clone();
            nm[s] -= 1;
            return D[s] * FastMath.exp(lgm(2, nm) - lgm(1, n));
        }

        /** class-s throughput G_m(n-1_s)/G_m(n) in the m-replica model. */
        double xput(int m, int[] n, int s) {
            if (n[s] < 1) {
                return 0.0;
            }
            int[] nm = n.clone();
            nm[s] -= 1;
            return FastMath.exp(lgm(m, nm) - lgm(m, n));
        }

        /** Q^{+1}_{1,s}(n): class-s queue at one replica of the doubled station. */
        double qplus(int[] n, int s) {
            if (n[s] < 1) {
                return 0.0;
            }
            int[] nm = n.clone();
            nm[s] -= 1;
            return D[s] * FastMath.exp(lgm(3, nm) - lgm(2, n));
        }
    }

    // =====================================================================
    // Differentiated-MVA kernel (general model)
    // =====================================================================
    private static Ret.pfqnSens sensMva(Matrix L, Matrix N, Matrix Z, Matrix mi) {
        int M = L.getNumRows();
        int R = L.getNumCols();

        // ---- parameter list: all L(i,r), then all Z(r) ------------------
        int P = M * R + R;
        int[] paramType = new int[P];
        int[] paramStation = new int[P];
        int[] paramClass = new int[P];
        int[][] pL = new int[M][R];
        int[] pZ = new int[R];
        int p = 0;
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                paramType[p] = 0;
                paramStation[p] = i;
                paramClass[p] = r;
                pL[i][r] = p;
                p++;
            }
        }
        for (int r = 0; r < R; r++) {
            paramType[p] = 1;
            paramStation[p] = -1;
            paramClass[p] = r;
            pZ[r] = p;
            p++;
        }

        Matrix X = new Matrix(1, R);
        Matrix Q = new Matrix(M, R);
        Matrix U = new Matrix(M, R);
        Matrix C = new Matrix(M, R);
        Matrix dX = new Matrix(R, P);
        Matrix[] dQ = zeros(P, M, R);
        Matrix[] dU = zeros(P, M, R);
        Matrix[] dC = zeros(P, M, R);

        if (!N.any()) {
            return new Ret.pfqnSens(X, Q, U, C, dX, dQ, dU, dC,
                    paramType, paramStation, paramClass);
        }

        // ---- population-lattice odometer, identical to pfqn_mva ---------
        Matrix prods = new Matrix(1, Math.max(0, R - 1));
        for (int w = 0; w < R - 1; w++) {
            double acc = 1.0;
            for (int i = 0; i < R - (w + 2) + 1; i++) {
                acc *= (1.0 + N.get(0, w + 1 + i));
            }
            prods.set(0, w, acc);
        }
        int firstNonEmpty = R - 1;
        while (firstNonEmpty >= 0 && N.get(0, firstNonEmpty) == 0.0) {
            firstNonEmpty--;
        }
        double totpop = 1.0;
        for (int r = 0; r < R; r++) {
            totpop *= (N.get(0, r) + 1.0);
        }
        int TP = (int) totpop;
        double ctr = totpop;
        Matrix Qtot = new Matrix(TP, M);
        double[][][] Qtotd = new double[P][TP][M];
        int currentpop = 1;

        Matrix n = new Matrix(1, R);
        n.set(0, firstNonEmpty, 1);

        double[] CNtotd = new double[P];
        double[] Xd = new double[P];
        while (ctr > 0) {
            int s = 0;
            while (s < R) {
                int pos = 0;
                if (n.get(0, s) > 0) {
                    n.set(0, s, n.get(0, s) - 1);
                    pos = (int) n.get(0, R - 1);
                    int w = 0;
                    while (w < R - 1) {
                        pos = (int) (pos + n.get(0, w) * prods.get(0, w));
                        w++;
                    }
                    n.set(0, s, n.get(0, s) + 1);
                }
                double CNtot = 0.0;
                for (int j = 0; j < P; j++) {
                    CNtotd[j] = 0.0;
                }
                double[][] CdAll = new double[M][P];
                for (int i = 0; i < M; i++) {
                    double base = mi.get(0, i) + Qtot.get(pos, i);
                    double Lis = L.get(i, s);
                    C.set(i, s, Lis * base);
                    for (int j = 0; j < P; j++) {
                        double v = Lis * Qtotd[j][pos][i];
                        CdAll[i][j] = v;
                        CNtotd[j] += v;
                    }
                    CdAll[i][pL[i][s]] += base;
                    CNtotd[pL[i][s]] += base;
                    CNtot += C.get(i, s);
                }
                double den = Z.get(0, s) + CNtot;
                double ns = n.get(0, s);
                X.set(0, s, den > 0 ? ns / den : 0.0);
                for (int j = 0; j < P; j++) {
                    Xd[j] = den > 0 ? (-ns * CNtotd[j] / (den * den)) : 0.0;
                }
                if (den > 0) {
                    Xd[pZ[s]] -= ns / (den * den);
                }
                for (int j = 0; j < P; j++) {
                    dX.set(s, j, Xd[j]);
                }
                double xs = X.get(0, s);
                for (int i = 0; i < M; i++) {
                    double cis = C.get(i, s);
                    Q.set(i, s, xs * cis);
                    for (int j = 0; j < P; j++) {
                        double qd = Xd[j] * cis + xs * CdAll[i][j];
                        dQ[j].set(i, s, qd);
                        dC[j].set(i, s, CdAll[i][j]);
                        Qtotd[j][currentpop][i] += qd;
                    }
                    Qtot.set(currentpop, i, Qtot.get(currentpop, i) + xs * cis);
                }
                s++;
            }
            s = R - 1;
            while ((s >= 0 && (n.get(0, s) == N.get(0, s))) || s > firstNonEmpty) {
                s--;
            }
            if (s == -1) {
                break;
            }
            n.set(0, s, n.get(0, s) + 1);
            s++;
            while (s < R) {
                n.set(0, s, 0);
                s++;
            }
            ctr--;
            currentpop++;
        }

        // ---- utilization and its derivatives ---------------------------
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                U.set(i, r, X.get(0, r) * L.get(i, r));
                for (int j = 0; j < P; j++) {
                    double ud = dX.get(r, j) * L.get(i, r);
                    dU[j].set(i, r, ud);
                }
                dU[pL[i][r]].set(i, r, dU[pL[i][r]].get(i, r) + X.get(0, r));
            }
        }

        return new Ret.pfqnSens(X, Q, U, C, dX, dQ, dU, dC,
                paramType, paramStation, paramClass);
    }

    // =====================================================================
    // Exact queue-length second moments. Both sources below rest on the same
    // product-form identity Cov[n_{i,r},n_{j,s}] = D_{j,s} dQ_{i,r}/dD_{j,s},
    // which follows from L_{j,s} d/dL_{j,s}(G Q_{i,r}) = G E[n_{i,r} n_{j,s}]:
    //  - same-station blocks (i==j) come from Pfqn_sens_mva, which evaluates them
    //    by the exact MVA-like recursion of de Souza e Silva and Muntz (1988),
    //    Corollary 1. That recursion is self-contained at a station and costs
    //    about 1/(M+1) of the Jacobian pass above, so sourcing the variances
    //    from it is essentially free and independent of the differentiated-MVA
    //    algorithm.
    //  - cross-station blocks (i!=j) are read off the Jacobian. They are not
    //    self-contained in the recursion, which for them would couple every
    //    station pair and cost as much as the Jacobian itself.
    private static void attachSecondMoments(Ret.pfqnSens s, Matrix L, Matrix N, Matrix Z, Matrix mi) {
        int M = L.getNumRows();
        int R = L.getNumCols();
        Ret.pfqnSensMva mom = Pfqn_sens_mva.pfqn_sens_mva(L, N, Z, mi);
        int[][] pLidx = new int[M][R];
        boolean[][] pLset = new boolean[M][R];
        for (int p = 0; p < s.paramType.length; p++) {
            if (s.paramType[p] == 0) {
                pLidx[s.paramStation[p]][s.paramClass[p]] = p;
                pLset[s.paramStation[p]][s.paramClass[p]] = true;
            }
        }
        Matrix[][] QCov = new Matrix[M][R];
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                Matrix cov = new Matrix(M, R);
                for (int j = 0; j < M; j++) {
                    for (int sc = 0; sc < R; sc++) {
                        if (i == j) {
                            cov.set(j, sc, mom.QCov[i].get(r, sc));
                        } else if (pLset[j][sc]) {
                            int p = pLidx[j][sc];
                            cov.set(j, sc, L.get(j, sc) * s.dQ[p].get(i, r));
                        }
                    }
                }
                QCov[i][r] = cov;
            }
        }
        s.QCov = QCov;
        s.QVar = mom.QVar;
        s.QTotVar = mom.QTotVar;
        s.QCovAsym = mom.QCovAsym;
    }

    private static Matrix[] zeros(int P, int M, int R) {
        Matrix[] arr = new Matrix[P];
        for (int j = 0; j < P; j++) {
            arr[j] = new Matrix(M, R);
        }
        return arr;
    }
}
