/**
 * @file de Souza e Silva-Lavenberg-Muntz Clustering Approximation (CA)
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.mva;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import jline.io.Ret;
import jline.util.matrix.Matrix;

/**
 * de Souza e Silva-Lavenberg-Muntz Clustering Approximation (CA).
 *
 * <p>E. de Souza e Silva, S. S. Lavenberg, R. R. Muntz, "A clustering
 * approximation technique for queueing network models with a large number of
 * chains", IEEE Trans. Computers C-35(5), 1986. The network is covered by
 * subnetworks whose union is the whole network but which need not be disjoint.
 * Every class visiting a subnetwork S is either LOCAL to S, and is then solved
 * inside it, or FOREIGN, and is then seen only through the utilization it
 * leaves behind. Each subnetwork is solved by an ordinary approximate MVA
 * algorithm with two replacements: the complement of S is collapsed into a
 * per-class delay P_c and the foreign classes into a per-centre utilization
 * U_k,</p>
 *
 * <pre>
 * X_c(N) = N_c / (sum_{k in S} R_ck(N) + Z_c + P_c),
 * Q_k(N) = [sum_{c in LC(S)} R_ck(N) X_c(N) + U_k] / (1 - U_k).
 * </pre>
 *
 * <p>Choosing the PE algorithm for every subnetwork reproduces global PE
 * exactly, so the useful setting is Linearizer inside, PE outside: the cost
 * then sits between {@link Pfqn_bs} and {@link Pfqn_linearizer}, which is the
 * point of the method.</p>
 *
 * <p>When no decomposition is supplied the criterion of the paper is applied
 * automatically: the cheap PAMB estimate of the centre utilizations is taken,
 * every class is attached to the centre where it loads the most, classes
 * sharing that centre form one cluster, and the subnetwork of a cluster is the
 * set of centres its classes visit.</p>
 */
public final class Pfqn_clust {
    private Pfqn_clust() {}

    public static Ret.pfqnAMVA pfqn_clust(Matrix L, Matrix N, Matrix Z) {
        return pfqn_clust(L, N, Z, null, null, "lin", 1e-6, 1000);
    }

    public static Ret.pfqnAMVA pfqn_clust(Matrix L, Matrix N, Matrix Z, double tol, int maxiter) {
        return pfqn_clust(L, N, Z, null, null, "lin", tol, maxiter);
    }

    public static Ret.pfqnAMVA pfqn_clust(Matrix L, Matrix N, Matrix Zin,
                                          List<int[]> subnets, List<int[]> localclasses,
                                          String inner, double tol, int maxiter) {
        int M = L.getNumRows();
        int R = L.getNumCols();
        Matrix Z = (Zin == null || Zin.isEmpty()) ? new Matrix(1, R) : Zin;
        String innerAlg = (inner == null) ? "lin" : inner;

        Ret.pfqnAMVA seed = Pfqn_pam.pfqn_pam(L, N, Z, Pfqn_pam.PAMB);
        Matrix XN = seed.X.copy();
        Matrix QN = seed.Q.copy();

        List<int[]> nets = subnets;
        List<int[]> locals = localclasses;
        if (nets == null || locals == null || nets.isEmpty() || locals.isEmpty()) {
            int[] bottleneck = new int[R];
            for (int r = 0; r < R; r++) {
                double best = -1.0;
                int arg = 0;
                for (int i = 0; i < M; i++) {
                    double u = L.get(i, r) * XN.get(r);
                    if (L.get(i, r) > 0 && u > best) {
                        best = u;
                        arg = i;
                    }
                }
                bottleneck[r] = arg;
            }
            // TreeMap, not insertion order: the outer sweep is Gauss-Seidel, so
            // the cluster order is part of the answer. MATLAB takes
            // unique(bottleneck), which is sorted by centre; first-encounter
            // order moved the 7th digit.
            Map<Integer, List<Integer>> byCentre = new TreeMap<Integer, List<Integer>>();
            for (int r = 0; r < R; r++) {
                List<Integer> lst = byCentre.get(bottleneck[r]);
                if (lst == null) {
                    lst = new ArrayList<Integer>();
                    byCentre.put(bottleneck[r], lst);
                }
                lst.add(r);
            }
            nets = new ArrayList<int[]>();
            locals = new ArrayList<int[]>();
            boolean[] covered = new boolean[M];
            for (Map.Entry<Integer, List<Integer>> e : byCentre.entrySet()) {
                List<Integer> cls = e.getValue();
                int[] lc = new int[cls.size()];
                for (int a = 0; a < cls.size(); a++) {
                    lc[a] = cls.get(a);
                }
                List<Integer> st = new ArrayList<Integer>();
                for (int i = 0; i < M; i++) {
                    for (int a = 0; a < lc.length; a++) {
                        if (L.get(i, lc[a]) > 0) {
                            st.add(i);
                            break;
                        }
                    }
                }
                if (st.isEmpty()) {
                    st.add(e.getKey());
                }
                int[] s = new int[st.size()];
                for (int a = 0; a < st.size(); a++) {
                    s[a] = st.get(a);
                    covered[s[a]] = true;
                }
                nets.add(s);
                locals.add(lc);
            }
            List<Integer> missing = new ArrayList<Integer>();
            for (int i = 0; i < M; i++) {
                if (!covered[i]) {
                    missing.add(i);
                }
            }
            if (!missing.isEmpty()) {
                int[] s = new int[missing.size()];
                for (int a = 0; a < missing.size(); a++) {
                    s[a] = missing.get(a);
                }
                nets.add(s);
                locals.add(new int[0]);
            }
        }
        int G = nets.size();
        int[] owner = new int[R];
        for (int r = 0; r < R; r++) {
            owner[r] = -1;
        }
        for (int g = 0; g < G; g++) {
            for (int a = 0; a < locals.get(g).length; a++) {
                owner[locals.get(g)[a]] = g;
            }
        }

        int it = 1;
        while (it <= maxiter) {
            Matrix QN_1 = QN.copy();
            double[] Qk = new double[M];
            for (int i = 0; i < M; i++) {
                double acc = 0.0;
                for (int r = 0; r < R; r++) {
                    acc += QN.get(i, r);
                }
                Qk[i] = acc;
            }
            for (int g = 0; g < G; g++) {
                int[] S = nets.get(g);
                int[] LC = locals.get(g);
                if (LC.length == 0 || S.length == 0) {
                    continue;
                }
                boolean[] inS = new boolean[M];
                for (int a = 0; a < S.length; a++) {
                    inS[S[a]] = true;
                }
                boolean[] isLocal = new boolean[R];
                for (int a = 0; a < LC.length; a++) {
                    isLocal[LC[a]] = true;
                }
                // foreign classes: visit S but are not local to it
                boolean[] isForeign = new boolean[R];
                for (int r = 0; r < R; r++) {
                    if (isLocal[r]) continue;
                    for (int a = 0; a < S.length; a++) {
                        if (L.get(S[a], r) > 0) {
                            isForeign[r] = true;
                            break;
                        }
                    }
                }
                // per-class delay in the complement of S
                Matrix Zeff = new Matrix(1, LC.length);
                for (int a = 0; a < LC.length; a++) {
                    int c = LC[a];
                    double p = 0.0;
                    if (N.get(c) > 0) {
                        for (int i = 0; i < M; i++) {
                            if (!inS[i]) {
                                p += L.get(i, c) * (1 + Qk[i]) / (1 + L.get(i, c) * XN.get(c) / N.get(c));
                            }
                        }
                    }
                    Zeff.set(a, Z.get(c) + p);
                }
                // utilization left in S by the foreign classes
                double[] Uk = new double[S.length];
                for (int b = 0; b < S.length; b++) {
                    int i = S[b];
                    double u = 0.0;
                    for (int r = 0; r < R; r++) {
                        if (isForeign[r] && N.get(r) > 0) {
                            u += L.get(i, r) * XN.get(r) / (1 + L.get(i, r) * XN.get(r) / N.get(r));
                        }
                    }
                    Uk[b] = Math.min(u, 1 - 1e-8);
                }
                Matrix Lsub = new Matrix(S.length, LC.length);
                Matrix Nsub = new Matrix(1, LC.length);
                for (int b = 0; b < S.length; b++) {
                    for (int a = 0; a < LC.length; a++) {
                        Lsub.set(b, a, L.get(S[b], LC[a]));
                    }
                }
                for (int a = 0; a < LC.length; a++) {
                    Nsub.set(a, N.get(LC[a]));
                }
                Matrix[] sub = subnetSolve(Lsub, Nsub, Zeff, Uk, innerAlg, tol, maxiter);
                Matrix Xs = sub[0];
                Matrix Qs = sub[1];
                for (int a = 0; a < LC.length; a++) {
                    XN.set(LC[a], Xs.get(a));
                    for (int b = 0; b < S.length; b++) {
                        QN.set(S[b], LC[a], Qs.get(b, a));
                    }
                }
                // the local classes still hold jobs outside S
                for (int a = 0; a < LC.length; a++) {
                    int c = LC[a];
                    double nc = Math.max(N.get(c), Double.MIN_NORMAL);
                    for (int i = 0; i < M; i++) {
                        if (!inS[i]) {
                            QN.set(i, c, XN.get(c) * L.get(i, c) * (1 + Qk[i])
                                    / (1 + L.get(i, c) * XN.get(c) / nc));
                        }
                    }
                }
            }
            for (int r = 0; r < R; r++) {
                if (owner[r] < 0) {
                    for (int i = 0; i < M; i++) {
                        QN.set(i, r, XN.get(r) * L.get(i, r));
                    }
                }
            }
            double maxdiff = 0.0;
            for (int i = 0; i < M; i++) {
                for (int r = 0; r < R; r++) {
                    maxdiff = Math.max(maxdiff, Math.abs(QN.get(i, r) - QN_1.get(i, r)));
                }
            }
            if (maxdiff < tol) {
                break;
            }
            it++;
        }

        Matrix UN = new Matrix(M, R);
        Matrix RN = new Matrix(M, R);
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                UN.set(i, r, XN.get(r) * L.get(i, r));
                double v = (N.get(r) == 0.0 || XN.get(r) == 0.0) ? 0.0 : QN.get(i, r) / XN.get(r);
                RN.set(i, r, Double.isFinite(v) ? v : 0.0);
            }
        }
        return new Ret.pfqnAMVA(QN, UN, RN, null, RN, XN, it);
    }

    /**
     * Approximate MVA restricted to the local classes of one subnetwork, with the
     * complement folded into Z and the foreign classes into Uk.
     */
    private static Matrix[] subnetSolve(Matrix L, Matrix N, Matrix Z, double[] Uk,
                                        String inner, double tol, int maxiter) {
        int M = L.getNumRows();
        int R = L.getNumCols();
        Matrix X = new Matrix(1, R);
        Matrix Q = new Matrix(M, R);
        if (M == 0 || R == 0) {
            return new Matrix[]{X, Q};
        }
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                Q.set(i, r, N.get(r) / M);
            }
        }
        if ("bs".equals(inner)) {
            for (int it = 0; it < maxiter; it++) {
                Matrix Qold = Q.copy();
                for (int r = 0; r < R; r++) {
                    if (N.get(r) <= 0) {
                        continue;
                    }
                    double[] W = new double[M];
                    double wsum = 0.0;
                    for (int i = 0; i < M; i++) {
                        double loc = 0.0;
                        for (int s = 0; s < R; s++) {
                            loc += Q.get(i, s);
                        }
                        loc -= Q.get(i, r) / N.get(r);
                        double A = (loc + Uk[i]) / (1 - Uk[i]);
                        W[i] = L.get(i, r) * (1 + A);
                        wsum += W[i];
                    }
                    X.set(r, N.get(r) / (Z.get(r) + wsum));
                    for (int i = 0; i < M; i++) {
                        Q.set(i, r, X.get(r) * W[i]);
                    }
                }
                if (maxAbsDiff(Q, Qold) < tol) {
                    break;
                }
            }
            return new Matrix[]{X, Q};
        }
        // Linearizer inside the subnetwork
        Matrix[] Qs = new Matrix[R + 1];
        for (int s = 0; s <= R; s++) {
            Qs[s] = Q.copy();
        }
        double[][][] Delta = new double[M][R][R];
        for (int pass = 0; pass < 3; pass++) {
            for (int s = 0; s <= R; s++) {
                Matrix Ns = N.copy();
                if (s > 0) {
                    Ns.set(s - 1, N.get(s - 1) - 1);
                }
                Qs[s] = subnetCore(L, Ns, Z, Uk, Qs[s], Delta, tol, maxiter);
            }
            for (int i = 0; i < M; i++) {
                for (int r = 0; r < R; r++) {
                    for (int s = 1; s <= R; s++) {
                        double ns = N.get(r) - (r == s - 1 ? 1.0 : 0.0);
                        if (N.get(r) > 0 && ns > 0) {
                            Delta[i][r][s - 1] = Qs[s].get(i, r) / ns - Qs[0].get(i, r) / N.get(r);
                        } else if (N.get(r) > 0) {
                            Delta[i][r][s - 1] = -Qs[0].get(i, r) / N.get(r);
                        } else {
                            Delta[i][r][s - 1] = 0.0;
                        }
                    }
                }
            }
        }
        Qs[0] = subnetCore(L, N, Z, Uk, Qs[0], Delta, tol, maxiter);
        return subnetForward(L, N, Z, Uk, Qs[0], Delta);
    }

    private static Matrix subnetCore(Matrix L, Matrix N, Matrix Z, double[] Uk,
                                     Matrix Q, double[][][] Delta, double tol, int maxiter) {
        Matrix cur = Q;
        for (int it = 0; it < maxiter; it++) {
            Matrix old = cur;
            cur = subnetForward(L, N, Z, Uk, cur, Delta)[1];
            if (maxAbsDiff(cur, old) < tol) {
                break;
            }
        }
        return cur;
    }

    private static Matrix[] subnetForward(Matrix L, Matrix N, Matrix Z, double[] Uk,
                                          Matrix Q, double[][][] Delta) {
        int M = L.getNumRows();
        int R = L.getNumCols();
        Matrix Qout = new Matrix(M, R);
        Matrix X = new Matrix(1, R);
        for (int r = 0; r < R; r++) {
            if (N.get(r) <= 0) {
                continue;
            }
            double[] W = new double[M];
            double wsum = 0.0;
            for (int i = 0; i < M; i++) {
                double qm = 0.0;
                for (int s = 0; s < R; s++) {
                    double ns = N.get(s) - (s == r ? 1.0 : 0.0);
                    if (N.get(s) > 0 && ns > 0) {
                        qm += ns * (Q.get(i, s) / N.get(s) + Delta[i][s][r]);
                    }
                }
                if (qm < 0) {
                    qm = 0;
                }
                double A = (qm + Uk[i]) / (1 - Uk[i]);
                W[i] = L.get(i, r) * (1 + A);
                wsum += W[i];
            }
            X.set(r, N.get(r) / (Z.get(r) + wsum));
            for (int i = 0; i < M; i++) {
                Qout.set(i, r, X.get(r) * W[i]);
            }
        }
        return new Matrix[]{X, Qout};
    }

    private static double maxAbsDiff(Matrix a, Matrix b) {
        double m = 0.0;
        for (int i = 0; i < a.getNumRows(); i++) {
            for (int j = 0; j < a.getNumCols(); j++) {
                m = Math.max(m, Math.abs(a.get(i, j) - b.get(i, j)));
            }
        }
        return m;
    }
}
