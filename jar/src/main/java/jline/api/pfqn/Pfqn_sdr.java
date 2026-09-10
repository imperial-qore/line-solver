/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.pfqn;

import java.util.ArrayList;
import java.util.List;

import jline.lang.StateDepRouting;
import jline.util.matrix.Matrix;

/**
 * Product-form state-dependent routing of Krzesinski (1987), "Multiclass
 * Queueing Networks with State-Dependent Routing", Performance Evaluation
 * 7(2):125-143, the multiclass generalization of Towsley (1980), J. ACM
 * 27(2):323-337.
 *
 * <p>Branch index 1 denotes the complement M-V and is unused; the SDR branches
 * are numbered 2..B, that is indices 1..B-1 of the zero-based arrays, following
 * the paper's own indexing so that d_tb transcribes straight from the text.</p>
 *
 * @see jline.lang.StateDepRouting
 */
public class Pfqn_sdr {

    /** Derived coefficients of an SDR structure, eqs. (11)-(14). */
    public static class Coeff {
        /** Number of levels of subnetwork nesting. */
        public int T;
        /** Number of branch indices, including the unused complement index 0. */
        public int B;
        /** level[b], the unique t with B_b in V_t - V_{t+1}. */
        public int[] level;
        /** Coefficients C_t. */
        public double[] C;
        /** Coefficients d_tb. */
        public double[][] d;
        /** inA[t] holds the branch indices b with level[b] &gt;= t, the set A_t. */
        public int[][] inA;
        /** D_tt = sum over A_t of d_tb, eq. (14); index 0 unused. */
        public double[] Dtt;
        /** D_{t-1,t} = sum over A_t of d_{t-1,b}, eq. (14); indices 0 and 1 unused. */
        public double[] Dprev;
        /** Largest branch population with delta nonnegative; infinite when C_t &gt;= 0. */
        public double[] mmax;
        /** Largest subnetwork population with omega nonnegative; infinite when C_t &gt;= 0. */
        public double[] vmax;
        /** The validated structure. */
        public StateDepRouting sdr;
    }

    /**
     * Validates an SDR structure and returns its derived coefficients.
     *
     * <p>The population bounds are consequences of the coefficients, not
     * independent inputs: with C_t negative the routing enforces
     * m_b &lt;= d_tb/(-C_t) and v_t &lt;= D_tt/(-C_t) by itself.</p>
     *
     * @param sdr the structure to validate
     * @return its derived coefficients
     */
    public static Coeff pfqn_sdrcoeff(StateDepRouting sdr) {
        if (sdr == null) {
            throw new RuntimeException("The SDR structure is null.");
        }
        int B = sdr.branch == null ? 0 : sdr.branch.length;
        int T = sdr.C == null ? 0 : sdr.C.length;
        if (B < 2) {
            throw new RuntimeException("An SDR structure must declare at least one branch "
                    + "(branch indices start at 2).");
        }
        if (T < 1) {
            throw new RuntimeException("An SDR structure must declare at least one level of subnetwork nesting.");
        }
        if (sdr.level == null || sdr.level.length != B) {
            throw new RuntimeException("SDR level must have one entry per branch index.");
        }
        for (int b = 1; b < B; b++) {
            if (sdr.level[b] < 1 || sdr.level[b] > T) {
                throw new RuntimeException("SDR branch levels must be integers in 1.." + T + ".");
            }
        }
        if (sdr.d == null || sdr.d.length < T) {
            throw new RuntimeException("SDR coefficient matrix d must be at least " + T + "x" + B + ".");
        }
        for (int t = 0; t < T; t++) {
            if (sdr.d[t] == null || sdr.d[t].length < B) {
                throw new RuntimeException("SDR coefficient matrix d must be at least " + T + "x" + B + ".");
            }
        }

        // The nesting V_1 > ... > V_T must be strict, else two hierarchically
        // adjacent subnetworks coincide and the ratio of eq. (10) is not the paper's
        for (int t = 1; t <= T; t++) {
            boolean found = false;
            for (int b = 1; b < B; b++) {
                if (sdr.level[b] == t) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                throw new RuntimeException("SDR level " + t + " carries no branch: the subnetwork nesting must be strict.");
            }
        }

        List<Integer> seen = new ArrayList<Integer>();
        for (int b = 1; b < B; b++) {
            int[] sb = sdr.branch[b];
            if (sb == null || sb.length == 0) {
                throw new RuntimeException("SDR branch " + (b + 1) + " is empty.");
            }
            boolean hasEntry = false;
            boolean hasDeparture = false;
            for (int k = 0; k < sb.length; k++) {
                if (seen.contains(Integer.valueOf(sb[k]))) {
                    throw new RuntimeException("SDR branches must be mutually disjoint.");
                }
                if (sb[k] == sdr.entry || sb[k] == sdr.departure) {
                    throw new RuntimeException("The entry and departure centers of Q(V,V) must not belong to any branch.");
                }
                if (sb[k] == sdr.entryOf[b]) {
                    hasEntry = true;
                }
                if (sb[k] == sdr.departureOf[b]) {
                    hasDeparture = true;
                }
                seen.add(Integer.valueOf(sb[k]));
            }
            if (!hasEntry || !hasDeparture) {
                throw new RuntimeException("The entry and departure centers of SDR branch " + (b + 1)
                        + " must belong to that branch.");
            }
        }

        Coeff c = new Coeff();
        c.T = T;
        c.B = B;
        c.level = sdr.level.clone();
        c.C = sdr.C.clone();
        c.d = sdr.d;
        c.sdr = sdr;
        c.inA = new int[T + 1][];
        c.Dtt = new double[T + 1];
        c.Dprev = new double[T + 1];
        for (int t = 1; t <= T; t++) {
            List<Integer> at = new ArrayList<Integer>();
            for (int b = 1; b < B; b++) {
                if (sdr.level[b] >= t) {
                    at.add(Integer.valueOf(b));
                }
            }
            int[] arr = new int[at.size()];
            for (int k = 0; k < at.size(); k++) {
                arr[k] = at.get(k).intValue();
                c.Dtt[t] += sdr.d[t - 1][arr[k]];
                if (t > 1) {
                    c.Dprev[t] += sdr.d[t - 2][arr[k]];
                }
            }
            c.inA[t] = arr;
        }
        c.mmax = new double[B];
        for (int b = 1; b < B; b++) {
            int t = sdr.level[b];
            c.mmax[b] = (c.C[t - 1] < 0) ? Math.floor(sdr.d[t - 1][b] / (-c.C[t - 1])) : Double.POSITIVE_INFINITY;
        }
        c.vmax = new double[T + 1];
        for (int t = 1; t <= T; t++) {
            c.vmax[t] = (c.C[t - 1] < 0) ? Math.floor(c.Dtt[t] / (-c.C[t - 1])) : Double.POSITIVE_INFINITY;
            if (t > 1 && c.C[t - 2] < 0) {
                c.vmax[t] = Math.min(c.vmax[t], Math.floor(c.Dprev[t] / (-c.C[t - 2])));
            }
        }
        return c;
    }

    /**
     * SDR routing probabilities of eq. (10).
     *
     * <p>Entry b of the returned array is the probability of proceeding from the
     * entry center e of Q(V,V) to the entry center of branch b; entry 0 is zero
     * because branch index 1 denotes the complement M-V. The residual mass
     * 1 - sum is the probability of proceeding directly to the departure center
     * d, that is of being denied entry into Q(V,V) and returned to e, which is
     * the busy form of waiting of Section 2.5.</p>
     *
     * <p>These probabilities are chain independent: they read the total branch
     * and subnetwork populations, not the per-chain ones. The chain-dependent
     * form of eq. (1) has no published product form and is not implemented. A
     * branch population beyond the bound that SDR enforces itself is
     * unreachable, and the probability returned there is zero.</p>
     *
     * @param c derived coefficients from {@link #pfqn_sdrcoeff}
     * @param n per-center total populations, indexed as the structure is
     * @return the probabilities to each branch entry, index 0 unused
     */
    public static double[] pfqn_sdrprob(Coeff c, double[] n) {
        double[] m = new double[c.B];
        for (int b = 1; b < c.B; b++) {
            double acc = 0;
            for (int k = 0; k < c.sdr.branch[b].length; k++) {
                acc += n[c.sdr.branch[b][k]];
            }
            m[b] = acc;
        }
        double[] v = new double[c.T + 1];
        for (int t = 1; t <= c.T; t++) {
            double acc = 0;
            for (int k = 0; k < c.inA[t].length; k++) {
                acc += m[c.inA[t][k]];
            }
            v[t] = acc;
        }
        double[] om = new double[c.T + 1];
        double[] omprev = new double[c.T + 1];
        for (int s = 1; s <= c.T; s++) {
            om[s] = c.C[s - 1] * v[s] + c.Dtt[s];
            omprev[s] = (s > 1) ? (c.C[s - 2] * v[s] + c.Dprev[s]) : 1.0;
        }
        double[] P = new double[c.B];
        for (int b = 1; b < c.B; b++) {
            int t = c.level[b];
            boolean closed = false;
            for (int s = 1; s <= t; s++) {
                if (om[s] <= 0) {
                    closed = true; // eq. (10): the branch is closed to new arrivals
                    break;
                }
            }
            if (closed) {
                continue;
            }
            double delta = c.C[t - 1] * m[b] + c.d[t - 1][b];
            if (delta <= 0) {
                continue;
            }
            double ratio = 1.0;
            for (int s = 1; s <= t; s++) {
                ratio *= omprev[s] / om[s];
            }
            P[b] = delta * ratio;
        }
        return P;
    }

    /**
     * Probability of being denied entry into Q(V,V) and routed straight to the
     * departure center, given the branch probabilities of {@link #pfqn_sdrprob}.
     *
     * @param P the branch probabilities
     * @return one minus their sum
     */
    public static double pfqn_sdrped(double[] P) {
        double s = 0;
        for (int b = 0; b < P.length; b++) {
            s += P[b];
        }
        return 1.0 - s;
    }

    /** Mean performance measures returned by {@link #pfqn_sdr}. */
    public static class Result {
        /** Mean queue lengths, centers by chains. */
        public Matrix Q;
        /** Per-center chain throughputs. */
        public Matrix X;
        /** Mean number in service, X elementwise times S. */
        public Matrix U;
        /** Mean response times at the center, Q elementwise over X. */
        public Matrix R;
        /** Logarithm of the normalizing constant. */
        public double lG;
    }

    /**
     * Exact product form of eq. (16), by summation over the reachable state space.
     *
     * <p>P(n) is G^-1 times the product over centers of f_i(n_i), the product
     * over levels of Omega_{t-1,t}(v_t)/Omega_tt(v_t), and the product over
     * branches of Delta_tb(m_b), with f_i(n_i) = [n_i!/beta_i(n_i)] times the
     * product over chains of gamma_ij^n_ij/n_ij! and gamma_ij = xi_ij/mu_ij.</p>
     *
     * <p>This is general in the branch topology: a branch may hold several
     * interconnected centers. Only the paper's Section 4 MVA and convolution
     * algorithm is restricted to single-center branches.</p>
     *
     * <p>S and xi are required separately rather than as their product because
     * under SDR the xi are not visit ratios, so the per-center throughputs
     * cannot be recovered from the demands alone.</p>
     *
     * @param S mean service times, centers by chains
     * @param xi coefficients of Section 3.2, centers by chains
     * @param N chain populations
     * @param sdr the routing structure, in center indices
     * @param alpha load-dependent rate scalings, alpha.get(i,k-1) = alpha_i(k); null for fixed rate
     * @return the mean performance measures
     */
    public static Result pfqn_sdr(Matrix S, Matrix xi, Matrix N, StateDepRouting sdr, Matrix alpha) {
        int M = S.getNumRows();
        int J = S.getNumCols();
        if (xi.getNumRows() != M || xi.getNumCols() != J) {
            throw new RuntimeException("S and xi must have the same size.");
        }
        int[] Nv = new int[J];
        int Ntot = 0;
        for (int j = 0; j < J; j++) {
            Nv[j] = (int) N.get(j);
            if (Nv[j] < 0 || Double.isInfinite(N.get(j))) {
                throw new RuntimeException("State-dependent routing is defined for closed networks only.");
            }
            Ntot += Nv[j];
        }
        Coeff c = pfqn_sdrcoeff(sdr);

        double[][] alp = new double[M][Math.max(1, Ntot)];
        for (int i = 0; i < M; i++) {
            for (int k = 0; k < alp[i].length; k++) {
                alp[i][k] = (alpha != null && k < alpha.getNumCols() && i < alpha.getNumRows())
                        ? alpha.get(i, k) : 1.0;
            }
        }
        double[][] logbeta = new double[M][Ntot + 1];
        for (int i = 0; i < M; i++) {
            for (int k = 1; k <= Ntot; k++) {
                logbeta[i][k] = logbeta[i][k - 1] + Math.log(alp[i][k - 1]);
            }
        }
        double[][] gamma = new double[M][J];
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < J; j++) {
                gamma[i][j] = xi.get(i, j) * S.get(i, j);
            }
        }

        double[][] logDelta = new double[c.B][];
        boolean[][] okDelta = new boolean[c.B][];
        for (int b = 1; b < c.B; b++) {
            int t = c.level[b];
            double[][] r = cumlog(c.C[t - 1], c.d[t - 1][b], Ntot);
            logDelta[b] = r[0];
            okDelta[b] = toBool(r[1]);
        }
        double[][] logOmTT = new double[c.T + 1][];
        boolean[][] okOmTT = new boolean[c.T + 1][];
        double[][] logOmPrev = new double[c.T + 1][];
        boolean[][] okOmPrev = new boolean[c.T + 1][];
        for (int t = 1; t <= c.T; t++) {
            double[][] r = cumlog(c.C[t - 1], c.Dtt[t], Ntot);
            logOmTT[t] = r[0];
            okOmTT[t] = toBool(r[1]);
            if (t > 1) {
                double[][] r2 = cumlog(c.C[t - 2], c.Dprev[t], Ntot);
                logOmPrev[t] = r2[0];
                okOmPrev[t] = toBool(r2[1]);
            }
        }

        List<int[]> states = states(M, Nv);
        int ns = states.size();
        double[] logw = new double[ns];
        double lmax = Double.NEGATIVE_INFINITY;
        for (int s = 0; s < ns; s++) {
            int[] flat = states.get(s);
            int[] ni = new int[M];
            for (int i = 0; i < M; i++) {
                int acc = 0;
                for (int j = 0; j < J; j++) {
                    acc += flat[j * M + i];
                }
                ni[i] = acc;
            }
            double lw = 0;
            boolean ok = true;
            for (int i = 0; i < M && ok; i++) {
                lw += lgamma(ni[i] + 1) - logbeta[i][ni[i]];
                for (int j = 0; j < J; j++) {
                    int nij = flat[j * M + i];
                    if (nij > 0) {
                        if (gamma[i][j] <= 0) {
                            ok = false;
                            break;
                        }
                        lw += nij * Math.log(gamma[i][j]) - lgamma(nij + 1);
                    }
                }
            }
            if (!ok) {
                logw[s] = Double.NEGATIVE_INFINITY;
                continue;
            }
            int[] m = new int[c.B];
            for (int b = 1; b < c.B; b++) {
                int acc = 0;
                for (int k = 0; k < c.sdr.branch[b].length; k++) {
                    acc += ni[c.sdr.branch[b][k]];
                }
                m[b] = acc;
                if (!okDelta[b][acc]) {
                    ok = false;
                    break;
                }
                lw += logDelta[b][acc];
            }
            if (!ok) {
                logw[s] = Double.NEGATIVE_INFINITY;
                continue;
            }
            for (int t = 1; t <= c.T && ok; t++) {
                int v = 0;
                for (int k = 0; k < c.inA[t].length; k++) {
                    v += m[c.inA[t][k]];
                }
                if (!okOmTT[t][v]) {
                    ok = false;
                    break;
                }
                lw -= logOmTT[t][v];
                if (t > 1) {
                    if (!okOmPrev[t][v]) {
                        ok = false;
                        break;
                    }
                    lw += logOmPrev[t][v];
                }
            }
            logw[s] = ok ? lw : Double.NEGATIVE_INFINITY;
            if (logw[s] > lmax) {
                lmax = logw[s];
            }
        }
        if (Double.isInfinite(lmax)) {
            throw new RuntimeException("The SDR network has no reachable state at the given populations: "
                    + "the routing coefficients forbid every state.");
        }

        double Gs = 0;
        double[] w = new double[ns];
        for (int s = 0; s < ns; s++) {
            w[s] = Math.exp(logw[s] - lmax);
            Gs += w[s];
        }
        Result res = new Result();
        res.lG = lmax + Math.log(Gs);
        res.Q = new Matrix(M, J);
        res.X = new Matrix(M, J);
        res.U = new Matrix(M, J);
        res.R = new Matrix(M, J);
        for (int s = 0; s < ns; s++) {
            double p = w[s] / Gs;
            if (p == 0) {
                continue;
            }
            int[] flat = states.get(s);
            for (int i = 0; i < M; i++) {
                int nitot = 0;
                for (int j = 0; j < J; j++) {
                    nitot += flat[j * M + i];
                }
                for (int j = 0; j < J; j++) {
                    int nij = flat[j * M + i];
                    if (nij > 0) {
                        res.Q.set(i, j, res.Q.get(i, j) + p * nij);
                    }
                    if (nitot > 0 && S.get(i, j) > 0) {
                        res.X.set(i, j, res.X.get(i, j)
                                + p * alp[i][nitot - 1] * ((double) nij / nitot) / S.get(i, j));
                    }
                }
            }
        }
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < J; j++) {
                res.U.set(i, j, res.X.get(i, j) * S.get(i, j));
                if (res.X.get(i, j) > 0) {
                    res.R.set(i, j, res.Q.get(i, j) / res.X.get(i, j));
                }
            }
        }
        return res;
    }

    /**
     * Coefficients xi of Section 3.2.
     *
     * <p>P is a list with one center-by-center matrix per chain, holding the
     * state-independent routing probabilities. Three rules fix the coefficients:
     * the complement M-V obeys the ordinary traffic equations with the whole SDR
     * subnetwork collapsed into a single e to d arc of probability one; every
     * branch obeys its own traffic equations driven by an injection of xi_e at
     * its entry center; and xi_e is one.</p>
     *
     * <p>The paper states xi_ij = xi_ej for the branch entry and departure
     * centers and works out only single-center branches. The traffic equations
     * above are the reading that extends it: they return xi at the branch
     * departure equal to xi_e because a customer leaves a branch only through
     * it, and xi at the branch entry equal to xi_e whenever that center takes no
     * internal feedback. They have been checked against a brute-force CTMC on a
     * branch that does take such feedback, where the literal rule fails.</p>
     *
     * <p>These xi are not relative visit counts: the rate at which customers
     * enter a branch is state dependent, so a ratio of two xi carries no flow
     * meaning.</p>
     *
     * @param sdr the routing structure, in center indices
     * @param P one center-by-center SIR routing matrix per chain
     * @return the coefficients, centers by chains
     */
    public static Matrix pfqn_sdrvisits(StateDepRouting sdr, List<Matrix> P) {
        Coeff c = pfqn_sdrcoeff(sdr);
        int M = P.get(0).getNumRows();
        int J = P.size();
        Matrix xi = new Matrix(M, J);

        boolean[] inV = new boolean[M];
        for (int b = 1; b < c.B; b++) {
            for (int k = 0; k < sdr.branch[b].length; k++) {
                inV[sdr.branch[b][k]] = true;
            }
        }
        List<Integer> mv = new ArrayList<Integer>();
        for (int i = 0; i < M; i++) {
            if (!inV[i]) {
                mv.add(Integer.valueOf(i));
            }
        }
        int ie = mv.indexOf(Integer.valueOf(sdr.entry));
        int id = mv.indexOf(Integer.valueOf(sdr.departure));
        if (ie < 0 || id < 0) {
            throw new RuntimeException("The entry and departure centers of Q(V,V) must lie outside every branch.");
        }

        for (int j = 0; j < J; j++) {
            Matrix Pj = P.get(j);
            int nm = mv.size();
            double[][] Pmv = new double[nm][nm];
            for (int a = 0; a < nm; a++) {
                for (int b = 0; b < nm; b++) {
                    Pmv[a][b] = Pj.get(mv.get(a).intValue(), mv.get(b).intValue());
                }
            }
            for (int b = 0; b < nm; b++) {
                Pmv[ie][b] = 0;
            }
            Pmv[ie][id] = 1.0;
            for (int a = 0; a < nm; a++) {
                double rs = 0;
                for (int b = 0; b < nm; b++) {
                    rs += Pmv[a][b];
                }
                if (Math.abs(rs - 1.0) > 1e-3) {
                    throw new RuntimeException("The SIR routing of chain " + j
                            + " does not keep customers inside the complement M-V.");
                }
            }
            double[] xmv = dtmcSolve(Pmv);
            if (xmv[ie] <= 0) {
                throw new RuntimeException("The entry center of Q(V,V) is unreachable in chain " + j + ".");
            }
            for (int a = 0; a < nm; a++) {
                xi.set(mv.get(a).intValue(), j, xmv[a] / xmv[ie]);
            }
            for (int b = 1; b < c.B; b++) {
                int[] sb = sdr.branch[b];
                int nb = sb.length;
                double[][] A = new double[nb][nb];
                double[] rhs = new double[nb];
                for (int a = 0; a < nb; a++) {
                    for (int k = 0; k < nb; k++) {
                        // transpose of (I - P_bb), so the solve is the row-vector equation
                        A[k][a] = (a == k ? 1.0 : 0.0) - Pj.get(sb[a], sb[k]);
                    }
                    if (sb[a] == sdr.entryOf[b]) {
                        rhs[a] = xi.get(sdr.entry, j);
                    }
                }
                double[] sol = solve(A, rhs);
                for (int a = 0; a < nb; a++) {
                    xi.set(sb[a], j, sol[a]);
                }
            }
        }
        return xi;
    }

    /**
     * Section 4 mean value analysis and convolution.
     *
     * <p>Same inputs and outputs as {@link #pfqn_sdr}, which evaluates eq. (16)
     * exactly by state enumeration, so the two are directly comparable. This
     * routine costs O(J T M (V_1...V_J)^2) rather than the size of the state
     * space, at the price of two restrictions the paper itself imposes: every
     * SDR branch must hold a single centre, and every C_t must be negative. A
     * C_t other than -1 is rescaled internally, which leaves eqs. (10) and (16)
     * unchanged because the level factors telescope.</p>
     *
     * <p>Two formulas of Section 4 are corrected here, both verified against
     * {@link #pfqn_sdr}. The initialise step of 4.2.2 divides by
     * T_j(V-1_j,V_T) where the convolution identity G(V)/G(V-1_j) = 1/T_j(V)
     * gives T_j(V,V_T); this implementation forms G = g_mva Omega_{T-1,T}/Omega_TT
     * directly instead. And 4.2.3's T_ij = xi_ij [d_1i - Q_i] T_j drops the
     * state-dependent omega ratios of eq. (10); the exact identity is
     * T_ij = xi_ij T_j(N,M) E_{N-1_j}[P_{e,e(i)}].</p>
     *
     * @param S mean service times, centres by chains
     * @param xi coefficients of Section 3.2, centres by chains
     * @param N chain populations
     * @param sdr the routing structure, in centre indices
     * @param alpha load-dependent rate scalings, or null for fixed rate
     * @return the mean performance measures
     */
    public static Result pfqn_sdrmva(Matrix S, Matrix xi, Matrix N, StateDepRouting sdr, Matrix alpha) {
        final int M = S.getNumRows();
        final int J = S.getNumCols();
        Coeff c0 = pfqn_sdrcoeff(sdr);
        for (int b = 1; b < c0.B; b++) {
            if (sdr.branch[b].length != 1) {
                throw new RuntimeException("pfqn_sdrmva requires every SDR branch to hold a single centre: "
                        + "the MVA and convolution of Krzesinski (1987) Section 4 is stated that way and its "
                        + "general case is in an unpublished technical report. Use pfqn_sdr, which evaluates "
                        + "eq. (16) exactly for any branch topology.");
            }
        }
        for (int t = 0; t < c0.T; t++) {
            if (c0.C[t] >= 0) {
                throw new RuntimeException("pfqn_sdrmva requires every C_t to be negative: Section 2.5 assumes "
                        + "it and Section 4 is written for C_t = -1.");
            }
        }
        // Rescale each level to C_t = -1, which leaves eqs. (10) and (16) unchanged
        StateDepRouting sdr1 = sdr.copy();
        for (int t = 0; t < c0.T; t++) {
            double k = -c0.C[t];
            sdr1.C[t] = -1.0;
            for (int b = 0; b < sdr1.d[t].length; b++) sdr1.d[t][b] = sdr.d[t][b] / k;
        }
        Coeff c = pfqn_sdrcoeff(sdr1);

        int[] Nv = new int[J];
        int Ntot = 0;
        for (int j = 0; j < J; j++) {
            Nv[j] = (int) N.get(j);
            Ntot += Nv[j];
        }
        double[][] gamma = new double[M][J];
        for (int i = 0; i < M; i++)
            for (int j = 0; j < J; j++) gamma[i][j] = xi.get(i, j) * S.get(i, j);
        double[][] alp = new double[M][Math.max(1, Ntot)];
        for (int i = 0; i < M; i++)
            for (int k = 0; k < alp[i].length; k++)
                alp[i][k] = (alpha != null && i < alpha.getNumRows() && k < alpha.getNumCols())
                        ? alpha.get(i, k) : 1.0;

        int[][] latt = lattice(Nv);
        int nl = latt.length;

        boolean[] inV = new boolean[M];
        double[] dvec = new double[M];
        int[] lvl = new int[M];
        java.util.Arrays.fill(dvec, Double.NaN);
        for (int b = 1; b < c.B; b++) {
            int i = sdr1.branch[b][0];
            inV[i] = true;
            lvl[i] = c.level[b];
            dvec[i] = c.d[c.level[b] - 1][b];
        }
        List<Integer> mv = new ArrayList<Integer>();
        for (int i = 0; i < M; i++) if (!inV[i]) mv.add(Integer.valueOf(i));

        // Sec. 4.1: the complement is an ordinary BCMP subnetwork, delta = 1
        int[] mvArr = new int[mv.size()];
        double[] mvD = new double[mv.size()];
        for (int k = 0; k < mv.size(); k++) {
            mvArr[k] = mv.get(k).intValue();
            mvD[k] = Double.NaN;
        }
        double[][][] Qc = new double[M][J][nl];
        double[][] Tc = new double[J][nl];
        double[] gc = new double[nl];
        setMva(mvArr, mvD, gamma, alp, Nv, latt, Ntot, M, J, Qc, Tc, gc);

        // Sec. 4.2: outermost level last
        double[] Gin = new double[nl];
        Gin[key(new int[J], Nv)] = 1.0;
        double[][][] Qin = new double[M][J][nl];
        double[][] Tin = new double[J][nl];
        double[][] Ain = new double[M][nl];
        for (int t = c.T; t >= 1; t--) {
            List<Integer> stL = new ArrayList<Integer>();
            List<Integer> inL = new ArrayList<Integer>();
            for (int i = 0; i < M; i++) {
                if (lvl[i] == t) stL.add(Integer.valueOf(i));
                else if (inV[i] && lvl[i] > t) inL.add(Integer.valueOf(i));
            }
            int[] St = new int[stL.size()];
            double[] Sd = new double[stL.size()];
            for (int k = 0; k < St.length; k++) {
                St[k] = stL.get(k).intValue();
                Sd[k] = dvec[St[k]];
            }
            double[][][] Qt = new double[M][J][nl];
            double[][] Tt = new double[J][nl];
            double[] gt = new double[nl];
            setMva(St, Sd, gamma, alp, Nv, latt, Ntot, M, J, Qt, Tt, gt);

            double[] Gnew = new double[nl];
            double[][][] Qnew = new double[M][J][nl];
            double[][] Tnew = new double[J][nl];
            double[][] Anew = new double[M][nl];
            for (int v = 0; v < nl; v++) {
                int[] V = latt[v];
                int vv = 0;
                for (int j = 0; j < J; j++) vv += V[j];
                double omr = omegaCum(c, t, vv);
                if (omr == 0.0) continue;
                double[] anum = new double[M];
                for (int[] L : sublattice(V)) {
                    int[] VL = new int[J];
                    for (int j = 0; j < J; j++) VL[j] = V[j] - L[j];
                    int iL = key(L, Nv), iVL = key(VL, Nv);
                    double pb = omr * gt[iVL] * Gin[iL];
                    if (pb == 0.0) continue;
                    Gnew[v] += pb;
                    for (int k = 0; k < St.length; k++)
                        for (int j = 0; j < J; j++) Qnew[St[k]][j][v] += Qt[St[k]][j][iVL] * pb;
                    for (int j = 0; j < J; j++) Tnew[j][v] += Tt[j][iVL] * pb;
                    for (int q = 0; q < inL.size(); q++) {
                        int i = inL.get(q).intValue();
                        for (int j = 0; j < J; j++) Qnew[i][j][v] += Qin[i][j][iL] * pb;
                        anum[i] += Ain[i][iL] * pb;
                    }
                }
                if (Gnew[v] > 0) {
                    for (int i = 0; i < M; i++)
                        for (int j = 0; j < J; j++) Qnew[i][j][v] /= Gnew[v];
                    for (int j = 0; j < J; j++) Tnew[j][v] /= Gnew[v];
                    // eq. (10) carries, besides delta_ti, the single-step ratios
                    // omega_{s-1,s}(v_s)/omega_ss(v_s) down to the centre's own level.
                    // Conditioning on the population of Q(V,V_t) fixes v_t, so that
                    // ratio leaves the expectation and the rest recurses through the
                    // same convolution as the queue lengths.
                    double st = omegaStep(c, t, vv);
                    for (int k = 0; k < St.length; k++) {
                        double qi = 0;
                        for (int j = 0; j < J; j++) qi += Qnew[St[k]][j][v];
                        Anew[St[k]][v] = st * (dvec[St[k]] - qi);
                    }
                    for (int q = 0; q < inL.size(); q++) {
                        int i = inL.get(q).intValue();
                        Anew[i][v] = st * anum[i] / Gnew[v];
                    }
                }
            }
            Gin = Gnew; Qin = Qnew; Tin = Tnew; Ain = Anew;
        }

        // Sec. 4.3: convolve against the complement at every population, because
        // the per-centre throughputs below read the network at N - 1_j
        double[][][] Qall = new double[M][J][nl];
        double[][] Tall = new double[J][nl];
        double[][] Aall = new double[M][nl];
        double[] Gall = new double[nl];
        for (int v = 0; v < nl; v++) {
            int[] Np = latt[v];
            for (int[] V : sublattice(Np)) {
                int[] C2 = new int[J];
                for (int j = 0; j < J; j++) C2[j] = Np[j] - V[j];
                int iV = key(V, Nv), iC = key(C2, Nv);
                double pb = gc[iC] * Gin[iV];
                if (pb == 0.0) continue;
                Gall[v] += pb;
                for (int k = 0; k < mvArr.length; k++)
                    for (int j = 0; j < J; j++) Qall[mvArr[k]][j][v] += Qc[mvArr[k]][j][iC] * pb;
                for (int i = 0; i < M; i++) {
                    if (!inV[i]) continue;
                    for (int j = 0; j < J; j++) Qall[i][j][v] += Qin[i][j][iV] * pb;
                    Aall[i][v] += Ain[i][iV] * pb;
                }
                for (int j = 0; j < J; j++) Tall[j][v] += Tc[j][iC] * pb;
            }
            if (Gall[v] > 0) {
                for (int i = 0; i < M; i++) {
                    for (int j = 0; j < J; j++) Qall[i][j][v] /= Gall[v];
                    Aall[i][v] /= Gall[v];
                }
                for (int j = 0; j < J; j++) Tall[j][v] /= Gall[v];
            }
        }

        int vN = key(Nv, Nv);
        Result res = new Result();
        res.Q = new Matrix(M, J);
        res.X = new Matrix(M, J);
        res.U = new Matrix(M, J);
        res.R = new Matrix(M, J);
        res.lG = Math.log(Gall[vN]);
        for (int i = 0; i < M; i++)
            for (int j = 0; j < J; j++) res.Q.set(i, j, Qall[i][j][vN]);
        for (int j = 0; j < J; j++) {
            if (Nv[j] == 0) continue;
            int[] Nm = Nv.clone();
            Nm[j]--;
            int vm = key(Nm, Nv);
            for (int i = 0; i < M; i++) {
                double x = inV[i] ? xi.get(i, j) * Tall[j][vN] * Aall[i][vm]
                                  : xi.get(i, j) * Tall[j][vN];
                res.X.set(i, j, x);
            }
        }
        for (int i = 0; i < M; i++)
            for (int j = 0; j < J; j++) {
                res.U.set(i, j, res.X.get(i, j) * S.get(i, j));
                if (res.X.get(i, j) > 0) res.R.set(i, j, res.Q.get(i, j) / res.X.get(i, j));
            }
        return res;
    }

    /**
     * Section 4.2.1: MVA of the centres cidx in isolation. dvec[k] is the SDR
     * admission coefficient of centre cidx[k], so delta_k(n) = dvec[k] - n; NaN
     * means delta = 1, the ordinary BCMP arrival theorem, which the complement
     * M-V uses.
     *
     * <p>The factor [d_ti - Q_i(V - 1_j)] that the paper places in both W_ij and
     * the denominator of T_j cancels between them and is not formed here: it can
     * vanish at the population bound, where it would divide by zero without
     * changing any queue length or throughput.</p>
     */
    private static void setMva(int[] cidx, double[] dvec, double[][] gamma, double[][] alp,
                               int[] N, int[][] latt, int Ntot, int M, int J,
                               double[][][] Qs, double[][] Ts, double[] gs) {
        final int nc = cidx.length;
        final int nl = latt.length;
        int z = key(new int[J], N);
        if (nc == 0) {
            for (int v = 0; v < nl; v++) {
                int tot = 0;
                for (int j = 0; j < J; j++) tot += latt[v][j];
                gs[v] = (tot == 0) ? 1.0 : 0.0;
            }
            return;
        }
        double[][][] Ps = new double[nc][Ntot + 1][nl];
        gs[z] = 1.0;
        for (int k = 0; k < nc; k++) Ps[k][0][z] = 1.0;

        Integer[] ord = new Integer[nl];
        for (int v = 0; v < nl; v++) ord[v] = Integer.valueOf(v);
        final int[][] lattF = latt;
        final int JF = J;
        java.util.Arrays.sort(ord, new java.util.Comparator<Integer>() {
            public int compare(Integer a, Integer b) {
                int sa = 0, sb = 0;
                for (int j = 0; j < JF; j++) { sa += lattF[a.intValue()][j]; sb += lattF[b.intValue()][j]; }
                return sa - sb;
            }
        });

        for (int oi = 0; oi < nl; oi++) {
            int v = ord[oi].intValue();
            int[] V = latt[v];
            int vv = 0;
            for (int j = 0; j < J; j++) vv += V[j];
            if (vv == 0) continue;
            double[][] A = new double[nc][J];
            int[] vm = new int[J];
            java.util.Arrays.fill(vm, -1);
            for (int j = 0; j < J; j++) {
                if (V[j] == 0) continue;
                int[] Vm = V.clone();
                Vm[j]--;
                vm[j] = key(Vm, N);
                for (int k = 0; k < nc; k++) {
                    double acc = 0;
                    for (int n = 1; n <= vv; n++) {
                        double df = Double.isNaN(dvec[k]) ? 1.0 : dvec[k] - (n - 1);
                        if (df <= 0) break;
                        acc += n * (df / alp[cidx[k]][n - 1]) * Ps[k][n - 1][vm[j]];
                    }
                    A[k][j] = acc;
                }
            }
            for (int j = 0; j < J; j++) {
                if (V[j] == 0) continue;
                double den = 0;
                for (int k = 0; k < nc; k++) den += gamma[cidx[k]][j] * A[k][j];
                if (den > 0) Ts[j][v] = V[j] / den;
            }
            for (int k = 0; k < nc; k++)
                for (int j = 0; j < J; j++) {
                    if (V[j] == 0) continue;
                    Qs[cidx[k]][j][v] = gamma[cidx[k]][j] * Ts[j][v] * A[k][j];
                }
            for (int k = 0; k < nc; k++) {
                double tot = 0;
                for (int n = 1; n <= vv; n++) {
                    double df = Double.isNaN(dvec[k]) ? 1.0 : dvec[k] - (n - 1);
                    if (df <= 0) break;
                    double acc = 0;
                    for (int j = 0; j < J; j++) {
                        if (V[j] == 0) continue;
                        acc += gamma[cidx[k]][j] * Ts[j][v] * Ps[k][n - 1][vm[j]];
                    }
                    Ps[k][n][v] = (df / alp[cidx[k]][n - 1]) * acc;
                    tot += Ps[k][n][v];
                }
                Ps[k][0][v] = 1.0 - tot;
            }
            for (int j = 0; j < J; j++)
                if (V[j] > 0 && Ts[j][v] > 0) {
                    gs[v] = gs[vm[j]] / Ts[j][v];
                    break;
                }
        }
    }

    /** Omega_{t-1,t}(v)/Omega_tt(v), the cumulative ratio of eq. (16). */
    private static double omegaCum(Coeff c, int t, int v) {
        double num = 1, den = 1;
        for (int k = 0; k < v; k++) {
            double f = -k + c.Dtt[t];
            if (f <= 0) return 0.0;
            den *= f;
            if (t > 1) {
                double g = -k + c.Dprev[t];
                if (g <= 0) return 0.0;
                num *= g;
            }
        }
        return num / den;
    }

    /**
     * omega_{t-1,t}(v)/omega_tt(v), the single-step ratio of eq. (10). Distinct
     * from {@link #omegaCum}: the routing probability carries the lowercase
     * omega, the normalizing constant the uppercase one.
     */
    private static double omegaStep(Coeff c, int t, int v) {
        double den = -v + c.Dtt[t];
        if (den <= 0) return 0.0;
        double num = 1.0;  // omega_{0,1} is one
        if (t > 1) {
            num = -v + c.Dprev[t];
            if (num <= 0) return 0.0;
        }
        return num / den;
    }

    /** Every population vector V with 0 <= V <= N, mixed-radix ordered. */
    private static int[][] lattice(int[] N) {
        int tot = 1;
        for (int j = 0; j < N.length; j++) tot *= N[j] + 1;
        int[][] out = new int[tot][N.length];
        for (int r = 0; r < tot; r++) {
            int rem = r;
            for (int j = 0; j < N.length; j++) {
                out[r][j] = rem % (N[j] + 1);
                rem /= (N[j] + 1);
            }
        }
        return out;
    }

    /** Mixed-radix row of V in the lattice of {@link #lattice}. */
    private static int key(int[] V, int[] N) {
        int k = 0, mul = 1;
        for (int j = 0; j < N.length; j++) {
            k += V[j] * mul;
            mul *= N[j] + 1;
        }
        return k;
    }

    /** Every L with 0 <= L <= V. */
    private static int[][] sublattice(int[] V) {
        int tot = 1;
        for (int j = 0; j < V.length; j++) tot *= V[j] + 1;
        int[][] out = new int[tot][V.length];
        for (int r = 0; r < tot; r++) {
            int rem = r;
            for (int j = 0; j < V.length; j++) {
                out[r][j] = rem % (V[j] + 1);
                rem /= (V[j] + 1);
            }
        }
        return out;
    }

    /**
     * Cumulative log-product of Cv*k + dv for k = 0..n-1, n = 0..nmax.
     * Row 0 holds the logs, row 1 a nonzero flag that turns off once a
     * nonpositive factor has been met, which is where the SDR population bound
     * closes the branch or the subnetwork.
     */
    private static double[][] cumlog(double Cv, double dv, int nmax) {
        double[] lc = new double[nmax + 1];
        double[] ok = new double[nmax + 1];
        ok[0] = 1.0;
        for (int n = 1; n <= nmax; n++) {
            double f = Cv * (n - 1) + dv;
            if (ok[n - 1] == 0.0 || f <= 0) {
                ok[n] = 0.0;
                lc[n] = Double.NEGATIVE_INFINITY;
            } else {
                ok[n] = 1.0;
                lc[n] = lc[n - 1] + Math.log(f);
            }
        }
        return new double[][]{lc, ok};
    }

    private static boolean[] toBool(double[] v) {
        boolean[] out = new boolean[v.length];
        for (int i = 0; i < v.length; i++) {
            out[i] = v[i] != 0.0;
        }
        return out;
    }

    /** Every population matrix with the given column sums, flattened chain-major. */
    private static List<int[]> states(int M, int[] N) {
        List<int[]> acc = new ArrayList<int[]>();
        acc.add(new int[0]);
        for (int j = 0; j < N.length; j++) {
            List<int[]> comps = compositions(N[j], M);
            List<int[]> next = new ArrayList<int[]>(acc.size() * comps.size());
            for (int a = 0; a < acc.size(); a++) {
                for (int b = 0; b < comps.size(); b++) {
                    int[] head = acc.get(a);
                    int[] tail = comps.get(b);
                    int[] joined = new int[head.length + tail.length];
                    System.arraycopy(head, 0, joined, 0, head.length);
                    System.arraycopy(tail, 0, joined, head.length, tail.length);
                    next.add(joined);
                }
            }
            acc = next;
        }
        return acc;
    }

    /** All nonnegative integer m-vectors summing to n. */
    private static List<int[]> compositions(int n, int m) {
        List<int[]> out = new ArrayList<int[]>();
        if (m == 1) {
            out.add(new int[]{n});
            return out;
        }
        for (int k = 0; k <= n; k++) {
            List<int[]> tail = compositions(n - k, m - 1);
            for (int t = 0; t < tail.size(); t++) {
                int[] row = new int[m];
                row[0] = k;
                System.arraycopy(tail.get(t), 0, row, 1, m - 1);
                out.add(row);
            }
        }
        return out;
    }

    /** Stationary distribution of a stochastic matrix, by solving the balance equations. */
    private static double[] dtmcSolve(double[][] Pm) {
        int n = Pm.length;
        double[][] A = new double[n][n];
        double[] rhs = new double[n];
        for (int a = 0; a < n - 1; a++) {
            for (int b = 0; b < n; b++) {
                A[a][b] = Pm[b][a] - (a == b ? 1.0 : 0.0);
            }
        }
        for (int b = 0; b < n; b++) {
            A[n - 1][b] = 1.0;
        }
        rhs[n - 1] = 1.0;
        return solve(A, rhs);
    }

    /** Dense linear solve by Gaussian elimination with partial pivoting. */
    private static double[] solve(double[][] Ain, double[] bin) {
        int n = bin.length;
        double[][] A = new double[n][];
        double[] b = bin.clone();
        for (int i = 0; i < n; i++) {
            A[i] = Ain[i].clone();
        }
        for (int col = 0; col < n; col++) {
            int piv = col;
            for (int r = col + 1; r < n; r++) {
                if (Math.abs(A[r][col]) > Math.abs(A[piv][col])) {
                    piv = r;
                }
            }
            if (Math.abs(A[piv][col]) < 1e-14) {
                throw new RuntimeException("Singular system while solving the SDR traffic equations.");
            }
            double[] tmp = A[col];
            A[col] = A[piv];
            A[piv] = tmp;
            double tb = b[col];
            b[col] = b[piv];
            b[piv] = tb;
            for (int r = col + 1; r < n; r++) {
                double f = A[r][col] / A[col][col];
                if (f == 0) {
                    continue;
                }
                for (int k = col; k < n; k++) {
                    A[r][k] -= f * A[col][k];
                }
                b[r] -= f * b[col];
            }
        }
        double[] x = new double[n];
        for (int r = n - 1; r >= 0; r--) {
            double acc = b[r];
            for (int k = r + 1; k < n; k++) {
                acc -= A[r][k] * x[k];
            }
            x[r] = acc / A[r][r];
        }
        return x;
    }

    /** Log-gamma, used for the factorials of f_i(n_i). */
    private static double lgamma(double x) {
        return org.apache.commons.math3.special.Gamma.logGamma(x);
    }
}
