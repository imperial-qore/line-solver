package jline.api.pfqn;

import jline.api.pfqn.ld.Pfqn_mvald;
import jline.api.pfqn.mva.Pfqn_mva;
import jline.io.Ret;
import jline.lang.processes.Coxian;
import jline.util.matrix.Matrix;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Marie's iterative aggregation-decomposition (Marie 1979/1980) for closed
 * queueing networks with FCFS non-exponential (Coxian) service.
 *
 * <p>Single class only: each station is analyzed in isolation as a
 * lambda(n)/Cox/1 queue whose state-dependent arrival rate is the complementary
 * network throughput; the resulting conditional throughputs mu_i(n) drive an
 * exact load-dependent product-form solve (Pfqn_mvald), iterated to a fixed
 * point. The method reduces to exact product form for exponential service.
 *
 * <p>Multiclass (R &gt; 1) is handled by {@link #pfqn_marie_multi}: a QD-AMVA
 * aggregate with class-dependent (cd) scaling beta_{i,r}(nvec) supplied by a
 * multiclass lambda_r/Cox/1 FCFS isolation sub-model, iterated to a fixed
 * point on the per-class throughput. beta==1 recovers standard FCFS AMVA, so
 * the cd-scaling carries only the non-exponential correction; exponential
 * class-independent service dispatches to exact MVA.
 *
 * Ported from matlab/src/api/pfqn/pfqn_marie.m (both single- and multiclass).
 */
public final class Pfqn_marie {

    private Pfqn_marie() {}

    /** Result of {@link #pfqn_marie}. */
    public static class Result {
        public double X;      // chain throughput (reference station)
        public Matrix Q;      // mean queue length per station (M x 1)
        public Matrix U;      // utilization per station (M x 1)
        public Matrix C;      // residence time per station (M x 1)
        public int iter;      // iterations performed
        public Matrix mu;     // converged LD rate-multiplier matrix (M x N)

        public Result(double X, Matrix Q, Matrix U, Matrix C, int iter, Matrix mu) {
            this.X = X; this.Q = Q; this.U = U; this.C = C; this.iter = iter; this.mu = mu;
        }
    }

    public static Result pfqn_marie(Matrix L, double N, double Z, Matrix scv) {
        return pfqn_marie(L, N, Z, scv, 1e-8, 1000, null);
    }

    /**
     * @param L   service demand column vector (M x 1)
     * @param N   total population (scalar)
     * @param Z   think time (total delay demand)
     * @param scv per-station squared coefficient of variation (M x 1); null =&gt; all 1
     * @param tol convergence tolerance on mu
     * @param maxiter maximum iterations
     * @param nservers per-station server count (M x 1); null =&gt; all 1
     */
    public static Result pfqn_marie(Matrix L, double N, double Z, Matrix scv,
                                    double tol, int maxiter, Matrix nservers) {
        int M = L.length();
        int Nint = (int) Math.round(N);

        double[] Lv = new double[M];
        for (int i = 0; i < M; i++) Lv[i] = L.get(i);

        double[] scvv = new double[M];
        for (int i = 0; i < M; i++) scvv[i] = (scv == null) ? 1.0 : scv.get(i);

        double[] mvec = new double[M];
        for (int i = 0; i < M; i++) mvec[i] = (nservers == null) ? 1.0 : nservers.get(i);

        // Per-station Coxian phase representation (mean = L(i), scv(i)).
        double[][] phRate = new double[M][];
        double[][] phCompl = new double[M][];
        for (int i = 0; i < M; i++) {
            Coxian cx = Coxian.fitMeanAndSCV(Lv[i], scvv[i]);
            Matrix mu_i = cx.getMu();
            Matrix phi_i = cx.getPhi();
            int P = mu_i.length();
            phRate[i] = new double[P];
            phCompl[i] = new double[P];
            for (int k = 0; k < P; k++) { phRate[i][k] = mu_i.get(k); phCompl[i][k] = phi_i.get(k); }
        }

        // Initial LD rate multipliers: min(n, m_i).
        double[][] mu = new double[M][Nint];
        for (int i = 0; i < M; i++)
            for (int n = 1; n <= Nint; n++) mu[i][n - 1] = Math.min(n, mvec[i]);

        Matrix Npop = new Matrix(1, 1); Npop.set(0, 0, N);
        Matrix Zmat = new Matrix(1, 1); Zmat.set(0, 0, Z);

        double Xchain = 0.0;
        Matrix QN = new Matrix(M, 1), UN = new Matrix(M, 1), CN = new Matrix(M, 1);
        int it = 0;
        while (it < maxiter) {
            it++;
            Matrix muMat = new Matrix(M, Nint);
            for (int i = 0; i < M; i++) for (int n = 0; n < Nint; n++) muMat.set(i, n, mu[i][n]);

            Ret.pfqnMVALD r = Pfqn_mvald.pfqn_mvald(L, Npop, Zmat, muMat);
            Matrix Pg = r.pi; // M x (N+1): Pg(i,k) = P(n_i = k)

            double[][] muNew = new double[M][Nint];
            double delta = 0.0;
            for (int i = 0; i < M; i++) {
                double[] lam = new double[Nint]; // lam[n] holds lambda_i(n), n=0..N-1
                for (int n = 0; n < Nint; n++) {
                    double pn = Pg.get(i, n);
                    lam[n] = (pn > 0) ? (mu[i][n] / Lv[i]) * Pg.get(i, n + 1) / pn : 0.0;
                }
                double[] muabs = isolCondTput(lam, phRate[i], phCompl[i], Nint, mvec[i]);
                for (int n = 0; n < Nint; n++) {
                    muNew[i][n] = muabs[n] * Lv[i];
                    delta = Math.max(delta, Math.abs(muNew[i][n] - mu[i][n]));
                }
            }
            for (int i = 0; i < M; i++) System.arraycopy(muNew[i], 0, mu[i], 0, Nint);

            Xchain = r.X.get(0);
            QN = r.Q.copy(); UN = r.U.copy(); CN = r.R.copy();
            if (delta < tol) break;
        }

        Matrix muMat = new Matrix(M, Nint);
        for (int i = 0; i < M; i++) for (int n = 0; n < Nint; n++) muMat.set(i, n, mu[i][n]);
        return new Result(Xchain, QN, UN, CN, it, muMat);
    }

    /**
     * Stationary analysis of a lambda(n)/Cox/1(-m) queue in isolation; returns
     * the conditional throughput mu(n) = departure rate given n present, for
     * n = 1..N (indexed 0..N-1).
     */
    private static double[] isolCondTput(double[] lam, double[] phRate, double[] phCompl, int N, double m) {
        int P = phRate.length;
        int S = 1 + N * P; // state 0 = empty; (n,k) -> 1 + (n-1)*P + (k)
        double[][] Gq = new double[S][S];

        // From empty: arrival starts a customer in phase 1.
        Gq[0][idx(1, 0, P)] += lam[0];

        for (int n = 1; n <= N; n++) {
            double sc = Math.min(n, m);
            for (int k = 0; k < P; k++) {
                int rrow = idx(n, k, P);
                if (n < N) Gq[rrow][idx(n + 1, k, P)] += lam[n];
                double compl = phRate[k] * phCompl[k] * sc;
                double adv = phRate[k] * (1.0 - phCompl[k]) * sc;
                if (adv > 0 && k < P - 1) Gq[rrow][idx(n, k + 1, P)] += adv;
                if (compl > 0) {
                    if (n > 1) Gq[rrow][idx(n - 1, 0, P)] += compl;
                    else Gq[rrow][0] += compl;
                }
            }
        }
        // Diagonal = -rowsum.
        for (int i = 0; i < S; i++) {
            double rowsum = 0.0;
            for (int j = 0; j < S; j++) rowsum += Gq[i][j];
            Gq[i][i] -= rowsum;
        }

        // Solve p * Gq = 0, sum(p) = 1: square system Gq' with last row replaced
        // by the normalization equation.
        Matrix A = new Matrix(S, S);
        for (int i = 0; i < S; i++)
            for (int j = 0; j < S; j++) A.set(i, j, Gq[j][i]); // transpose
        for (int j = 0; j < S; j++) A.set(S - 1, j, 1.0);
        Matrix b = new Matrix(S, 1);
        b.set(S - 1, 0, 1.0);
        Matrix pv = A.leftMatrixDivide(b);

        double[] muvec = new double[N];
        double invSum = 0.0;
        for (int k = 0; k < P; k++) invSum += 1.0 / phRate[k];
        for (int n = 1; n <= N; n++) {
            double Pn = 0.0, dep = 0.0;
            for (int k = 0; k < P; k++) {
                double pk = pv.get(idx(n, k, P));
                Pn += pk;
                dep += pk * phRate[k] * phCompl[k] * Math.min(n, m);
            }
            muvec[n - 1] = (Pn > 0) ? dep / Pn : Math.min(n, m) / invSum;
        }
        return muvec;
    }

    private static int idx(int n, int k, int P) {
        return 1 + (n - 1) * P + k;
    }

    // ================= Multiclass Marie (QD-AMVA + cd-scaling) =================

    /** Result of {@link #pfqn_marie_multi}. */
    public static class MultiResult {
        public Matrix X;   // 1 x R chain throughput
        public Matrix Q;   // M x R mean queue length
        public Matrix U;   // M x R utilization (busy fraction)
        public Matrix C;   // M x R residence (work-based W)
        public int iter;   // outer iterations performed

        public MultiResult(Matrix X, Matrix Q, Matrix U, Matrix C, int iter) {
            this.X = X; this.Q = Q; this.U = U; this.C = C; this.iter = iter;
        }
    }

    public static MultiResult pfqn_marie_multi(Matrix L, Matrix N, Matrix Z, Matrix scv) {
        return pfqn_marie_multi(L, N, Z, scv, 1e-8, 1000);
    }

    /**
     * Marie's method for multiclass FCFS Coxian closed networks. The aggregate
     * is QD-AMVA with class-dependent scaling beta_{i,r}(nvec) = (Coxian
     * isolation conditional throughput)/(exponential isolation conditional
     * throughput), iterated to a fixed point on the per-class throughput X.
     *
     * @param L   M x R service demand matrix
     * @param N   1 x R (or R x 1) per-class population
     * @param Z   1 x R (or R x 1) per-class think time; null =&gt; zeros
     * @param scv M x R per-station per-class SCV; null =&gt; all 1
     */
    public static MultiResult pfqn_marie_multi(Matrix L, Matrix N, Matrix Z, Matrix scv,
                                               double tol, int maxiter) {
        int M = L.getNumRows();
        int R = L.getNumCols();

        int[] Nvec = new int[R];
        for (int r = 0; r < R; r++) Nvec[r] = (int) Math.round(N.get(r));
        double[] Zv = new double[R];
        for (int r = 0; r < R; r++) Zv[r] = (Z == null) ? 0.0 : Z.get(r);
        double[][] Lv = new double[M][R];
        double[][] scvv = new double[M][R];
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                Lv[i][r] = L.get(i, r);
                scvv[i][r] = (scv == null) ? 1.0 : scv.get(i, r);
            }
        }

        // Exact product-form dispatch: exponential + class-independent -> exact MVA.
        boolean isPF = true;
        for (int i = 0; i < M && isPF; i++)
            for (int r = 0; r < R; r++) if (Math.abs(scvv[i][r] - 1.0) > 1e-12) { isPF = false; break; }
        if (isPF) {
            for (int i = 0; i < M; i++) {
                double lo = Double.POSITIVE_INFINITY, hi = Double.NEGATIVE_INFINITY;
                for (int r = 0; r < R; r++) { lo = Math.min(lo, Lv[i][r]); hi = Math.max(hi, Lv[i][r]); }
                if (hi - lo > 1e-12) { isPF = false; break; }
            }
        }
        if (isPF) {
            Matrix Nrow = new Matrix(1, R); for (int r = 0; r < R; r++) Nrow.set(0, r, Nvec[r]);
            Matrix Zrow = new Matrix(1, R); for (int r = 0; r < R; r++) Zrow.set(0, r, Zv[r]);
            Ret.pfqnMVA mr = Pfqn_mva.pfqn_mva(L, Nrow, Zrow);
            Matrix Xr = new Matrix(1, R);
            for (int r = 0; r < R; r++) Xr.set(0, r, mr.X.get(r));
            return new MultiResult(Xr, mr.Q.copy(), mr.U.copy(), mr.R.copy(), 0);
        }

        // Coxian phase representation per (i,r), plus an exponential reference.
        double[][][] phR = new double[M][R][];
        double[][][] phP = new double[M][R][];
        double[][][] eR  = new double[M][R][];
        double[][][] eP  = new double[M][R][];
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                Coxian cx = Coxian.fitMeanAndSCV(Lv[i][r], scvv[i][r]);
                Matrix mm = cx.getMu(); Matrix pp = cx.getPhi();
                int P = mm.length();
                phR[i][r] = new double[P]; phP[i][r] = new double[P];
                for (int k = 0; k < P; k++) { phR[i][r][k] = mm.get(k); phP[i][r][k] = pp.get(k); }
                eR[i][r] = new double[]{1.0 / Lv[i][r]};
                eP[i][r] = new double[]{1.0};
            }
        }

        int npops = 1; for (int r = 0; r < R; r++) npops *= (Nvec[r] + 1);
        // cd-scaling lattices per station: muCox[i][r][p], muExp[i][r][p]; null => beta=1.
        double[][][] muCox = new double[M][][];
        double[][][] muExp = new double[M][][];

        double[] X = new double[R];
        double[] Xprev = new double[R]; for (int r = 0; r < R; r++) Xprev[r] = Double.POSITIVE_INFINITY;
        double[][] Q = new double[M][R], U = new double[M][R], W = new double[M][R];
        int it = 0;
        while (it < maxiter) {
            it++;
            amvaQd(Lv, Nvec, Zv, muCox, muExp, npops, X, Q, U, W);
            for (int i = 0; i < M; i++) {
                muCox[i] = isolMc(X, phR[i], phP[i], Nvec, npops);
                muExp[i] = isolMc(X, eR[i],  eP[i],  Nvec, npops);
            }
            double d = 0.0;
            for (int r = 0; r < R; r++) d = Math.max(d, Math.abs(X[r] - Xprev[r]));
            if (d < tol) break;
            System.arraycopy(X, 0, Xprev, 0, R);
        }

        Matrix Xr = new Matrix(1, R); for (int r = 0; r < R; r++) Xr.set(0, r, X[r]);
        Matrix Qm = new Matrix(M, R), Um = new Matrix(M, R), Cm = new Matrix(M, R);
        for (int i = 0; i < M; i++)
            for (int r = 0; r < R; r++) { Qm.set(i, r, Q[i][r]); Um.set(i, r, U[i][r]); Cm.set(i, r, W[i][r]); }
        return new MultiResult(Xr, Qm, Um, Cm, it);
    }

    /**
     * Multiclass Schweitzer AMVA with class-dependent service-rate scaling. The
     * effective class-r demand at station i is L(i,r)/beta_{i,r}(nvec_arrival),
     * beta supplied by the cd-scaling lattices (null =&gt; beta=1). Writes X, Q,
     * U, W in place.
     */
    private static void amvaQd(double[][] L, int[] Nvec, double[] Z,
                               double[][][] muCox, double[][][] muExp, int npops,
                               double[] X, double[][] Q, double[][] U, double[][] W) {
        int M = L.length, R = Nvec.length;
        double denomM = Math.max(M, 1);
        for (int i = 0; i < M; i++) for (int r = 0; r < R; r++) Q[i][r] = Nvec[r] / denomM;
        double[][] Qprev = new double[M][R];
        for (int i = 0; i < M; i++) for (int r = 0; r < R; r++) Qprev[i][r] = Q[i][r] + 1.0;
        double tol = 1e-9; int it = 0;
        double[] nv = new double[R];
        while (it < 5000) {
            double diff = 0.0;
            for (int i = 0; i < M; i++) for (int r = 0; r < R; r++) diff = Math.max(diff, Math.abs(Q[i][r] - Qprev[i][r]));
            if (diff <= tol) break;
            it++;
            for (int i = 0; i < M; i++) for (int r = 0; r < R; r++) Qprev[i][r] = Q[i][r];
            for (int r = 0; r < R; r++) {
                for (int i = 0; i < M; i++) {
                    for (int s = 0; s < R; s++) nv[s] = Q[i][s];
                    if (Nvec[r] > 0) nv[r] = Q[i][r] * (Nvec[r] - 1.0) / Nvec[r];
                    double[] be = betaEval(muCox[i], muExp[i], nv, Nvec, npops);
                    double wir = L[i][r] / be[r];
                    for (int s = 0; s < R; s++) wir += (L[i][s] / be[s]) * nv[s];
                    W[i][r] = wir;
                }
                double denom = Z[r];
                for (int i = 0; i < M; i++) denom += W[i][r];
                X[r] = (denom > 0) ? Nvec[r] / denom : 0.0;
                for (int i = 0; i < M; i++) Q[i][r] = X[r] * W[i][r];
            }
        }
        for (int i = 0; i < M; i++)
            for (int r = 0; r < R; r++) U[i][r] = X[r] * L[i][r];
    }

    /** beta_{r}(nv) = muCox_r(nv)/muExp_r(nv), clamped; beta=1 when lattices null. */
    private static double[] betaEval(double[][] muCox, double[][] muExp, double[] nv, int[] Nvec, int npops) {
        int R = Nvec.length;
        double[] be = new double[R];
        for (int r = 0; r < R; r++) be[r] = 1.0;
        if (muCox == null || muExp == null) return be;
        for (int r = 0; r < R; r++) {
            double num = ndLinInterp(muCox[r], nv, Nvec);
            double den = ndLinInterp(muExp[r], nv, Nvec);
            if (den > 0 && num > 0 && isFinite(num) && isFinite(den)) be[r] = num / den;
            be[r] = Math.min(Math.max(be[r], 1e-3), 1e3);
        }
        return be;
    }

    private static boolean isFinite(double v) { return !Double.isNaN(v) && !Double.isInfinite(v); }

    // Row-major linear index over the population box [0..Nvec], r=0 most significant.
    private static int linIndex(int[] nvec, int[] Nvec) {
        int idx = 0;
        for (int r = 0; r < Nvec.length; r++) idx = idx * (Nvec[r] + 1) + nvec[r];
        return idx;
    }

    private static int[] decode(int p, int[] Nvec) {
        int R = Nvec.length;
        int[] nvec = new int[R];
        for (int r = R - 1; r >= 0; r--) { int b = Nvec[r] + 1; nvec[r] = p % b; p /= b; }
        return nvec;
    }

    /** Multilinear interpolation of the flat lattice A at real point x, clamped to [0,Nvec]. */
    private static double ndLinInterp(double[] A, double[] x, int[] Nvec) {
        int R = Nvec.length;
        double[] xc = new double[R];
        int[] lo = new int[R], hi = new int[R];
        double[] fr = new double[R];
        for (int r = 0; r < R; r++) {
            xc[r] = Math.min(Math.max(x[r], 0.0), Nvec[r]);
            lo[r] = (int) Math.floor(xc[r]);
            hi[r] = Math.min(lo[r] + 1, Nvec[r]);
            fr[r] = xc[r] - lo[r];
        }
        double v = 0.0;
        int[] sub = new int[R];
        for (int mask = 0; mask < (1 << R); mask++) {
            double w = 1.0;
            for (int r = 0; r < R; r++) {
                if (((mask >> r) & 1) == 1) { sub[r] = hi[r]; w *= fr[r]; }
                else { sub[r] = lo[r]; w *= (1.0 - fr[r]); }
            }
            if (w == 0.0) continue;
            v += w * A[linIndex(sub, Nvec)];
        }
        return v;
    }

    /**
     * Stationary analysis of a multiclass lambda_r/Cox/1 FCFS queue in isolation
     * over the joint per-class population box [0..Nvec]. The head-of-line job is
     * tracked as (class, phase); on a departure the next head class is drawn in
     * random order (prob n_c/sum n). Returns per-class conditional throughput
     * mu_r(nvec) = (class-r departure rate in states with population nvec) /
     * P(nvec) on the flat population lattice.
     */
    private static double[][] isolMc(double[] lam, double[][] phRrow, double[][] phProw,
                                     int[] Nvec, int npops) {
        int R = lam.length;
        int[] Pc = new int[R];
        for (int r = 0; r < R; r++) Pc[r] = phRrow[r].length;

        // Enumerate states: id 0 = empty; (p,c,k) for occupied lattice points.
        List<int[]> ids = new ArrayList<int[]>();
        ids.add(new int[]{-1, -1, -1}); // empty placeholder
        Map<Long, Integer> key2id = new HashMap<Long, Integer>();
        for (int p = 0; p < npops; p++) {
            int[] nvec = decode(p, Nvec);
            int sum = 0; for (int r = 0; r < R; r++) sum += nvec[r];
            if (sum == 0) continue;
            for (int c = 0; c < R; c++) {
                if (nvec[c] > 0) {
                    for (int k = 0; k < Pc[c]; k++) {
                        int id = ids.size();
                        ids.add(new int[]{p, c, k});
                        key2id.put(stateKey(p, c, k, R), id);
                    }
                }
            }
        }
        int S = ids.size();

        double[][] Gq = new double[S][S];
        for (int s = 0; s < S; s++) {
            int[] info = ids.get(s);
            int p = info[0], c = info[1], k = info[2];
            int[] nvec = (c < 0) ? new int[R] : decode(p, Nvec);
            // arrivals
            for (int r = 0; r < R; r++) {
                if (nvec[r] < Nvec[r] && lam[r] > 0) {
                    int[] nnew = nvec.clone(); nnew[r]++;
                    int pn = linIndex(nnew, Nvec);
                    if (c < 0) Gq[s][key2id.get(stateKey(pn, r, 0, R))] += lam[r];
                    else Gq[s][key2id.get(stateKey(pn, c, k, R))] += lam[r];
                }
            }
            if (c < 0) continue;
            double rate = phRrow[c][k];
            double compl = rate * phProw[c][k];
            double adv = rate * (1.0 - phProw[c][k]);
            if (adv > 0 && k < Pc[c] - 1) Gq[s][key2id.get(stateKey(p, c, k + 1, R))] += adv;
            if (compl > 0) {
                int[] nnew = nvec.clone(); nnew[c]--;
                int tot = 0; for (int r = 0; r < R; r++) tot += nnew[r];
                if (tot == 0) {
                    Gq[s][0] += compl;
                } else {
                    int pn = linIndex(nnew, Nvec);
                    for (int cp = 0; cp < R; cp++) {
                        if (nnew[cp] > 0) Gq[s][key2id.get(stateKey(pn, cp, 0, R))] += compl * nnew[cp] / (double) tot;
                    }
                }
            }
        }
        for (int i = 0; i < S; i++) {
            double rowsum = 0.0;
            for (int j = 0; j < S; j++) rowsum += Gq[i][j];
            Gq[i][i] -= rowsum;
        }

        // Solve p*Gq=0, sum(p)=1: transpose, replace last row with normalization.
        Matrix A = new Matrix(S, S);
        for (int i = 0; i < S; i++) for (int j = 0; j < S; j++) A.set(i, j, Gq[j][i]);
        for (int j = 0; j < S; j++) A.set(S - 1, j, 1.0);
        Matrix b = new Matrix(S, 1); b.set(S - 1, 0, 1.0);
        Matrix pv = A.leftMatrixDivide(b);

        double[] Ppop = new double[npops];
        double[][] dep = new double[R][npops];
        for (int s = 0; s < S; s++) {
            int[] info = ids.get(s);
            int p = info[0], c = info[1], k = info[2];
            if (c < 0) continue;
            double ps = pv.get(s);
            Ppop[p] += ps;
            dep[c][p] += ps * phRrow[c][k] * phProw[c][k];
        }
        double[][] mumat = new double[R][npops];
        for (int r = 0; r < R; r++)
            for (int p = 0; p < npops; p++) mumat[r][p] = (Ppop[p] > 0) ? dep[r][p] / Ppop[p] : 0.0;
        return mumat;
    }

    private static long stateKey(int p, int c, int k, int R) {
        // p < npops, c < R, k < maxPhases(bounded, small); pack safely in a long.
        return (((long) p * (R + 1) + (c + 1)) << 8) + k;
    }
}
