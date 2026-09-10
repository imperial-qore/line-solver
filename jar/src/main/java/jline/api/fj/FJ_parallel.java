/**
 * @file Parallel-processing, team-service and independent-server models
 *
 * Ports of matlab/src/api/fj/fj_respt_nosplit.m, fj_respt_bulk.m,
 * fj_ism_green.m, fj_tsm_capacity.m, fj_serialization.m and fj_dag_makespan.m,
 * Sections 6 and 7 of A. Thomasian, "Analysis of Fork/Join and Related Queueing
 * Systems", ACM Computing Surveys 47(2), Article 17, 2014
 * (Eqs. (57)-(66) and Section 7.1).
 *
 * @since LINE 3.0
 */
package jline.api.fj;

import java.util.ArrayList;
import java.util.List;

public final class FJ_parallel {
    private FJ_parallel() {}

    /** [R, rho] of fj_respt_nosplit. */
    public static final class FJResptNosplitResult {
        public final double R;
        public final double rho;

        public FJResptNosplitResult(double R, double rho) {
            this.R = R;
            this.rho = rho;
        }
    }

    /**
     * Distributed no splitting: a job of K tasks goes in one piece to a single
     * server chosen uniformly among the K, which is an M/E_K/1 queue and reduces
     * to R = [ K - (K-1) rho/2 ] / (mu - lambda).
     */
    public static FJResptNosplitResult fj_respt_nosplit(int K, double lambda, double mu) {
        if (K < 1) {
            throw new IllegalArgumentException("K must be a positive integer. Got K=" + K + ".");
        }
        if (!(lambda > 0)) {
            throw new IllegalArgumentException("The arrival rate must be positive.");
        }
        if (!(mu > 0)) {
            throw new IllegalArgumentException("The service rate must be positive.");
        }
        double rho = lambda / mu;
        if (rho >= 1) {
            throw new IllegalArgumentException(
                    "System is unstable: rho = lambda/mu = " + String.format("%.4f", rho) + " >= 1.");
        }
        return new FJResptNosplitResult((K - (K - 1) * rho / 2) / (mu - lambda), rho);
    }

    /** [Rreq, Rtask, Q, p] of fj_respt_bulk. */
    public static final class FJResptBulkResult {
        public final double Rreq;
        public final double Rtask;
        public final double Q;
        public final double[] p;

        public FJResptBulkResult(double Rreq, double Rtask, double Q, double[] p) {
            this.Rreq = Rreq;
            this.Rtask = Rtask;
            this.Q = Q;
            this.p = p;
        }
    }

    /**
     * Centralized splitting as an M[K]/M/c bulk arrival system, solved by
     * truncating the level chain. The request response time is the completion of
     * the LAST of the K tasks: by PASTA the batch finds n tasks in system, its
     * last task is the (n+K)-th in line, and under first come first served with
     * c exponential servers it starts after max(0, n+K-c) departures, each an
     * exponential of rate c mu.
     */
    public static FJResptBulkResult fj_respt_bulk(int K, double lambda, double mu, int c,
                                                  int nmax) {
        if (K < 1) {
            throw new IllegalArgumentException("K must be a positive integer. Got K=" + K + ".");
        }
        if (c < 1) {
            throw new IllegalArgumentException("c must be a positive integer. Got c=" + c + ".");
        }
        if (!(lambda > 0) || !(mu > 0)) {
            throw new IllegalArgumentException("lambda and mu must be positive.");
        }
        double rho = lambda * K / (c * mu);
        if (rho >= 1) {
            throw new IllegalArgumentException(
                    "System is unstable: rho = lambda*K/(c*mu) = " + String.format("%.4f", rho)
                            + " >= 1.");
        }
        if (nmax <= 0) {
            nmax = (int) Math.max(200, Math.ceil(K + 40.0 * c / (1 - rho)));
        }
        int ns = nmax + 1;
        // Balance equations of the truncated chain, solved as a dense system with
        // the normalization replacing the last column
        double[][] A = new double[ns][ns + 1];
        for (int i = 0; i < ns; i++) {
            double out = 0;
            if (i > 0) {
                double srv = Math.min(i, c) * mu;
                // flow into i-1 from i
                A[i - 1][i] += srv;
                out += srv;
            }
            int j = i + K;
            if (j < ns) {
                A[j][i] += lambda;
                out += lambda;
            }
            A[i][i] -= out;
        }
        // Replace the last equation by the normalization
        for (int i = 0; i < ns; i++) {
            A[ns - 1][i] = 1;
        }
        A[ns - 1][ns] = 1;
        double[] p = gauss(A, ns);
        double sum = 0;
        for (int i = 0; i < ns; i++) {
            if (p[i] < 0) {
                p[i] = 0;
            }
            sum += p[i];
        }
        for (int i = 0; i < ns; i++) {
            p[i] /= sum;
        }
        double Q = 0;
        for (int n = 0; n < ns; n++) {
            Q += n * p[n];
        }
        double Rtask = Q / (lambda * K);
        double Rreq = 0;
        for (int n = 0; n < ns; n++) {
            int ahead = n + K;
            double wait = (ahead > c) ? (ahead - c) / (c * mu) : 0;
            Rreq += p[n] * (wait + 1 / mu);
        }
        return new FJResptBulkResult(Rreq, Rtask, Q, p);
    }

    public static FJResptBulkResult fj_respt_bulk(int K, double lambda, double mu, int c) {
        return fj_respt_bulk(K, lambda, mu, c, 0);
    }

    /** Everything Green's cycle decomposition produces. */
    public static final class FJIsmGreenResult {
        public final double W;
        public final double R;
        public final double EB;
        public final double EB2;
        public final double ED;
        public final double ED2;
        public final double EQ;
        public final double EQbar;
        public final double pq;
        public final double pd;
        public final double pi0;
        public final double rho;
        public final double ES;
        public final double[] q;

        public FJIsmGreenResult(double W, double R, double EB, double EB2, double ED, double ED2,
                                double EQ, double EQbar, double pq, double pd, double pi0,
                                double rho, double ES, double[] q) {
            this.W = W;
            this.R = R;
            this.EB = EB;
            this.EB2 = EB2;
            this.ED = ED;
            this.ED2 = ED2;
            this.EQ = EQ;
            this.EQbar = EQbar;
            this.pq = pq;
            this.pd = pd;
            this.pi0 = pi0;
            this.rho = rho;
            this.ES = ES;
            this.q = q;
        }
    }

    /**
     * Green's independent server model: a customer needs j servers at once with
     * probability c(j) and releases them asynchronously as each of its j tasks
     * completes at rate mu, so its own service is the maximum of j exponentials.
     *
     * E[B] is the j-th order statistic of s exponentials, because all s servers
     * are busy whenever a customer enters service during a queueing period; E[D]
     * is the initial delay of the customer that starts one; and the waiting-time
     * transform factors into the equilibrium transform of D and the
     * Pollaczek-Khinchine transform with service B, so the means add.
     *
     * Eq. (65) of the survey prints the inner sum of E[D] as starting at
     * 1/(s mu) even though only i servers are busy; it is started at 1/(i mu)
     * here, which is what the accompanying text prescribes.
     */
    public static FJIsmGreenResult fj_ism_green(double lambda, double mu, int s, double[] c) {
        if (s < 1) {
            throw new IllegalArgumentException("s must be a positive integer. Got s=" + s + ".");
        }
        if (c.length != s) {
            throw new IllegalArgumentException("c must have one entry per server requirement.");
        }
        if (!(lambda > 0) || !(mu > 0)) {
            throw new IllegalArgumentException("lambda and mu must be positive.");
        }
        double tot = 0;
        for (int i = 0; i < s; i++) {
            if (c[i] < 0) {
                throw new IllegalArgumentException(
                        "The server-requirement probabilities must be non-negative.");
            }
            tot += c[i];
        }
        if (Math.abs(tot - 1) > 1e-9) {
            throw new IllegalArgumentException(
                    "The server-requirement probabilities must sum to one. Got " + tot + ".");
        }
        double EB = 0, EB2 = 0;
        for (int j = 1; j <= s; j++) {
            double mj = 0, vj = 0;
            for (int i = 0; i < j; i++) {
                double st = 1.0 / ((s - i) * mu);
                mj += st;
                vj += st * st;
            }
            EB += c[j - 1] * mj;
            EB2 += c[j - 1] * (vj + mj * mj);
        }
        double rho = lambda * EB;
        if (rho >= 1) {
            throw new IllegalArgumentException(
                    "System is unstable: rho = lambda*E[B] = " + String.format("%.4f", rho) + " >= 1.");
        }
        int n = s + 1;
        double[][] T = new double[n][n];
        for (int i = 0; i < n; i++) {
            double den = lambda + i * mu;
            if (i > 0) {
                T[i][i - 1] = i * mu / den;
            }
            for (int j = 1; j + i <= s; j++) {
                T[i][i + j] += lambda * c[j - 1] / den;
            }
        }
        // v = row s of (I - T)^-1, obtained from (I-T)' x = e_s
        double[][] G = new double[n][n + 1];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                G[i][j] = (i == j ? 1 : 0) - T[j][i];
            }
            G[i][n] = (i == s) ? 1 : 0;
        }
        double[] v = gauss(G, n);
        double[] hold = new double[n];
        double EQbar = 0;
        for (int i = 0; i < n; i++) {
            hold[i] = 1.0 / (lambda + i * mu);
            EQbar += v[i] * hold[i];
        }
        double[] q = new double[n];
        for (int i = 0; i < n; i++) {
            q[i] = v[i] * hold[i] / EQbar;
        }
        double pd = 0;
        for (int i = 0; i < n; i++) {
            for (int j = s - i + 1; j <= s; j++) {
                pd += q[i] * c[j - 1];
            }
        }
        if (!(pd > 0)) {
            throw new IllegalArgumentException(
                    "No arrival can ever be delayed; the model degenerates to M/M/" + s + ".");
        }
        double ED = 0, ED2 = 0;
        for (int i = 1; i <= s; i++) {
            for (int k = 1; k <= i; k++) {
                int j = s - i + k;
                if (j < 1 || j > s) {
                    continue;
                }
                double wgt = q[i] * c[j - 1] / pd;
                if (wgt == 0) {
                    continue;
                }
                double mk = 0, vk = 0;
                for (int m = 0; m < k; m++) {
                    double st = 1.0 / ((i - m) * mu);
                    mk += st;
                    vk += st * st;
                }
                ED += wgt * mk;
                ED2 += wgt * (vk + mk * mk);
            }
        }
        double EQ = ED / (1 - rho);
        double pq = EQ / (EQ + EQbar);
        double pi0 = (1 - rho) / (1 - lambda * (EB - ED));
        double W = (1 - pi0) * (ED2 / (2 * ED) + lambda * EB2 / (2 * (1 - rho)));
        double ES = 0;
        for (int j = 1; j <= s; j++) {
            ES += c[j - 1] * FJ_harmonic.fj_harmonic(j) / mu;
        }
        return new FJIsmGreenResult(W, W + ES, EB, EB2, ED, ED2, EQ, EQbar, pq, pd, pi0, rho, ES, q);
    }

    /** [Lmax, Llp, Lfcfs, states, prob] of fj_tsm_capacity. */
    public static final class FJTsmCapacityResult {
        public final double Lmax;
        public final double Llp;
        public final double Lfcfs;
        public final boolean fcfsAvailable;
        public final int[][] states;
        public final double[] prob;

        public FJTsmCapacityResult(double Lmax, double Llp, double Lfcfs, boolean fcfsAvailable,
                                   int[][] states, double[] prob) {
            this.Lmax = Lmax;
            this.Llp = Llp;
            this.Lfcfs = Lfcfs;
            this.fcfsAvailable = fcfsAvailable;
            this.states = states;
            this.prob = prob;
        }
    }

    /**
     * Saturation throughput of the team service model. The apparent rate
     * Lambda_max = s / sum_k f(k) r(k) x(k) is attainable only when the
     * scheduler can pack jobs into execution states that leave no server idle;
     * the attainable capacity is the optimum of
     *
     * max Lambda s.t. sum_j p_j n(j,k)/x(k) = Lambda f(k), sum_j p_j = 1, p &gt;= 0
     *
     * over the multisets of jobs whose total server demand is at most s. For the
     * two-server two-class case with r = (1,2), strict first come first served
     * cannot pack at all and reaches only the printed lambda_FCFS.
     */
    public static FJTsmCapacityResult fj_tsm_capacity(int s, double[] f, int[] r, double[] x) {
        int K = f.length;
        if (s < 1) {
            throw new IllegalArgumentException("s must be a positive integer. Got s=" + s + ".");
        }
        if (r.length != K || x.length != K) {
            throw new IllegalArgumentException("f, r and x must have the same length.");
        }
        double tot = 0;
        for (int k = 0; k < K; k++) {
            if (f[k] < 0) {
                throw new IllegalArgumentException("The class frequencies must be non-negative.");
            }
            if (r[k] < 1 || r[k] > s) {
                throw new IllegalArgumentException("The server requirements must lie in 1.." + s + ".");
            }
            if (!(x[k] > 0)) {
                throw new IllegalArgumentException("The mean service times must be positive.");
            }
            tot += f[k];
        }
        if (Math.abs(tot - 1) > 1e-9) {
            throw new IllegalArgumentException(
                    "The class frequencies must sum to one. Got " + tot + ".");
        }
        double den = 0;
        for (int k = 0; k < K; k++) {
            den += f[k] * r[k] * x[k];
        }
        double Lmax = s / den;

        List<int[]> stateList = new ArrayList<int[]>();
        enumerateStates(r, s, 0, new int[K], stateList);
        int ns = stateList.size();
        if (ns == 0) {
            throw new IllegalArgumentException("No feasible execution state.");
        }
        int[][] states = stateList.toArray(new int[0][]);

        // Variables [p_1..p_ns, Lambda], minimise -Lambda
        List<Integer> active = new ArrayList<Integer>();
        for (int k = 0; k < K; k++) {
            if (f[k] > 0) {
                active.add(k);
            }
        }
        int rows = active.size() + 1;
        int cols = ns + 1;
        double[][] Aeq = new double[rows][cols];
        double[] beq = new double[rows];
        for (int a = 0; a < active.size(); a++) {
            int k = active.get(a);
            for (int j = 0; j < ns; j++) {
                Aeq[a][j] = states[j][k] / x[k];
            }
            Aeq[a][ns] = -f[k];
            beq[a] = 0;
        }
        for (int j = 0; j < ns; j++) {
            Aeq[rows - 1][j] = 1;
        }
        beq[rows - 1] = 1;
        double[] cost = new double[cols];
        cost[ns] = -1;
        double[] sol = simplex(Aeq, beq, cost);
        double Llp = sol[ns];
        double[] prob = new double[ns];
        System.arraycopy(sol, 0, prob, 0, ns);

        double Lfcfs = Double.NaN;
        boolean have = false;
        if (s == 2 && K == 2 && ((r[0] == 1 && r[1] == 2) || (r[0] == 2 && r[1] == 1))) {
            int a = (r[0] == 1) ? 0 : 1;
            int b = 1 - a;
            double f1 = f[a], f2 = f[b], mu1 = 1 / x[a], mu2 = 1 / x[b];
            Lfcfs = 2 * mu1 * mu2
                    / (f1 * f1 * mu2 + 2 * f2 * f2 * mu1 + 2 * f1 * f2 * (mu1 + mu2));
            have = true;
        }
        return new FJTsmCapacityResult(Lmax, Llp, Lfcfs, have, states, prob);
    }

    private static void enumerateStates(int[] r, int left, int k, int[] stack, List<int[]> out) {
        if (k == r.length) {
            for (int i = 0; i < stack.length; i++) {
                if (stack[i] > 0) {
                    out.add(stack.clone());
                    return;
                }
            }
            return;
        }
        int nmax = left / r[k];
        for (int n = 0; n <= nmax; n++) {
            stack[k] = n;
            enumerateStates(r, left - n * r[k], k + 1, stack, out);
        }
        stack[k] = 0;
    }

    /** [P, delay, Rtot] of fj_serialization. */
    public static final class FJSerializationResult {
        public final double[] P;
        public final double[] delay;
        public final double Rtot;

        public FJSerializationResult(double[] P, double[] delay, double Rtot) {
            this.P = P;
            this.delay = delay;
            this.Rtot = Rtot;
        }
    }

    /**
     * Blocking probability and pseudoserver delay of serialization phases:
     * P_s(M) = 1 - [ 1 - R_s(M)/R(M) ]^(M-1), with the delay charged at the
     * pseudoserver equal to alpha R_s(M).
     */
    public static FJSerializationResult fj_serialization(double[] Rs, double R0, int M,
                                                         double alpha) {
        int S = Rs.length;
        if (S < 1) {
            throw new IllegalArgumentException("At least one serialization phase is required.");
        }
        for (int i = 0; i < S; i++) {
            if (Rs[i] < 0) {
                throw new IllegalArgumentException(
                        "The phase residence times must be non-negative.");
            }
        }
        if (R0 < 0) {
            throw new IllegalArgumentException(
                    "The nonserialized residence time must be non-negative.");
        }
        if (M < 1) {
            throw new IllegalArgumentException("M must be a positive integer. Got M=" + M + ".");
        }
        if (alpha < 0 || alpha > 1) {
            throw new IllegalArgumentException("alpha must lie in [0,1]. Got " + alpha + ".");
        }
        double R = R0;
        for (int i = 0; i < S; i++) {
            R += Rs[i];
        }
        if (!(R > 0)) {
            throw new IllegalArgumentException("The total residence time vanished.");
        }
        double[] P = new double[S];
        double[] delay = new double[S];
        double Rtot = R;
        for (int i = 0; i < S; i++) {
            P[i] = 1 - Math.pow(1 - Rs[i] / R, M - 1);
            delay[i] = P[i] * alpha * Rs[i];
            Rtot += delay[i];
        }
        return new FJSerializationResult(P, delay, Rtot);
    }

    public static FJSerializationResult fj_serialization(double[] Rs, double R0, int M) {
        return fj_serialization(Rs, R0, M, 0.5);
    }

    /** [C, I, Cend, E] of fj_dag_makespan. */
    public static final class FJDagMakespanResult {
        public final double C;
        public final double[] I;
        public final double[] Cend;
        public final double[] E;

        public FJDagMakespanResult(double C, double[] I, double[] Cend, double[] E) {
            this.C = C;
            this.I = I;
            this.Cend = Cend;
            this.E = E;
        }
    }

    /**
     * Makespan of a task system with precedence constraints. The chain whose
     * state is the SET of completed tasks is acyclic, so it is swept level by
     * level: task i among the k eligible completes at rate rate[i][k-1], the
     * state is held for 1/T(S), and
     *
     * p(R) = sum p(S) b(S,R),  D(R) = M(R) p(R) + sum b(S,R) D(S),
     *
     * started at p(empty) = 1. Making the rate depend on the concurrency is what
     * couples the task system to the queueing network underneath it.
     */
    public static FJDagMakespanResult fj_dag_makespan(boolean[][] pred, double[][] rate) {
        int n = pred.length;
        if (n < 1) {
            throw new IllegalArgumentException("At least one task is required.");
        }
        if (n > 20) {
            throw new IllegalArgumentException("The completed-set sweep enumerates 2^n states.");
        }
        for (int i = 0; i < n; i++) {
            if (pred[i].length != n) {
                throw new IllegalArgumentException("pred must be square.");
            }
            if (rate[i].length != n) {
                throw new IllegalArgumentException("The rate table must be n by n.");
            }
            for (int k = 0; k < n; k++) {
                if (!(rate[i][k] > 0)) {
                    throw new IllegalArgumentException("All completion rates must be positive.");
                }
            }
        }
        int[] predmask = new int[n];
        int[] indeg = new int[n];
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < n; i++) {
                if (pred[i][j]) {
                    predmask[j] |= (1 << i);
                    indeg[j]++;
                }
            }
        }
        boolean[] seen = new boolean[n];
        int remaining = n;
        for (int pass = 0; pass < n; pass++) {
            int pick = -1;
            for (int i = 0; i < n; i++) {
                if (!seen[i] && indeg[i] == 0) {
                    pick = i;
                    break;
                }
            }
            if (pick < 0) {
                break;
            }
            seen[pick] = true;
            indeg[pick] = -1;
            for (int j = 0; j < n; j++) {
                if (pred[pick][j] && indeg[j] > 0) {
                    indeg[j]--;
                }
            }
            remaining--;
        }
        if (remaining > 0) {
            throw new IllegalArgumentException("The precedence relation contains a cycle.");
        }
        int nmask = 1 << n;
        int[] eligmask = new int[nmask];
        boolean[] closed = new boolean[nmask];
        for (int mask = 0; mask < nmask; mask++) {
            boolean ok = true;
            int em = 0;
            for (int i = 0; i < n; i++) {
                int bit = 1 << i;
                if ((mask & bit) != 0) {
                    if ((predmask[i] & mask) != predmask[i]) {
                        ok = false;
                        break;
                    }
                } else if ((predmask[i] & mask) == predmask[i]) {
                    em |= bit;
                }
            }
            closed[mask] = ok;
            if (ok) {
                eligmask[mask] = em;
            }
        }
        double[] p = new double[nmask];
        double[] D = new double[nmask];
        double[] I = new double[n];
        double[] Cend = new double[n];
        p[0] = 1;
        for (int mask = 0; mask < nmask; mask++) {
            if (!closed[mask]) {
                continue;
            }
            int em = eligmask[mask];
            int k = Integer.bitCount(em);
            if (k == 0) {
                continue;
            }
            double Ttot = 0;
            for (int i = 0; i < n; i++) {
                if ((em & (1 << i)) != 0) {
                    Ttot += rate[i][k - 1];
                }
            }
            D[mask] += p[mask] / Ttot;
            for (int i = 0; i < n; i++) {
                if ((em & (1 << i)) == 0) {
                    continue;
                }
                double b = rate[i][k - 1] / Ttot;
                int nxt = mask | (1 << i);
                double contrib = b * D[mask];
                p[nxt] += p[mask] * b;
                D[nxt] += contrib;
                Cend[i] += contrib;
                int fresh = eligmask[nxt] & ~em;
                for (int j = 0; j < n; j++) {
                    if ((fresh & (1 << j)) != 0) {
                        I[j] += contrib;
                    }
                }
            }
        }
        double[] E = new double[n];
        for (int i = 0; i < n; i++) {
            E[i] = Cend[i] - I[i];
        }
        return new FJDagMakespanResult(D[nmask - 1], I, Cend, E);
    }

    /** Gaussian elimination with partial pivoting on an augmented n by n+1 system. */
    private static double[] gauss(double[][] A, int n) {
        for (int col = 0; col < n; col++) {
            int piv = col;
            for (int r = col + 1; r < n; r++) {
                if (Math.abs(A[r][col]) > Math.abs(A[piv][col])) {
                    piv = r;
                }
            }
            if (A[piv][col] == 0) {
                throw new IllegalArgumentException("The linear system is singular.");
            }
            double[] tmp = A[col];
            A[col] = A[piv];
            A[piv] = tmp;
            for (int r = 0; r < n; r++) {
                if (r == col) {
                    continue;
                }
                double f = A[r][col] / A[col][col];
                for (int j = col; j <= n; j++) {
                    A[r][j] -= f * A[col][j];
                }
            }
        }
        double[] xs = new double[n];
        for (int i = 0; i < n; i++) {
            xs[i] = A[i][n] / A[i][i];
        }
        return xs;
    }

    /**
     * Two-phase dense simplex with Bland's rule on min c'z, A z = b, z &gt;= 0.
     *
     * Bland's rule terminates without cycling and without an anti-cycling
     * perturbation, at the cost of a slower pivot sequence; the capacity
     * programs here have at most a few hundred columns, so the dense tableau is
     * the right trade and keeps this class free of any solver dependency.
     */
    private static double[] simplex(double[][] Ain, double[] bin, double[] cin) {
        int m = Ain.length;
        int n = cin.length;
        double tol = 1e-10;
        double[][] A = new double[m][n];
        double[] b = new double[m];
        for (int i = 0; i < m; i++) {
            System.arraycopy(Ain[i], 0, A[i], 0, n);
            b[i] = bin[i];
            if (b[i] < 0) {
                for (int j = 0; j < n; j++) {
                    A[i][j] = -A[i][j];
                }
                b[i] = -b[i];
            }
        }
        int width = n + m + 1;
        double[][] T = new double[m][width];
        int[] basis = new int[m];
        for (int i = 0; i < m; i++) {
            System.arraycopy(A[i], 0, T[i], 0, n);
            T[i][n + i] = 1;
            T[i][width - 1] = b[i];
            basis[i] = n + i;
        }
        double[] obj = new double[width];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                obj[j] -= T[i][j];
            }
            obj[width - 1] -= T[i][width - 1];
        }
        pivotLoop(T, basis, obj, tol);
        if (-obj[width - 1] > 1e-7) {
            throw new IllegalArgumentException("The linear program is infeasible.");
        }
        // Drive the artificials out of the basis
        for (int row = 0; row < T.length; row++) {
            if (basis[row] >= n) {
                int piv = -1;
                for (int j = 0; j < n; j++) {
                    if (Math.abs(T[row][j]) > tol) {
                        piv = j;
                        break;
                    }
                }
                if (piv >= 0) {
                    doPivot(T, row, piv);
                    basis[row] = piv;
                }
            }
        }
        // Phase two on the original columns
        double[][] T2 = new double[m][n + 1];
        for (int i = 0; i < m; i++) {
            System.arraycopy(T[i], 0, T2[i], 0, n);
            T2[i][n] = T[i][width - 1];
        }
        double[] obj2 = new double[n + 1];
        System.arraycopy(cin, 0, obj2, 0, n);
        for (int i = 0; i < m; i++) {
            if (basis[i] < n) {
                double f = obj2[basis[i]];
                for (int j = 0; j <= n; j++) {
                    obj2[j] -= f * T2[i][j];
                }
            }
        }
        pivotLoop(T2, basis, obj2, tol);
        double[] z = new double[n];
        for (int i = 0; i < m; i++) {
            if (basis[i] < n) {
                z[basis[i]] = T2[i][n];
            }
        }
        return z;
    }

    private static void pivotLoop(double[][] T, int[] basis, double[] obj, double tol) {
        int ncol = T[0].length - 1;
        int maxit = 200 * (T.length + ncol) + 1000;
        for (int it = 0; it < maxit; it++) {
            int enter = -1;
            for (int j = 0; j < ncol; j++) {
                if (obj[j] < -tol) {
                    enter = j;
                    break;
                }
            }
            if (enter < 0) {
                return;
            }
            int leave = -1;
            double best = Double.POSITIVE_INFINITY;
            for (int i = 0; i < T.length; i++) {
                if (T[i][enter] > tol) {
                    double ratio = T[i][T[i].length - 1] / T[i][enter];
                    if (ratio < best - tol
                            || (Math.abs(ratio - best) <= tol && (leave < 0 || basis[i] < basis[leave]))) {
                        best = ratio;
                        leave = i;
                    }
                }
            }
            if (leave < 0) {
                throw new IllegalArgumentException("The linear program is unbounded.");
            }
            doPivot(T, leave, enter);
            double f = obj[enter];
            for (int j = 0; j < T[leave].length; j++) {
                obj[j] -= f * T[leave][j];
            }
            basis[leave] = enter;
        }
        throw new IllegalArgumentException("The simplex hit its iteration limit.");
    }

    private static void doPivot(double[][] T, int row, int col) {
        double piv = T[row][col];
        for (int j = 0; j < T[row].length; j++) {
            T[row][j] /= piv;
        }
        for (int i = 0; i < T.length; i++) {
            if (i != row && T[i][col] != 0) {
                double f = T[i][col];
                for (int j = 0; j < T[i].length; j++) {
                    T[i][j] -= f * T[row][j];
                }
            }
        }
    }
}
