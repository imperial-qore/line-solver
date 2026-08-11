package jline.api.pfqn;

import org.apache.commons.math3.random.MersenneTwister;

import jline.util.matrix.Matrix;
import jline.util.RandomManager;

/**
 * Perfect/approximate stationary state sampler for closed single-class
 * multiserver product-form networks.
 *
 * <p>Draws states exactly distributed according to the product-form stationary
 * distribution of a closed single-class Jackson network with multiple servers,
 * using monotone Coupling From The Past (Propp-Wilson) as proposed by Kijima
 * and Matsui.</p>
 *
 * <p>Reference: S. Kijima and T. Matsui, "Approximate/Perfect Samplers for
 * Closed Jackson Networks", Proc. Winter Simulation Conference, 2005.</p>
 *
 * @since LINE 3.0
 */
public final class Pfqn_cftp {

    private Pfqn_cftp() {}

    /**
     * Result holder mirroring the MATLAB [Q,X,T] return of pfqn_cftp.
     */
    public static final class PfqnCftpReturn {
        /** Empirical mean queue length per station (1 x M). */
        public final Matrix Q;
        /** Sampled states, one per row (nsamples x M), each row sums to N. */
        public final Matrix X;
        /** Per-sample coalescence horizon ('cftp') or mixing steps ('approx') (nsamples x 1). */
        public final Matrix T;

        public PfqnCftpReturn(Matrix Q, Matrix X, Matrix T) {
            this.Q = Q;
            this.X = X;
            this.T = T;
        }
    }

    /**
     * Exact (perfect) stationary state sampling with default single servers,
     * one sample, and the 'cftp' method.
     *
     * @param L service demands (M entries), L(i) = theta_i/mu_i
     * @param N total closed population K
     * @return the [Q,X,T] result
     */
    public static PfqnCftpReturn pfqn_cftp(Matrix L, int N) {
        return pfqn_cftp(L, N, null, 1, "cftp");
    }

    /**
     * Exact (perfect) or approximate stationary state sampling for closed
     * single-class multiserver product-form networks via monotone Coupling From
     * The Past.
     *
     * @param L        service demands (M entries), L(i) = theta_i/mu_i
     * @param N        total closed population K
     * @param S        servers per station (M entries), Double.POSITIVE_INFINITY
     *                 for infinite server; null defaults to one server each
     * @param nsamples number of independent samples to draw
     * @param method   'cftp' (exact, default) or 'approx' (rapidly mixing M_A)
     * @return the [Q,X,T] result
     */
    public static PfqnCftpReturn pfqn_cftp(Matrix L, int N, Matrix S, int nsamples, String method) {
        double[] Lv = flatten(L);
        int M = Lv.length;

        double[] Sv;
        if (S == null) {
            Sv = new double[M];
            for (int i = 0; i < M; i++) {
                Sv[i] = 1.0;
            }
        } else {
            Sv = flatten(S);
        }

        int K = (int) Math.round((double) N);

        if (M < 2) {
            throw new RuntimeException("pfqn_cftp: at least two stations are required.");
        }
        for (int i = 0; i < M; i++) {
            if (Lv[i] <= 0) {
                throw new RuntimeException("pfqn_cftp: all demands L must be strictly positive.");
            }
        }
        if (method == null || method.length() == 0) {
            method = "cftp";
        }

        // see _kb/03-api-layer.md for rationale
        double[][] logfac = new double[M][K + 1];
        for (int i = 0; i < M; i++) {
            double acc = 0.0;
            for (int m = 1; m <= K; m++) {
                acc += Math.log(Math.min(m, Sv[i]));
                logfac[i][m] = acc;
            }
        }
        double[] logL = new double[M];
        for (int i = 0; i < M; i++) {
            logL[i] = Math.log(Lv[i]);
        }

        MersenneTwister rand = RandomManager.getThreadRandom();

        Matrix X = new Matrix(nsamples, M);
        Matrix T = new Matrix(nsamples, 1);
        String m = method.toLowerCase();
        for (int smp = 0; smp < nsamples; smp++) {
            int[] x;
            int horizon;
            if (m.equals("cftp")) {
                int[] out = new int[M + 1];
                horizon = drawCftp(logL, logfac, K, M, rand, out);
                x = out;
            } else if (m.equals("approx")) {
                int[] out = new int[M + 1];
                horizon = drawApprox(logL, logfac, K, M, rand, out);
                x = out;
            } else {
                throw new RuntimeException("pfqn_cftp: unknown method '" + method + "'.");
            }
            for (int i = 0; i < M; i++) {
                X.set(smp, i, x[i]);
            }
            T.set(smp, 0, horizon);
        }

        // Q = column mean of X (1 x M).
        Matrix Q = new Matrix(1, M);
        for (int i = 0; i < M; i++) {
            double s = 0.0;
            for (int smp = 0; smp < nsamples; smp++) {
                s += X.get(smp, i);
            }
            Q.set(0, i, s / nsamples);
        }

        return new PfqnCftpReturn(Q, X, T);
    }

    /**
     * One exact draw via monotone CFTP. Writes the coalesced state into out[0..M-1]
     * and returns the coalescence horizon.
     */
    private static int drawCftp(double[] logL, double[][] logfac, int K, int M,
                                MersenneTwister rand, int[] out) {
        // top state x_U = (K,0,...,0), bottom state x_L = (0,...,0,K)
        double[] u = new double[0];
        int Tback = 1;
        while (true) {
            int old = u.length;
            // prepend randomness for the newly exposed (older) steps; older steps
            // sit at the front and are applied first, recent steps are reused.
            double[] nu = new double[Tback];
            for (int i = 0; i < Tback - old; i++) {
                nu[i] = rand.nextDouble();
            }
            System.arraycopy(u, 0, nu, Tback - old, old);
            u = nu;

            int[] xU = new int[M];
            xU[0] = K;
            int[] xL = new int[M];
            xL[M - 1] = K;
            for (int t = 0; t < Tback; t++) {        // oldest step (front) first
                monotoneUpdate(xU, u[t], logL, logfac, M);
                monotoneUpdate(xL, u[t], logL, logfac, M);
            }
            boolean coalesced = true;
            for (int i = 0; i < M; i++) {
                if (xU[i] != xL[i]) {
                    coalesced = false;
                    break;
                }
            }
            if (coalesced) {
                System.arraycopy(xU, 0, out, 0, M);
                return Tback;
            }
            Tback *= 2;
        }
    }

    /**
     * Approximate rapidly-mixing sampler M_A. Writes the state into out[0..M-1]
     * and returns the number of mixing steps used.
     */
    private static int drawApprox(double[] logL, double[][] logfac, int K, int M,
                                  MersenneTwister rand, int[] out) {
        double eps = 1e-2;
        int steps = (int) Math.ceil(M * (M - 1) / 2.0 * Math.log(K / eps)); // mixing-time bound
        int[] x = new int[M];
        x[0] = K;                                    // arbitrary feasible start
        for (int t = 0; t < steps; t++) {
            int i = rand.nextInt(M);
            int j = rand.nextInt(M - 1);             // distinct pair (not necessarily adjacent)
            if (j >= i) {
                j++;
            }
            int k = x[i] + x[j];
            int l = splitIndex(logL, logfac, i, j, k, rand.nextDouble());
            x[i] = l;
            x[j] = k - l;
        }
        System.arraycopy(x, 0, out, 0, M);
        return steps;
    }

    /**
     * Monotone update on a consecutive pair. A single uniform u in [0,1) encodes
     * the pair index (integer part) and the split Lambda (fractional part).
     */
    private static void monotoneUpdate(int[] x, double u, double[] logL, double[][] logfac, int M) {
        double lam = 1 + u * (M - 1);                // in [1, M)
        int jj = (int) Math.floor(lam);              // 1-based pair index (jj, jj+1), 1 <= jj <= M-1
        if (jj > M - 1) {
            jj = M - 1;
        }
        double Lambda = lam - jj;                    // fractional part in [0, 1)
        int p = jj - 1;                              // 0-based left station of the adjacent pair
        int k = x[p] + x[p + 1];
        int l = splitIndex(logL, logfac, p, p + 1, k, Lambda);
        x[p] = l;
        x[p + 1] = k - l;
    }

    /**
     * Inverse-CDF split: smallest l with Lambda &lt;= g^k_{ij}(l), where the
     * split weight w(s) is proportional to alpha_i(s)*alpha_j(k-s), s = 0..k.
     */
    private static int splitIndex(double[] logL, double[][] logfac, int i, int j, int k, double Lambda) {
        double[] lw = new double[k + 1];
        double mx = Double.NEGATIVE_INFINITY;
        for (int s = 0; s <= k; s++) {
            lw[s] = (s * logL[i] - logfac[i][s]) + ((k - s) * logL[j] - logfac[j][k - s]);
            if (lw[s] > mx) {
                mx = lw[s];
            }
        }
        double cum = 0.0;
        double[] cdf = new double[k + 1];
        for (int s = 0; s <= k; s++) {
            cum += Math.exp(lw[s] - mx);
            cdf[s] = cum;
        }
        double total = cdf[k];
        for (int s = 0; s <= k; s++) {
            if (Lambda <= cdf[s] / total) {
                return s;
            }
        }
        return k;
    }

    /** Flatten a Matrix (row/column vector) into a double[] in element order. */
    private static double[] flatten(Matrix v) {
        int n = v.getNumRows() * v.getNumCols();
        double[] out = new double[n];
        int idx = 0;
        for (int r = 0; r < v.getNumRows(); r++) {
            for (int c = 0; c < v.getNumCols(); c++) {
                out[idx++] = v.get(r, c);
            }
        }
        return out;
    }
}
