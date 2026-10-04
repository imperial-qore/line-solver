/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.mc;

import java.util.List;

import org.apache.commons.math3.complex.Complex;

import jline.io.InputOutput;
import jline.util.matrix.Matrix;

/**
 * Saddlepoint approximation of Pr{N(t)=k} for the counting process of a MAP.
 *
 * <p>The probability that the Markovian arrival process (D0,D1) records exactly
 * k events in (0,t], obtained by steepest-descent inversion of the counting
 * generating function instead of by forming the k-th superdiagonal block of
 * expm(t*X).
 *
 * <p>The counting generating function is the matrix exponential
 *
 * <pre>
 *     sum_k P(k,t) z^k = expm(t*(D0 + z*D1)),
 * </pre>
 *
 * so the cumulant generating function of N(t) is eta(theta) = spectral abscissa
 * of A(theta) = D0 + exp(theta)*D1, the Perron root of an irreducible Metzler
 * matrix: real, simple, strictly convex in theta, with eta(0)=0 and
 * eta'(0)=lambda. Inverting by steepest descent gives Daniels (1954),
 *
 * <pre>
 *     Pr{N(t)=k} ~ g(theta*) * exp(t*eta(theta*) - k*theta*)
 *                            / sqrt(2*pi*t*eta''(theta*)),
 * </pre>
 *
 * with the saddle theta* solving eta'(theta*) = k/t and g the amplitude of the
 * Perron projection, g(theta) = (pi0*v)*(u*1), u and v the left and right
 * Perron vectors normalised by u*v = 1.
 *
 * <p>THE EXPANSION PARAMETER IS K2 = t*eta''(theta*), THE VARIANCE OF THE
 * COUNT, not its mean and not t. Measured error laws, with the constants flat
 * to two digits over Erlang orders 1..8 and horizons 10..160:
 *
 * <pre>
 *     err(DANIELS) = 0.083 / K2        err(DANIELS2) = 0.017 / K2^2
 * </pre>
 *
 * For a renewal Erlang(r) the count variance rate is lambda/r, so
 * K2 = lambda*t/r and an Erlang-4 at t=50 is as accurate as a Poisson at
 * t=12.5: low variability shrinks the parameter, it does not break the method.
 * Below K2 = 5 the expansion is out of its regime and the call warns.
 *
 * <p>This is an asymptotic method, not a quadrature: use it for rare-event and
 * large-deviation coefficients, where k/t is away from lambda or where the
 * probability underflows. For the bulk of the transient distribution, i.e.
 * every block k=0..N-1 at once at moderate t, uniformization
 * ({@link Ctmc_uniformization}, {@link Ctmc_foxglynn}) is both exact and faster.
 *
 * <p>ATTRIBUTION. The first-order form is Daniels (1954). The amplitude g and the
 * whole 'daniels2' bracket are NOT a rederivation: they are Jensen, "Saddlepoint
 * Expansions for Sums of Markov Dependent Variables on a Continuous State Space",
 * Probab. Th. Rel. Fields 89, 1991, Eq. (4.4) with the coefficients on p.191. His
 * gamma_0(s) = (sum_i c_i)(sum_i r_i P(Y_0=i)) is exactly g under his own
 * normalisation sum_i r_i c_i = 1, and expanding his
 * alpha_0 + (1/n){-alpha_3/2 + alpha_4/8 - 5*alpha_5/24} reproduces
 * g*(1 + lam4/8 - 5*lam3^2/24) - g''/(2*K2) + g'*K3/(2*K2^2) term for term; his
 * Theorem 4.1 gives the O(n^-2) error measured here as 0.017/K2^2. Jensen works
 * with discrete-n sums over a Markov chain, so the continuous-time MAP counting
 * process is that result transcribed, n -> t and the kernel eigenvalue -> the
 * Perron root of D0+exp(theta)*D1.
 *
 * @since LINE 3.0
 */
public final class Ctmc_saddlepoint {

    private Ctmc_saddlepoint() {
    }

    /**
     * Below this value of K2 = t*eta''(theta*) the expansion is out of its
     * regime. Do NOT threshold on lambda*t: for Erlang(r) the count variance
     * rate is lambda/r, so K2 = lambda*t/r, and lambda*t over-warns on
     * Poisson-like processes while under-warning on low-variability ones.
     */
    public static final double K2_MIN = 5.0;

    /** Second-order saddlepoint, error O(1/K2^2). The default. */
    public static final String DANIELS2 = "daniels2";
    /** First-order saddlepoint with the Perron amplitude, error O(1/K2). */
    public static final String DANIELS = "daniels";
    /** Bare first-order form with the amplitude set to 1. */
    public static final String PLAIN = "plain";

    /**
     * Perron root of A(theta) with its first two derivatives in theta and the
     * amplitude of the Perron projection between pi0 and 1.
     */
    public static final class PerronState {
        public final double eta;
        public final double deta;
        public final double d2eta;
        public final double ampl;

        PerronState(double eta, double deta, double d2eta, double ampl) {
            this.eta = eta;
            this.deta = deta;
            this.d2eta = d2eta;
            this.ampl = ampl;
        }
    }

    /**
     * Result of a saddlepoint evaluation, one entry per (t,k) pair.
     */
    public static final class SaddlepointResult {
        /** Approximation of Pr{N(t)=k}. */
        public final double[] p;
        /** Its natural logarithm, accurate below the smallest positive double. */
        public final double[] logp;
        /** The saddle theta*, -Infinity where k=0. */
        public final double[] theta;
        /** eta(theta*). */
        public final double[] eta;
        /** eta'(theta*), equal to k/t at convergence. */
        public final double[] deta;
        /** eta''(theta*). */
        public final double[] d2eta;
        /** The Perron amplitude g(theta*). */
        public final double[] ampl;
        /** The bracket multiplying the leading term. */
        public final double[] corr;
        /** K2 = t*eta''(theta*), the expansion parameter. */
        public final double[] k2;
        /** Newton steps taken. */
        public final int[] iter;
        /** True where the value was computed exactly rather than approximated. */
        public final boolean[] exact;
        /** The stationary event rate eta'(0). */
        public final double lambda;

        SaddlepointResult(int n, double lambda) {
            this.p = new double[n];
            this.logp = new double[n];
            this.theta = new double[n];
            this.eta = new double[n];
            this.deta = new double[n];
            this.d2eta = new double[n];
            this.ampl = new double[n];
            this.corr = new double[n];
            this.k2 = new double[n];
            this.iter = new int[n];
            this.exact = new boolean[n];
            this.lambda = lambda;
        }
    }

    /**
     * Pr{N(t)=k} at a single (t,k), with the default method.
     *
     * @param D0 generator of the phase process with the counted transitions removed
     * @param D1 rates of the counted transitions; D0+D1 must be an irreducible generator
     * @param t  time horizon
     * @param k  event count, a nonnegative integer
     * @return the approximation
     */
    public static double ctmc_saddlepoint(Matrix D0, Matrix D1, double t, int k) {
        return ctmc_saddlepoint(D0, D1, new double[]{t}, new int[]{k}, DANIELS2, null).p[0];
    }

    /**
     * Pr{N(t)=k} at a single (t,k).
     *
     * @param D0     generator of the phase process with the counted transitions removed
     * @param D1     rates of the counted transitions
     * @param t      time horizon
     * @param k      event count
     * @param method {@link #DANIELS2} (default), {@link #DANIELS} or {@link #PLAIN}
     * @return the approximation
     */
    public static double ctmc_saddlepoint(Matrix D0, Matrix D1, double t, int k, String method) {
        return ctmc_saddlepoint(D0, D1, new double[]{t}, new int[]{k}, method, null).p[0];
    }

    /**
     * Pr{N(t)=k} over arrays of horizons and counts.
     *
     * @param D0     generator of the phase process with the counted transitions removed
     * @param D1     rates of the counted transitions
     * @param t      time horizons; length 1 broadcasts against k
     * @param k      event counts; length 1 broadcasts against t
     * @param method {@link #DANIELS2} (default), {@link #DANIELS} or {@link #PLAIN};
     *               null selects the default
     * @param pi0    initial phase distribution; null selects the stationary
     *               distribution of D0+D1
     * @return the per-point approximations and diagnostics
     */
    public static SaddlepointResult ctmc_saddlepoint(Matrix D0, Matrix D1, double[] t, int[] k,
                                                     String method, double[] pi0) {
        final int nph = D0.getNumRows();
        if (D0.getNumCols() != nph || D1.getNumRows() != nph || D1.getNumCols() != nph) {
            throw new IllegalArgumentException(
                    "ctmc_saddlepoint: D0 and D1 must be square matrices of the same order.");
        }
        double maxrate = 0.0;
        for (int i = 0; i < nph; i++) {
            for (int j = 0; j < nph; j++) {
                double d1ij = D1.get(i, j);
                if (d1ij < 0) {
                    throw new IllegalArgumentException("ctmc_saddlepoint: D1 must be nonnegative.");
                }
                if (d1ij > maxrate) {
                    maxrate = d1ij;
                }
            }
        }
        if (maxrate <= 0.0) {
            throw new IllegalArgumentException("ctmc_saddlepoint: D1 has no counted transitions, "
                    + "the counting process is identically zero.");
        }
        double[][] Q = new double[nph][nph];
        double maxq = 0.0;
        for (int i = 0; i < nph; i++) {
            for (int j = 0; j < nph; j++) {
                Q[i][j] = D0.get(i, j) + D1.get(i, j);
                maxq = Math.max(maxq, Math.abs(Q[i][j]));
            }
        }
        for (int i = 0; i < nph; i++) {
            double rowsum = 0.0;
            for (int j = 0; j < nph; j++) {
                rowsum += Q[i][j];
            }
            if (Math.abs(rowsum) > 1e-8 * Math.max(1.0, maxq)) {
                throw new IllegalArgumentException("ctmc_saddlepoint: D0+D1 must be an "
                        + "infinitesimal generator (zero row sums).");
            }
        }

        String key = (method == null) ? DANIELS2 : method.trim().toLowerCase();
        int order;
        boolean useampl;
        if (DANIELS2.equals(key) || "sp2".equals(key)) {
            order = 2;
            useampl = true;
        } else if (DANIELS.equals(key) || "sp1".equals(key)) {
            order = 1;
            useampl = true;
        } else if (PLAIN.equals(key) || "bare".equals(key)) {
            order = 1;
            useampl = false;
        } else {
            throw new IllegalArgumentException("ctmc_saddlepoint: unknown method '" + key
                    + "', expected daniels2, daniels or plain.");
        }

        double[] pi = pi0;
        if (pi == null) {
            Matrix Qm = new Matrix(nph, nph);
            for (int i = 0; i < nph; i++) {
                for (int j = 0; j < nph; j++) {
                    Qm.set(i, j, Q[i][j]);
                }
            }
            Matrix piRow = Ctmc_solve.ctmc_solve(Qm);
            pi = new double[nph];
            for (int i = 0; i < nph; i++) {
                pi[i] = piRow.get(i);           // row or column vector, read linearly
            }
        }
        if (pi.length != nph) {
            throw new IllegalArgumentException("ctmc_saddlepoint: pi0 must have one entry per phase.");
        }
        double pisum = 0.0;
        for (int i = 0; i < nph; i++) {
            pisum += pi[i];
        }
        if (Math.abs(pisum - 1.0) > 1e-8) {
            throw new IllegalArgumentException("ctmc_saddlepoint: pi0 must sum to one.");
        }

        // Broadcast the horizons against the counts
        int n = Math.max(t.length, k.length);
        if (t.length != n && t.length != 1) {
            throw new IllegalArgumentException(
                    "ctmc_saddlepoint: t and k must be scalars or arrays of the same size.");
        }
        if (k.length != n && k.length != 1) {
            throw new IllegalArgumentException(
                    "ctmc_saddlepoint: t and k must be scalars or arrays of the same size.");
        }
        double[] tv = new double[n];
        int[] kv = new int[n];
        for (int i = 0; i < n; i++) {
            tv[i] = t.length == 1 ? t[0] : t[i];
            kv[i] = k.length == 1 ? k[0] : k[i];
            if (tv[i] < 0) {
                throw new IllegalArgumentException("ctmc_saddlepoint: the horizon t must be nonnegative.");
            }
            if (kv[i] < 0) {
                throw new IllegalArgumentException("ctmc_saddlepoint: the count k must be a nonnegative integer.");
            }
        }

        double[][] d0 = toArray(D0, nph);
        double[][] d1 = toArray(D1, nph);
        double[] ones = new double[nph];
        for (int i = 0; i < nph; i++) {
            ones[i] = 1.0;
        }

        // exp(theta) multiplies D1, so the saddle is confined to the range over
        // which A(theta) is representable; never active for a feasible k/t
        final double thmax = Math.log(Double.MAX_VALUE / 1e6) - Math.log(maxrate);
        final double thmin = Math.log(Double.MIN_NORMAL * 1e6) - Math.log(maxrate);

        double lambda = perronstate(d0, d1, pi, ones, 0.0, nph).deta;
        SaddlepointResult res = new SaddlepointResult(n, lambda);
        for (int i = 0; i < n; i++) {
            res.p[i] = 0.0;
            res.logp[i] = Double.NEGATIVE_INFINITY;
            res.theta[i] = Double.NEGATIVE_INFINITY;
            res.eta[i] = Double.NaN;
            res.deta[i] = Double.NaN;
            res.d2eta[i] = Double.NaN;
            res.ampl[i] = Double.NaN;
            res.corr[i] = Double.NaN;
            res.k2[i] = Double.NaN;
        }

        // Sorting by the rate k/t lets each Newton solve warm-start from the
        // previous saddle, the saddle being a monotone function of that rate alone
        Integer[] ord = new Integer[n];
        final double[] rate = new double[n];
        for (int i = 0; i < n; i++) {
            ord[i] = Integer.valueOf(i);
            rate[i] = tv[i] > 0 ? kv[i] / tv[i] : 0.0;
        }
        java.util.Arrays.sort(ord, new java.util.Comparator<Integer>() {
            public int compare(Integer a, Integer b) {
                return Double.compare(rate[a.intValue()], rate[b.intValue()]);
            }
        });

        double worstK2 = Double.POSITIVE_INFINITY;
        double worstT = 0.0;
        int worstK = 0;
        double thprev = 0.0;

        for (int idx = 0; idx < n; idx++) {
            int i = ord[idx].intValue();
            double ti = tv[i];
            int ki = kv[i];
            if (ti == 0.0) {
                // No time has elapsed, so the count is zero with probability one
                res.exact[i] = true;
                if (ki == 0) {
                    res.p[i] = 1.0;
                    res.logp[i] = 0.0;
                }
                continue;
            }
            if (ki == 0) {
                // The saddle runs off to -Infinity; the exact value is one matrix
                // exponential of the taboo generator and costs no more than a
                // step of the approximation itself
                res.exact[i] = true;
                Matrix Et = D0.scale(ti).expm();
                double acc = 0.0;
                for (int a = 0; a < nph; a++) {
                    for (int b = 0; b < nph; b++) {
                        acc += pi[a] * Et.get(a, b);
                    }
                }
                res.p[i] = acc;
                res.logp[i] = acc > 0 ? Math.log(acc) : Double.NEGATIVE_INFINITY;
                continue;
            }

            double[] saddle = solvesaddle(d0, d1, pi, ones, nph, ki / ti, thprev, thmin, thmax);
            double th = saddle[0];
            thprev = th;
            res.theta[i] = th;
            res.iter[i] = (int) saddle[1];

            PerronState s = perronstate(d0, d1, pi, ones, th, nph);
            res.eta[i] = s.eta;
            res.deta[i] = s.deta;
            res.d2eta[i] = s.d2eta;

            double K2 = ti * s.d2eta;
            res.k2[i] = K2;
            if (K2 < worstK2) {
                worstK2 = K2;
                worstT = ti;
                worstK = ki;
            }
            if (!(K2 > 0.0)) {
                throw new RuntimeException("ctmc_saddlepoint: the cumulant generating function is "
                        + "not strictly convex at the saddle (t=" + ti + ", k=" + ki + "): eta''="
                        + s.d2eta + ". D0+D1 is probably reducible.");
            }
            double base = ti * s.eta - ki * th - 0.5 * Math.log(2.0 * Math.PI * K2);
            double ampl = useampl ? s.ampl : 1.0;
            res.ampl[i] = s.ampl;

            double corr;
            if (order == 1) {
                corr = ampl;
            } else {
                // The higher cumulants and the derivatives of the amplitude come
                // from central differences of the analytic eta'' and g, both of
                // which carry full precision at each evaluation point
                double h = 1e-3 * Math.max(1.0, Math.abs(th));
                PerronState sp = perronstate(d0, d1, pi, ones, th + h, nph);
                PerronState sm = perronstate(d0, d1, pi, ones, th - h, nph);
                double d3 = (sp.d2eta - sm.d2eta) / (2.0 * h);
                double d4 = (sp.d2eta - 2.0 * s.d2eta + sm.d2eta) / (h * h);
                double K3 = ti * d3;
                double K4 = ti * d4;
                double lam3sq = K3 * K3 / (K2 * K2 * K2);
                double lam4 = K4 / (K2 * K2);
                double gp = useampl ? (sp.ampl - sm.ampl) / (2.0 * h) : 0.0;
                double gpp = useampl ? (sp.ampl - 2.0 * s.ampl + sm.ampl) / (h * h) : 0.0;
                // Steepest descent to O(1/K2), Jensen (1991) Eq. (4.4): the
                // Daniels bracket on the amplitude, plus the two terms the
                // amplitude contributes through its own curvature along the contour
                corr = ampl * (1.0 + lam4 / 8.0 - 5.0 * lam3sq / 24.0)
                        - gpp / (2.0 * K2) + gp * K3 / (2.0 * K2 * K2);
                if (corr <= 0.0) {
                    InputOutput.line_warning("Ctmc_saddlepoint", "The second-order correction is "
                            + "nonpositive at t=%s, k=%d; the expansion has broken down, "
                            + "returning the first-order value.%n", Double.valueOf(ti),
                            Integer.valueOf(ki));
                    corr = ampl;
                }
            }
            res.corr[i] = corr;
            res.logp[i] = base + Math.log(corr);
            res.p[i] = Math.exp(res.logp[i]);
        }

        if (worstK2 < K2_MIN) {
            // Once per call, not once per point: a vectorised call spans hundreds
            // of counts and the caller needs the worst one, not a page of repetitions
            InputOutput.line_warning("Ctmc_saddlepoint", "K2 = t*eta''(theta*) = %.2f at t=%s, "
                    + "k=%d is below %s, so the saddlepoint expansion is outside its asymptotic "
                    + "regime there and the result is unreliable (expect a relative error near "
                    + "%.0e). K2 is the variance of the count, not its mean: a low-variability "
                    + "process needs a longer horizon than its rate suggests. Take the exact "
                    + "value from the block chain instead.%n",
                    Double.valueOf(worstK2), Double.valueOf(worstT), Integer.valueOf(worstK),
                    Double.valueOf(K2_MIN), Double.valueOf(0.017 / (worstK2 * worstK2)));
        }
        return res;
    }

    private static double[][] toArray(Matrix M, int nph) {
        double[][] a = new double[nph][nph];
        for (int i = 0; i < nph; i++) {
            for (int j = 0; j < nph; j++) {
                a[i][j] = M.get(i, j);
            }
        }
        return a;
    }

    /**
     * Perron root of A(th) = D0 + exp(th)*D1 with deta, d2eta and the amplitude.
     *
     * <p>Only EIGENVALUES are taken from the eigensolver; the Perron vectors come
     * from bordered solves, the idiom Ctmc_solve already uses. That keeps all
     * four codebases on one algorithm: commons-math hands back Schur blocks
     * rather than eigenvectors as soon as a complex pair appears, and the C++
     * eig.h exposes values only, so neither can supply a left eigenvector.
     */
    static PerronState perronstate(double[][] d0, double[][] d1, double[] pi0, double[] ones,
                                   double th, int nph) {
        double ex = Math.exp(th);
        double[][] W = new double[nph][nph];        // A'(th) = A''(th) = exp(th)*D1
        double[][] A = new double[nph][nph];
        for (int i = 0; i < nph; i++) {
            for (int j = 0; j < nph; j++) {
                W[i][j] = ex * d1[i][j];
                A[i][j] = d0[i][j] + W[i][j];
            }
        }
        Matrix Am = new Matrix(nph, nph);
        for (int i = 0; i < nph; i++) {
            for (int j = 0; j < nph; j++) {
                Am.set(i, j, A[i][j]);
            }
        }
        List<Complex> ev = Am.eig();
        double eta = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < ev.size(); i++) {
            double re = ev.get(i).getReal();
            if (re > eta) {
                eta = re;
            }
        }

        double[][] Ashift = new double[nph][nph];
        for (int i = 0; i < nph; i++) {
            for (int j = 0; j < nph; j++) {
                Ashift[i][j] = A[i][j] - (i == j ? eta : 0.0);
            }
        }
        double[] rhs = new double[nph];
        rhs[nph - 1] = 1.0;
        // (A-eta*I)v = 0 with the last row replaced by sum(v)=1. A row may be
        // dropped because A-eta*I is a singular irreducible M-matrix, every
        // proper principal submatrix of which is nonsingular
        double[][] M = new double[nph][nph];
        for (int i = 0; i < nph; i++) {
            for (int j = 0; j < nph; j++) {
                M[i][j] = (i == nph - 1) ? 1.0 : Ashift[i][j];
            }
        }
        double[] v = gauss(M, rhs, nph);
        // u(A-eta*I) = 0 by the same construction on the transpose
        double[][] Mt = new double[nph][nph];
        for (int i = 0; i < nph; i++) {
            for (int j = 0; j < nph; j++) {
                Mt[i][j] = (i == nph - 1) ? 1.0 : Ashift[j][i];
            }
        }
        double[] u = gauss(Mt, rhs, nph);
        double uv = 0.0;
        for (int i = 0; i < nph; i++) {
            uv += u[i] * v[i];
        }
        for (int i = 0; i < nph; i++) {
            u[i] /= uv;                              // u*v = 1 fixes the residual scale
        }
        double deta = quad(u, W, v, nph);

        // First-order eigenvector perturbation (A-eta*I)v' = (eta'*I-W)v taken
        // with u*v'=0; the bordered system is nonsingular because the Perron
        // root of an irreducible Metzler matrix is simple
        double[][] B = new double[nph + 1][nph + 1];
        double[] r = new double[nph + 1];
        for (int i = 0; i < nph; i++) {
            for (int j = 0; j < nph; j++) {
                B[i][j] = Ashift[i][j];
            }
            B[i][nph] = v[i];
            B[nph][i] = u[i];
            double acc = deta * v[i];
            for (int j = 0; j < nph; j++) {
                acc -= W[i][j] * v[j];
            }
            r[i] = acc;
        }
        double[] sol = gauss(B, r, nph + 1);
        double[] vp = new double[nph];
        for (int i = 0; i < nph; i++) {
            vp[i] = sol[i];
        }
        double d2eta = deta + 2.0 * quad(u, W, vp, nph);

        double pv = 0.0;
        double u1 = 0.0;
        for (int i = 0; i < nph; i++) {
            pv += pi0[i] * v[i];
            u1 += u[i] * ones[i];
        }
        return new PerronState(eta, deta, d2eta, pv * u1);
    }

    private static double quad(double[] u, double[][] W, double[] v, int nph) {
        double acc = 0.0;
        for (int i = 0; i < nph; i++) {
            double inner = 0.0;
            for (int j = 0; j < nph; j++) {
                inner += W[i][j] * v[j];
            }
            acc += u[i] * inner;
        }
        return acc;
    }

    /**
     * Dense solve with partial pivoting. The systems here are of order K or K+1
     * with K the number of phases, so a self-contained elimination is cheaper
     * than routing small dense problems through the Matrix backends, and it
     * keeps this file free of storage-format concerns.
     */
    private static double[] gauss(double[][] Ain, double[] bin, int n) {
        double[][] a = new double[n][n + 1];
        for (int i = 0; i < n; i++) {
            System.arraycopy(Ain[i], 0, a[i], 0, n);
            a[i][n] = bin[i];
        }
        for (int c = 0; c < n; c++) {
            int piv = c;
            double best = Math.abs(a[c][c]);
            for (int i = c + 1; i < n; i++) {
                if (Math.abs(a[i][c]) > best) {
                    best = Math.abs(a[i][c]);
                    piv = i;
                }
            }
            if (best == 0.0) {
                throw new RuntimeException("ctmc_saddlepoint: singular system at the Perron pair; "
                        + "D0+D1 is probably reducible.");
            }
            if (piv != c) {
                double[] tmp = a[piv];
                a[piv] = a[c];
                a[c] = tmp;
            }
            for (int i = c + 1; i < n; i++) {
                double f = a[i][c] / a[c][c];
                if (f == 0.0) {
                    continue;
                }
                for (int j = c; j <= n; j++) {
                    a[i][j] -= f * a[c][j];
                }
            }
        }
        double[] x = new double[n];
        for (int i = n - 1; i >= 0; i--) {
            double acc = a[i][n];
            for (int j = i + 1; j < n; j++) {
                acc -= a[i][j] * x[j];
            }
            x[i] = acc / a[i][i];
        }
        return x;
    }

    /**
     * Saddle of the counting cumulant generating function at rate r, the root of
     * eta'(th) = r. eta' is continuous and strictly increasing from 0 to
     * +Infinity, so the root exists and is unique for every r&gt;0; it is
     * bracketed by geometric expansion from th0 and refined by Newton on
     * log(eta'), safeguarded by bisection.
     *
     * @return two entries, the saddle and the Newton step count
     */
    static double[] solvesaddle(double[][] d0, double[][] d1, double[] pi0, double[] ones, int nph,
                                double r, double th0, double thmin, double thmax) {
        final double TOL = 1e-13;
        final int MAXIT = 200;
        double th = Math.min(Math.max(th0, thmin), thmax);
        double d1v = perronstate(d0, d1, pi0, ones, th, nph).deta;
        double lo = th;
        double hi = th;
        double dlo = d1v;
        double dhi = d1v;
        double step = 1.0;
        while (dlo > r) {
            hi = lo;
            dhi = dlo;
            lo = lo - step;
            if (lo <= thmin) {
                lo = thmin;
                dlo = perronstate(d0, d1, pi0, ones, lo, nph).deta;
                if (dlo > r) {
                    throw new RuntimeException("ctmc_saddlepoint: the rate k/t=" + r
                            + " is below the representable range of eta'.");
                }
                break;
            }
            dlo = perronstate(d0, d1, pi0, ones, lo, nph).deta;
            step *= 2.0;
        }
        step = 1.0;
        while (dhi < r) {
            lo = hi;
            dlo = dhi;
            hi = hi + step;
            if (hi >= thmax) {
                hi = thmax;
                dhi = perronstate(d0, d1, pi0, ones, hi, nph).deta;
                if (dhi < r) {
                    throw new RuntimeException("ctmc_saddlepoint: the rate k/t=" + r
                            + " is above the representable range of eta'.");
                }
                break;
            }
            dhi = perronstate(d0, d1, pi0, ones, hi, nph).deta;
            step *= 2.0;
        }
        th = Math.min(Math.max(th, lo), hi);
        double logr = Math.log(r);
        int iters = 0;
        for (int it = 1; it <= MAXIT; it++) {
            iters = it;
            PerronState si = perronstate(d0, d1, pi0, ones, th, nph);
            double f = Math.log(si.deta) - logr;
            if (Math.abs(f) <= TOL) {
                break;
            }
            if (f > 0.0) {
                hi = th;
            } else {
                lo = th;
            }
            double thn = th - f * si.deta / si.d2eta;
            if (Double.isNaN(thn) || Double.isInfinite(thn) || thn <= lo || thn >= hi) {
                thn = 0.5 * (lo + hi);
            }
            if (Math.abs(thn - th) <= TOL * Math.max(1.0, Math.abs(th))) {
                th = thn;
                break;
            }
            th = thn;
        }
        return new double[]{th, iters};
    }
}
