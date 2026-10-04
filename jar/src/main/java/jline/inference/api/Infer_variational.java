/*
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
package jline.inference.api;

import org.apache.commons.math3.special.Gamma;

/**
 * Variational inference for Markovian queueing networks, following I. Perez,
 * G. Casale, "Variational Inference for Markovian Queueing Networks", Advances
 * in Applied Probability 53(3), 2021.
 *
 * <p>The network trajectory is reparameterised by the transition counts
 * Y^eta, eta=(i,j,c), so that the station marginals decouple:
 * x_{i,c}(t) = x_{i,c}(0) + sum_{In(i,c)} Y^eta - sum_{Out(i,c)} Y^eta. The
 * variational family is a product of inhomogeneous pure-birth processes, one
 * per transition, with rate nu^eta(t,y), times a product of Gamma densities
 * over the unknown service rates. The state space is expanded by adding DELTA
 * to every feasible rate, so that queue lengths may go negative and the
 * approximating measure stays mutually absolutely continuous with the target;
 * the original model is recovered as DELTA to 0.</p>
 *
 * <p>Each iteration performs, per transition, a backward pass for the Lagrange
 * multipliers r^eta with multiplicative jumps at the observation epochs, the
 * rate update nu^eta(t,y) = exp(E log Xi^eta(t,y)) r^eta(t,y+1)/r^eta(t,y),
 * and a forward pass of the master equation for the marginal. The conjugate
 * Gamma posteriors are then refreshed from the expected number of firings and
 * the expected exposure time of each station-class pair.</p>
 *
 * <p>Expectations over the other transitions are taken on a deterministic
 * Halton lattice mapped through the inverse marginal c.d.f., so the estimator
 * carries no random-number stream and reproduces the MATLAB, Python and C++
 * implementations digit for digit.</p>
 */
public class Infer_variational {

    private Infer_variational() {}

    /** Load factor Upsilon of a transition leaving a station-class pair. */
    static double ups(double xic, double xis, double nservers, int sched,
                      double cap, double capstat) {
        if (sched == 2) {
            return 1.0;
        }
        double a = Math.min(cap, Math.max(0.0, xic));
        if (sched == 0) {
            return a;
        }
        double b = Math.min(capstat, Math.max(0.0, xis));
        if (b <= 0) {
            return 0.0;
        }
        return a / b * Math.min(nservers, b);
    }

    /** k-th prime, k >= 1. */
    static int prime(int k) {
        int n = 0;
        int c = 1;
        int p = 2;
        while (n < k) {
            c++;
            boolean isp = true;
            for (int d = 2; d * d <= c; d++) {
                if (c % d == 0) {
                    isp = false;
                    break;
                }
            }
            if (isp) {
                n++;
                p = c;
            }
        }
        return p;
    }

    /** Van der Corput radical inverse of i in the given base. */
    static double radicalInverse(int i, int base) {
        double r = 0.0;
        double f = 1.0 / base;
        while (i > 0) {
            r += f * (i % base);
            i /= base;
            f /= base;
        }
        return r;
    }

    /**
     * Inverse-c.d.f. samples of a marginal on a Halton lattice. Each transition
     * uses its own prime base, so the samples of distinct transitions are
     * jointly equidistributed rather than comonotone.
     */
    static double[][] sample(double[][] q, int S, int e) {
        int G = q.length;
        int ny = q[0].length;
        int base = prime(e + 1);
        double[] u = new double[S];
        Integer[] ord = new Integer[S];
        for (int s = 0; s < S; s++) {
            u[s] = radicalInverse(s + 1, base);
            ord[s] = s;
        }
        final double[] uu = u;
        java.util.Arrays.sort(ord, new java.util.Comparator<Integer>() {
            public int compare(Integer a, Integer b) {
                return Double.compare(uu[a], uu[b]);
            }
        });
        double[] us = new double[S];
        for (int s = 0; s < S; s++) {
            us[s] = u[ord[s]];
        }
        double[][] ys = new double[G][S];
        double[] c = new double[ny];
        for (int g = 0; g < G; g++) {
            double acc = 0.0;
            for (int y = 0; y < ny; y++) {
                acc += q[g][y];
                c[y] = acc;
            }
            if (c[ny - 1] > 0) {
                for (int y = 0; y < ny; y++) {
                    c[y] /= c[ny - 1];
                }
            }
            c[ny - 1] = 1.0;
            int j = 0;
            for (int s = 0; s < S; s++) {
                while (j < ny - 1 && c[j] < us[s]) {
                    j++;
                }
                ys[g][ord[s]] = j;
            }
        }
        return ys;
    }

    static double[][][] sampleAll(double[][][] Y, int S) {
        int narcs = Y.length;
        double[][][] ys = new double[narcs][][];
        for (int e = 0; e < narcs; e++) {
            ys[e] = sample(Y[e], S, e);
        }
        return ys;
    }

    /** Slack multiplier of the rate cap; unity when the cap is inactive. */
    static double damp(double sl, double ye, double floor) {
        if (sl == 0) {
            return 1.0;
        }
        double z = sl / Math.max(floor, ye);
        double d = (1.0 + z) / Math.exp(z);
        if (!(d >= 0) || Double.isInfinite(d)) {
            return 0.0;
        }
        return d;
    }

    static void rescale(double[] v) {
        double m = 0.0;
        for (int i = 0; i < v.length; i++) {
            if (v[i] > m) {
                m = v[i];
            }
        }
        if (m > 0 && !Double.isInfinite(m)) {
            for (int i = 0; i < v.length; i++) {
                v[i] /= m;
            }
        }
    }

    /** One uniformization step of the backward sub-generator. */
    static double[] backUniformize(double[] v0, double[] pd, double[] pu, double lt,
                                   VariationalOptions opt) {
        int ny = v0.length;
        double w = Math.exp(-lt);
        double[] v = new double[ny];
        double[] u = new double[ny];
        for (int i = 0; i < ny; i++) {
            v[i] = w * v0[i];
            u[i] = v0[i];
        }
        double cum = w;
        int n = 1;
        double[] un = new double[ny];
        while ((1.0 - cum) > opt.unifTol && n < opt.unifMaxTerms) {
            for (int i = 0; i < ny; i++) {
                un[i] = u[i] * (1.0 - pd[i]);
            }
            for (int i = 0; i < ny - 1; i++) {
                un[i] += u[i + 1] * pu[i];
            }
            System.arraycopy(un, 0, u, 0, ny);
            w = w * lt / n;
            for (int i = 0; i < ny; i++) {
                v[i] += w * u[i];
            }
            cum += w;
            n++;
        }
        return v;
    }

    /**
     * Backward pass for the Lagrange multipliers. The equation is linear in r,
     * so on a grid cell with frozen coefficients it is the action of a matrix
     * exponential. The generator has non-positive row sums by Jensen, so
     * uniformization evaluates it without the stiffness that an explicit rule
     * suffers when exp(E log Xi) falls orders of magnitude below E[Xi]. Only
     * the ratios r(y+1)/r(y) are used downstream, so r is rescaled at each step.
     */
    static double[][] backward(double[][] ge, double[][] he, double[][] sl, double[][] Ye,
                               int[] obsIdx, double[][] obsw, double dt,
                               VariationalOptions opt) {
        int G = ge.length;
        int ny = ge[0].length;
        double[][] r = new double[G][ny];
        double[] v = new double[ny];
        for (int y = 0; y < ny; y++) {
            v[y] = 1.0;
        }
        for (int q = 0; q < obsIdx.length; q++) {
            if (obsIdx[q] == G - 1) {
                for (int y = 0; y < ny; y++) {
                    v[y] *= Math.max(0.0, obsw[q][y]);
                }
            }
        }
        rescale(v);
        System.arraycopy(v, 0, r[G - 1], 0, ny);
        double[] pd = new double[ny];
        double[] pu = new double[ny];
        for (int g = G - 2; g >= 0; g--) {
            double lam = 0.0;
            for (int y = 0; y < ny; y++) {
                double gv = Math.max(0.0, ge[g][y]);
                double hv = Math.max(0.0, he[g][y] * damp(sl[g][y], Ye[g][y], opt.floor));
                pd[y] = gv;
                pu[y] = Math.min(hv, gv);
                if (gv > lam) {
                    lam = gv;
                }
            }
            if (lam > 0) {
                int ncell = Math.max(1, (int) Math.ceil(lam * dt / opt.unifmax));
                double h = dt / ncell;
                double[] pdn = new double[ny];
                double[] pun = new double[ny];
                for (int y = 0; y < ny; y++) {
                    pdn[y] = pd[y] / lam;
                    pun[y] = pu[y] / lam;
                }
                for (int c = 0; c < ncell; c++) {
                    v = backUniformize(v, pdn, pun, lam * h, opt);
                }
                rescale(v);
            }
            for (int q = 0; q < obsIdx.length; q++) {
                if (obsIdx[q] == g) {
                    for (int y = 0; y < ny; y++) {
                        v[y] *= Math.max(0.0, obsw[q][y]);
                    }
                    rescale(v);
                }
            }
            System.arraycopy(v, 0, r[g], 0, ny);
        }
        return r;
    }

    /** One uniformization step of the pure-birth chain. */
    static double[] uniformize(double[] v0, double[] p, double lt, VariationalOptions opt) {
        int ny = v0.length;
        double w = Math.exp(-lt);
        double[] v = new double[ny];
        double[] u = new double[ny];
        for (int i = 0; i < ny; i++) {
            v[i] = w * v0[i];
            u[i] = v0[i];
        }
        double cum = w;
        int n = 1;
        double[] un = new double[ny];
        while ((1.0 - cum) > opt.unifTol && n < opt.unifMaxTerms) {
            for (int i = 0; i < ny; i++) {
                un[i] = u[i] * (1.0 - p[i]);
            }
            for (int i = ny - 1; i >= 1; i--) {
                un[i] += u[i - 1] * p[i - 1];
            }
            System.arraycopy(un, 0, u, 0, ny);
            w = w * lt / n;
            for (int i = 0; i < ny; i++) {
                v[i] += w * u[i];
            }
            cum += w;
            n++;
        }
        return v;
    }

    /** Forward master equation of an inhomogeneous pure-birth process. */
    static double[][] forward(double[][] nue, double dt, VariationalOptions opt) {
        int G = nue.length;
        int ny = nue[0].length;
        double[][] q = new double[G][ny];
        double[] v = new double[ny];
        v[0] = 1.0;
        System.arraycopy(v, 0, q[0], 0, ny);
        double[] p = new double[ny];
        for (int g = 0; g < G - 1; g++) {
            double lam = 0.0;
            for (int y = 0; y < ny; y++) {
                double rate = Math.max(0.0, nue[g][y]);
                p[y] = rate;
                if (rate > lam) {
                    lam = rate;
                }
            }
            if (lam <= 0) {
                System.arraycopy(v, 0, q[g + 1], 0, ny);
                continue;
            }
            int ncell = Math.max(1, (int) Math.ceil(lam * dt / opt.unifmax));
            double h = dt / ncell;
            double[] pn = new double[ny];
            for (int y = 0; y < ny; y++) {
                pn[y] = p[y] / lam;
            }
            for (int c = 0; c < ncell; c++) {
                v = uniformize(v, pn, lam * h, opt);
            }
            double s = 0.0;
            for (int y = 0; y < ny; y++) {
                if (v[y] < 0) {
                    v[y] = 0;
                }
                s += v[y];
            }
            if (s > 0) {
                for (int y = 0; y < ny; y++) {
                    v[y] /= s;
                }
            }
            System.arraycopy(v, 0, q[g + 1], 0, ny);
        }
        return q;
    }

    /** Trapezoidal integral of a grid function. */
    static double trapz(double[] f, double dt) {
        int n = f.length;
        if (n < 2) {
            return 0.0;
        }
        double s = 0.0;
        for (int i = 0; i < n; i++) {
            s += f[i];
        }
        return dt * (s - 0.5 * f[0] - 0.5 * f[n - 1]);
    }

    /** KL(Gamma(a,b) || Gamma(a0,b0)) with rate parameterisation. */
    static double klGamma(double a, double b, double a0, double b0) {
        return (a - a0) * Gamma.digamma(a) - Gamma.logGamma(a) + Gamma.logGamma(a0)
                + a0 * (Math.log(b) - Math.log(b0)) + a * (b0 - b) / b;
    }

    /** Run the variational inference procedure. */
    public static VariationalResult infer_variational(VariationalSpec spec,
                                                      VariationalOptions options) {
        VariationalOptions opt = options == null ? new VariationalOptions() : options;
        int M = spec.nstations();
        int R = spec.nclasses();
        int narcs = spec.narcs();
        int P = spec.nparams();

        validate(spec, M, R, narcs);
        if (spec.capacity == null) {
            spec.capacity = new double[M * R];
            for (int k = 0; k < M * R; k++) {
                spec.capacity[k] = Double.POSITIVE_INFINITY;
            }
        }

        int[] src = new int[narcs];
        int[] dst = new int[narcs];
        int[] cls = new int[narcs];
        for (int e = 0; e < narcs; e++) {
            src[e] = spec.arcs[e][0];
            dst[e] = spec.arcs[e][1];
            cls[e] = spec.arcs[e][2];
        }

        double[][] sgnClass = new double[narcs][M * R];
        double[][] sgnStat = new double[narcs][M];
        for (int e = 0; e < narcs; e++) {
            if (dst[e] > 0) {
                sgnClass[e][(cls[e] - 1) * M + dst[e] - 1] += 1.0;
                sgnStat[e][dst[e] - 1] += 1.0;
            }
            if (src[e] > 0) {
                sgnClass[e][(cls[e] - 1) * M + src[e] - 1] -= 1.0;
                sgnStat[e][src[e] - 1] -= 1.0;
            }
        }

        double[] x0v = VariationalSpec.flatten(spec.x0);
        double[] x0s = new double[M];
        for (int m = 0; m < M; m++) {
            for (int r = 0; r < R; r++) {
                x0s[m] += spec.x0[m][r];
            }
        }
        double[] capStat = new double[M];
        for (int m = 0; m < M; m++) {
            for (int r = 0; r < R; r++) {
                capStat[m] += spec.capacity[r * M + m];
            }
        }

        // mean occupancy used to size the truncation and the initial rates
        double[] xbar = new double[M * R];
        System.arraycopy(x0v, 0, xbar, 0, M * R);
        for (int k = 0; k < M * R; k++) {
            double s = 0.0;
            int n = 0;
            for (int q = 0; q < spec.obsTimes.length; q++) {
                if (!Double.isNaN(spec.obsData[q][k])) {
                    s += spec.obsData[q][k];
                    n++;
                }
            }
            if (n > 0) {
                xbar[k] = s / n;
            }
        }
        double[] xbars = new double[M];
        for (int m = 0; m < M; m++) {
            for (int r = 0; r < R; r++) {
                xbars[m] += xbar[r * M + m];
            }
        }

        if (opt.tmax == null) {
            if (spec.obsTimes.length == 0) {
                throw new IllegalArgumentException(
                        "options.tmax is required when there are no observations.");
            }
            double t = 0.0;
            for (int k = 0; k < spec.obsTimes.length; k++) {
                t = Math.max(t, spec.obsTimes[k]);
            }
            opt.tmax = t;
        }
        if (opt.tmax <= 0) {
            throw new IllegalArgumentException("options.tmax must be positive.");
        }
        if (opt.ngrid == null && opt.dt == null) {
            opt.ngrid = 201;
        }
        if (opt.ngrid == null) {
            opt.ngrid = (int) Math.round(opt.tmax / opt.dt) + 1;
        }
        opt.ngrid = Math.max(2, opt.ngrid);
        opt.dt = opt.tmax / (opt.ngrid - 1);

        if (opt.ymax == null) {
            double fmax = 0.0;
            for (int e = 0; e < narcs; e++) {
                double lam;
                if (spec.arcparam[e] > 0) {
                    int p = spec.arcparam[e] - 1;
                    lam = spec.routeprob[e] * spec.alpha0[p] / spec.beta0[p];
                } else {
                    lam = spec.routeprob[e] * spec.arcrate[e];
                }
                double u = 1.0;
                int i = src[e];
                if (i > 0) {
                    int kc = (cls[e] - 1) * M + i - 1;
                    u = ups(xbar[kc], xbars[i - 1], spec.nservers[i - 1], spec.sched[i - 1],
                            spec.capacity[kc], capStat[i - 1]);
                }
                fmax = Math.max(fmax, lam * u * opt.tmax);
            }
            opt.ymax = Math.max(20, (int) Math.ceil(2 * fmax + 5 * Math.sqrt(Math.max(1.0, fmax))));
        }
        opt.ymax = Math.max(2, opt.ymax);
        if (opt.rateMax == null) {
            opt.rateMax = opt.rateCapFactor * opt.ymax / opt.tmax;
        }

        int G = opt.ngrid;
        double dt = opt.dt;
        int ymax = opt.ymax;
        int ny = ymax + 1;
        int S = opt.nsamples;
        int K = spec.obsTimes.length;
        double[] yvec = new double[ny];
        for (int y = 0; y < ny; y++) {
            yvec[y] = y;
        }
        double[] tgrid = new double[G];
        for (int g = 0; g < G; g++) {
            tgrid[g] = g * dt;
        }
        int[] obsIdx = new int[K];
        for (int k = 0; k < K; k++) {
            obsIdx[k] = Math.min(G - 1, Math.max(0, (int) Math.round(spec.obsTimes[k] / dt)));
        }

        int[] arcSched = new int[narcs];
        double[] arcServers = new double[narcs];
        for (int e = 0; e < narcs; e++) {
            if (src[e] > 0) {
                arcSched[e] = spec.sched[src[e] - 1];
                arcServers[e] = spec.nservers[src[e] - 1];
            } else {
                arcSched[e] = 2;
                arcServers[e] = 1.0;
            }
        }

        double[] alpha = new double[P];
        double[] beta = new double[P];
        System.arraycopy(spec.alpha0, 0, alpha, 0, P);
        System.arraycopy(spec.beta0, 0, beta, 0, P);

        double[][][] Y = new double[narcs][G][ny];
        double[][][] nu = new double[narcs][G][ny];
        double[][][] slack = new double[narcs][G][ny];
        double[][][] gexp = new double[narcs][G][ny];
        double[][][] hexp = new double[narcs][G][ny];

        for (int e = 0; e < narcs; e++) {
            double lam = rateMean(spec, alpha, beta, e);
            double u0 = 1.0;
            if (src[e] > 0) {
                int kc = (cls[e] - 1) * M + src[e] - 1;
                u0 = ups(xbar[kc], xbars[src[e] - 1], arcServers[e], arcSched[e],
                        spec.capacity[kc], capStat[src[e] - 1]);
            }
            double nu0 = Math.max(opt.delta, lam * u0);
            double[][] nue = new double[G][ny];
            for (int g = 0; g < G; g++) {
                for (int y = 0; y < ymax; y++) {
                    nue[g][y] = nu0;
                }
            }
            nu[e] = nue;
            Y[e] = forward(nue, dt, opt);
        }

        double[] bound = new double[opt.iterMax];
        double[][] alphaTrace = new double[P][opt.iterMax];
        double[][] betaTrace = new double[P][opt.iterMax];
        boolean converged = false;
        int iter = 0;

        double[][][] Ys = null;
        for (int it = 1; it <= opt.iterMax; it++) {
            iter = it;
            for (int e = 0; e < narcs; e++) {
                Ys = sampleAll(Y, S);
                double[][][] gh = rateMoments(spec, opt, e, Ys, Y, sgnClass, sgnStat,
                        x0v, x0s, capStat, arcSched, arcServers, alpha, beta, obsIdx, yvec, true);
                double[][] ge = gh[0];
                double[][] he = gh[1];
                double[][] obsw = gh[2];
                double[][] Ye = Y[e];
                double[][] r = backward(ge, he, slack[e], Ye, obsIdx, obsw, dt, opt);

                // Eq. (15). A vanishing multiplier marks a count the future
                // observations rule out; the rate there is zero, which is what
                // keeps the forward pass from placing mass on it.
                double[][] nue = new double[G][ny];
                double[][] sl = new double[G][ny];
                for (int g = 0; g < G; g++) {
                    for (int y = 0; y < ymax; y++) {
                        double den = r[g][y];
                        double val = den > 0 ? he[g][y] * r[g][y + 1] / den : 0.0;
                        if (!(val > 0) || Double.isInfinite(val)) {
                            val = 0.0;
                        }
                        if (val > opt.rateMax) {
                            sl[g][y] = Math.max(opt.floor, Ye[g][y]) * Math.log(val / opt.rateMax);
                            val = opt.rateMax;
                        }
                        nue[g][y] = val;
                    }
                }
                nu[e] = nue;
                slack[e] = sl;
                Y[e] = forward(nue, dt, opt);
            }

            // conjugate Gamma updates: the shape gains the expected number of
            // firings, the rate the expected exposure time of the station-class
            // pair that the parameter governs
            Ys = sampleAll(Y, S);
            double[] firings = new double[P];
            double[] exposure = new double[P];
            boolean[][] seen = new boolean[P][M * R];
            for (int e = 0; e < narcs; e++) {
                int p = spec.arcparam[e];
                if (p == 0) {
                    continue;
                }
                // expected number of firings over the horizon, taken from the
                // marginal itself, which is exact, rather than by quadrature of
                // the intensity, which a near-deterministic marginal makes
                // inaccurate
                double m1 = 0.0;
                double m0 = 0.0;
                for (int y = 0; y < ny; y++) {
                    m1 += Y[e][G - 1][y] * y;
                    m0 += Y[e][0][y] * y;
                }
                firings[p - 1] += m1 - m0;
                int kclass = (cls[e] - 1) * M + src[e] - 1;
                if (!seen[p - 1][kclass]) {
                    seen[p - 1][kclass] = true;
                    double[] ue = new double[G];
                    for (int g = 0; g < G; g++) {
                        double acc = 0.0;
                        for (int s = 0; s < S; s++) {
                            double a = x0v[kclass];
                            double b = x0s[src[e] - 1];
                            for (int f = 0; f < narcs; f++) {
                                a += sgnClass[f][kclass] * Ys[f][g][s];
                                b += sgnStat[f][src[e] - 1] * Ys[f][g][s];
                            }
                            acc += ups(a, b, spec.nservers[src[e] - 1], spec.sched[src[e] - 1],
                                    spec.capacity[kclass], capStat[src[e] - 1]);
                        }
                        ue[g] = acc / S;
                    }
                    exposure[p - 1] += trapz(ue, dt);
                }
            }
            for (int p = 0; p < P; p++) {
                alpha[p] = spec.alpha0[p] + firings[p];
                beta[p] = spec.beta0[p] + exposure[p];
                alphaTrace[p][it - 1] = alpha[p];
                betaTrace[p][it - 1] = beta[p];
            }

            // the bound is evaluated at the state the iteration ended in, so
            // the rate moments are recomputed against the updated marginals
            // rather than reused from the sweep that produced them
            for (int e = 0; e < narcs; e++) {
                double[][][] gh = rateMoments(spec, opt, e, Ys, Y, sgnClass, sgnStat,
                        x0v, x0s, capStat, arcSched, arcServers, alpha, beta, obsIdx, yvec, false);
                gexp[e] = gh[0];
                hexp[e] = gh[1];
            }

            bound[it - 1] = bound(spec, opt, alpha, beta, Y, nu, gexp, hexp, Ys,
                    sgnClass, x0v, obsIdx, dt);
            if (opt.verbose > 0) {
                double maxnu = 0.0;
                for (int e = 0; e < narcs; e++) {
                    for (int g = 0; g < G; g++) {
                        for (int y = 0; y < ny; y++) {
                            maxnu = Math.max(maxnu, nu[e][g][y]);
                        }
                    }
                }
                System.out.printf("infer_variational: iteration %d, lower bound %.6f, max rate %.3f%n",
                        it, bound[it - 1], maxnu);
            }
            // The rate update solves a stationarity condition rather than
            // maximising the bound in a block, so the bound need not ascend;
            // convergence is judged on the bound AND on the rate posteriors.
            if (it > 1) {
                double crit = Math.abs(bound[it - 1] - bound[it - 2])
                        / Math.max(1.0, Math.abs(bound[it - 2]));
                for (int p = 0; p < P; p++) {
                    double prev = alphaTrace[p][it - 2] / betaTrace[p][it - 2];
                    crit = Math.max(crit, Math.abs(alpha[p] / beta[p] - prev)
                            / Math.max(1e-12, prev));
                }
                if (crit <= opt.tol) {
                    converged = true;
                    break;
                }
            }
        }

        double tailmass = 0.0;
        for (int e = 0; e < narcs; e++) {
            for (int g = 0; g < G; g++) {
                tailmass = Math.max(tailmass, Y[e][g][ny - 1]);
            }
        }
        if (tailmass > 1e-6) {
            System.err.printf("Warning [Infer_variational]: transition-count truncation ymax=%d "
                    + "carries mass %.3e, increase options.ymax.%n", ymax, tailmass);
        }

        double[][] qlen = new double[G][M * R];
        for (int g = 0; g < G; g++) {
            System.arraycopy(x0v, 0, qlen[g], 0, M * R);
        }
        for (int e = 0; e < narcs; e++) {
            for (int g = 0; g < G; g++) {
                double my = 0.0;
                for (int y = 0; y < ny; y++) {
                    my += Y[e][g][y] * y;
                }
                for (int k = 0; k < M * R; k++) {
                    qlen[g][k] += my * sgnClass[e][k];
                }
            }
        }

        VariationalResult out = new VariationalResult();
        out.alpha = alpha;
        out.beta = beta;
        out.rates = new double[P];
        out.meanServiceTime = new double[P];
        for (int p = 0; p < P; p++) {
            out.rates[p] = alpha[p] / beta[p];
            out.meanServiceTime[p] = beta[p] / alpha[p];
        }
        out.bound = java.util.Arrays.copyOf(bound, iter);
        out.alphaTrace = new double[P][iter];
        out.betaTrace = new double[P][iter];
        for (int p = 0; p < P; p++) {
            System.arraycopy(alphaTrace[p], 0, out.alphaTrace[p], 0, iter);
            System.arraycopy(betaTrace[p], 0, out.betaTrace[p], 0, iter);
        }
        out.Y = Y;
        out.nu = nu;
        out.tgrid = tgrid;
        out.qlen = qlen;
        out.iter = iter;
        out.converged = converged;
        out.tailmass = tailmass;
        return out;
    }

    private static void validate(VariationalSpec spec, int M, int R, int narcs) {
        for (int e = 0; e < narcs; e++) {
            if (spec.arcs[e][0] == 0 && spec.arcs[e][1] == 0) {
                throw new IllegalArgumentException("A transition cannot be external at both ends.");
            }
            if (spec.arcparam[e] == 0 && !(spec.arcrate[e] > 0)) {
                throw new IllegalArgumentException("Transition " + e
                        + " has no parameter and no positive known rate.");
            }
        }
        if (spec.sched.length != M || spec.nservers.length != M) {
            throw new IllegalArgumentException(
                    "sched and nservers must have one entry per station.");
        }
        if (spec.routeprob.length != narcs || spec.arcparam.length != narcs
                || spec.arcrate.length != narcs) {
            throw new IllegalArgumentException(
                    "routeprob, arcparam and arcrate must have one entry per transition.");
        }
        if (spec.alpha0.length != spec.beta0.length) {
            throw new IllegalArgumentException("alpha0 and beta0 must have the same length.");
        }
        if (spec.obsData.length != spec.obsTimes.length) {
            throw new IllegalArgumentException(
                    "obsData must have one row per observation epoch.");
        }
        for (int k = 0; k < spec.obsData.length; k++) {
            if (spec.obsData[k].length != M * R) {
                throw new IllegalArgumentException("obsData must have M*R columns.");
            }
        }
        if (spec.obsRange.length != M * R) {
            throw new IllegalArgumentException("obsRange must have M*R entries.");
        }
        if (spec.capacity != null && spec.capacity.length != M * R) {
            throw new IllegalArgumentException("capacity must have M*R entries.");
        }
    }

    static double rateMean(VariationalSpec spec, double[] alpha, double[] beta, int e) {
        int p = spec.arcparam[e];
        if (p == 0) {
            return spec.routeprob[e] * spec.arcrate[e];
        }
        return spec.routeprob[e] * alpha[p - 1] / beta[p - 1];
    }

    static double rateLogMean(VariationalSpec spec, double[] alpha, double[] beta, int e) {
        int p = spec.arcparam[e];
        if (p == 0) {
            return Math.log(spec.routeprob[e] * spec.arcrate[e]);
        }
        return Math.log(spec.routeprob[e]) + Gamma.digamma(alpha[p - 1]) - Math.log(beta[p - 1]);
    }

    /**
     * Conditional rate moments of one transition, and its observation jumps.
     * Returns {E[Xi|Y^eta=y], exp(E[log Xi|Y^eta=y]), observation weights},
     * the first two taken under Q with the transition's own contribution
     * removed.
     */
    static double[][][] rateMoments(VariationalSpec spec, VariationalOptions opt, int e,
                                    double[][][] Ys, double[][][] Y, double[][] sgnClass,
                                    double[][] sgnStat, double[] x0v, double[] x0s,
                                    double[] capStat, int[] arcSched, double[] arcServers,
                                    double[] alpha, double[] beta, int[] obsIdx, double[] yvec,
                                    boolean wantObs) {
        int narcs = Y.length;
        int G = Y[0].length;
        int ny = Y[0][0].length;
        int S = Ys[0][0].length;
        int M = spec.nstations();
        int MR = x0v.length;
        int src = spec.arcs[e][0];
        int cls = spec.arcs[e][2];
        double lam = rateMean(spec, alpha, beta, e);
        double loglam = rateLogMean(spec, alpha, beta, e);
        int kclass = src > 0 ? (cls - 1) * M + src - 1 : -1;
        double[][] ge = new double[G][ny];
        double[][] he = new double[G][ny];
        double[][] obsw = new double[spec.obsTimes.length][ny];
        for (int q = 0; q < obsw.length; q++) {
            for (int y = 0; y < ny; y++) {
                obsw[q][y] = 1.0;
            }
        }
        double sgnOwnClass = kclass >= 0 ? sgnClass[e][kclass] : 0.0;
        double sgnOwnStat = src > 0 ? sgnStat[e][src - 1] : 0.0;
        double[] aAll = new double[MR];
        for (int g = 0; g < G; g++) {
            boolean hasObs = false;
            if (wantObs) {
                for (int q = 0; q < obsIdx.length; q++) {
                    if (obsIdx[q] == g) {
                        hasObs = true;
                        break;
                    }
                }
            }
            double[][] aStore = hasObs ? new double[MR][S] : null;
            double sumUps;
            double sumLog;
            double[] gy = ge[g];
            double[] hy = he[g];
            for (int y = 0; y < ny; y++) {
                gy[y] = 0.0;
                hy[y] = 0.0;
            }
            for (int s = 0; s < S; s++) {
                double a = kclass >= 0 ? x0v[kclass] : 0.0;
                double b = src > 0 ? x0s[src - 1] : 0.0;
                for (int f = 0; f < narcs; f++) {
                    if (f == e) {
                        continue;
                    }
                    if (kclass >= 0) {
                        a += sgnClass[f][kclass] * Ys[f][g][s];
                        b += sgnStat[f][src - 1] * Ys[f][g][s];
                    }
                }
                if (hasObs) {
                    for (int k = 0; k < MR; k++) {
                        double acc = x0v[k];
                        for (int f = 0; f < narcs; f++) {
                            if (f != e) {
                                acc += sgnClass[f][k] * Ys[f][g][s];
                            }
                        }
                        aStore[k][s] = acc;
                    }
                }
                for (int y = 0; y < ny; y++) {
                    double u;
                    if (kclass >= 0) {
                        u = ups(a + sgnOwnClass * yvec[y], b + sgnOwnStat * yvec[y],
                                arcServers[e], arcSched[e], spec.capacity[kclass],
                                capStat[src - 1]);
                    } else {
                        u = 1.0;
                    }
                    gy[y] += u;
                    hy[y] += Math.log(u + opt.delta / lam);
                }
            }
            for (int y = 0; y < ny; y++) {
                sumUps = gy[y] / S;
                sumLog = hy[y] / S;
                gy[y] = opt.delta + lam * sumUps;
                hy[y] = Math.exp(loglam + sumLog);
            }
            if (hasObs) {
                for (int q = 0; q < obsIdx.length; q++) {
                    if (obsIdx[q] != g) {
                        continue;
                    }
                    obsWeight(spec, opt, spec.obsData[q], aStore, sgnClass[e], yvec, obsw[q]);
                }
            }
        }
        double[][][] out = new double[3][][];
        out[0] = ge;
        out[1] = he;
        out[2] = obsw;
        return out;
    }

    /** Multiplicative jump carried by an observation in the backward pass. */
    static void obsWeight(VariationalSpec spec, VariationalOptions opt, double[] obsRow,
                          double[][] aAll, double[] sgnE, double[] yvec, double[] w) {
        int ny = yvec.length;
        int S = aAll[0].length;
        double[] acc = new double[ny];
        for (int k = 0; k < obsRow.length; k++) {
            if (Double.isNaN(obsRow[k])) {
                continue;
            }
            double range = Math.max(1.0, spec.obsRange[k]);
            for (int y = 0; y < ny; y++) {
                double a = 0.0;
                for (int s = 0; s < S; s++) {
                    double x = aAll[k][s] + sgnE[k] * yvec[y];
                    double p;
                    if (x == obsRow[k]) {
                        p = 1.0 - spec.epsilon;
                    } else if (x >= 0 && x <= spec.obsRange[k]) {
                        p = spec.epsilon / range;
                    } else {
                        p = 0.0;
                    }
                    a += Math.log(opt.floor + p);
                }
                acc[y] += a / S;
            }
        }
        for (int y = 0; y < ny; y++) {
            w[y] = Math.exp(acc[y]);
        }
    }

    /**
     * Evidence lower bound: path term, observation term and the divergence of
     * the rate posteriors from their priors.
     */
    static double bound(VariationalSpec spec, VariationalOptions opt, double[] alpha,
                        double[] beta, double[][][] Y, double[][][] nu, double[][][] gexp,
                        double[][][] hexp, double[][][] Ys, double[][] sgnClass, double[] x0v,
                        int[] obsIdx, double dt) {
        int narcs = Y.length;
        int G = Y[0].length;
        int ny = Y[0][0].length;
        int S = Ys[0][0].length;
        int MR = x0v.length;
        double b = 0.0;
        double[] acc = new double[G];
        for (int e = 0; e < narcs; e++) {
            for (int g = 0; g < G; g++) {
                double s = 0.0;
                for (int y = 0; y < ny; y++) {
                    double n = nu[e][g][y];
                    double term = n - gexp[e][g][y];
                    if (n > 0) {
                        term -= n * Math.log(n / Math.max(opt.floor, hexp[e][g][y]));
                    }
                    s += Y[e][g][y] * term;
                }
                acc[g] = s;
            }
            b += trapz(acc, dt);
        }
        for (int k = 0; k < spec.obsTimes.length; k++) {
            int g = obsIdx[k];
            double tot = 0.0;
            for (int s = 0; s < S; s++) {
                double a = 0.0;
                for (int j = 0; j < MR; j++) {
                    if (Double.isNaN(spec.obsData[k][j])) {
                        continue;
                    }
                    double x = x0v[j];
                    for (int f = 0; f < narcs; f++) {
                        x += sgnClass[f][j] * Ys[f][g][s];
                    }
                    double p;
                    if (x == spec.obsData[k][j]) {
                        p = 1.0 - spec.epsilon;
                    } else if (x >= 0 && x <= spec.obsRange[j]) {
                        p = spec.epsilon / Math.max(1.0, spec.obsRange[j]);
                    } else {
                        p = 0.0;
                    }
                    a += Math.log(opt.floor + p);
                }
                tot += a;
            }
            b += tot / S;
        }
        for (int p = 0; p < alpha.length; p++) {
            b -= klGamma(alpha[p], beta[p], spec.alpha0[p], spec.beta0[p]);
        }
        return b;
    }
}
