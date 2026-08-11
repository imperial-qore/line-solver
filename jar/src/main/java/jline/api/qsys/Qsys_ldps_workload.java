/**
 * @file Quantity of work in a load-dependent processor sharing station with blocking.
 *
 * Stationary distribution of the workload of the single-stage generalized
 * processor sharing model with Poisson arrivals and blocking, per J.W. Cohen,
 * "The multiple phase service network with generalized processor sharing",
 * Acta Informatica 12, 245-284 (1979), Sect. 9, eqs. (9.1)-(9.3).
 *
 * @since LINE 3.0
 */
package jline.api.qsys;

import jline.lang.processes.Distribution;

public final class Qsys_ldps_workload {
    private Qsys_ldps_workload() {}

    /** Default number of points of the internal quadrature grid. */
    private static final int DEFAULT_NGRID = 2001;

    /**
     * Stationary distribution of the quantity of work in a single-stage
     * load-dependent processor sharing station with Poisson arrivals and
     * blocking.
     *
     * Assumed model (Cohen 1979, Sect. 9; the model of Sect. 7 with one stage):
     * a single service stage fed by a Poisson arrival stream of rate lambda; a
     * blocking capacity N, so that a request arriving when N requests are
     * already present is lost and leaves no trace on the state; generalized
     * processor sharing, so that when x requests are present each accrues
     * service at rate f(x) and the stage completes work at total rate x*f(x);
     * and required service times i.i.d. with absolutely continuous distribution
     * B of finite mean beta. The station is parametrized by the LINE
     * load-dependent total rate scaling alpha(x)=x*f(x), the argument of
     * setLoadDependence at a PS station.
     *
     * With psi the total amount of service still to be given to the requests
     * present, Cohen eqs. (9.1)-(9.3) give
     *
     *   Pr{psi &lt; y} = sum_{h=0}^{N} p_h Psi^{h*}(y),
     *   p_h = (rho^h/h!) phi(h) / sum_k (rho^k/k!) phi(k),   rho = lambda*beta,
     *   phi(h) = 1/prod_{k=1}^{h} f(k),   phi(0)=1,
     *   Psi(y) = int_0^y (1-B(v))/beta dv,
     *
     * with Psi^{h*} the h-fold convolution of Psi and Psi^{0*} degenerate at
     * zero. Substituting f(k)=alpha(k)/k the factorial cancels, leaving
     * p_h proportional to rho^h/prod_{k=1}^{h} alpha(k), the familiar
     * load-dependent birth-death form. Psi is the equilibrium (residual life)
     * distribution of B, so psi is a mixture of h-fold convolutions of residual
     * service times with an atom p_0 at zero.
     *
     * This is the model of Cohen (1979) Sect. 9 only. It is not the weighted
     * GPS/DPS discipline of SchedStrategy.GPS, whose per-class weights this
     * formula does not represent.
     *
     * @param lambda rate of the Poisson arrival stream (finite, positive)
     * @param B      required service time distribution (continuous, finite
     *               positive mean)
     * @param alpha  rate scaling alpha(n)=n*f(n) for n=1..N (finite, positive)
     * @param N      blocking capacity (finite positive integer)
     * @param t      grid at which the CDF is returned, or null for an
     *               automatically sized grid
     * @param ngrid  number of points of the internal uniform quadrature grid,
     *               or 0 for the default of 2001. Accuracy is second order in
     *               the step for an absolutely continuous B, the case Cohen
     *               assumes, and falls back to first order when B has an atom
     *               so that 1-B is discontinuous; raise ngrid for those.
     * @return the workload CDF, the grid it is reported on, and the stationary
     *         number in system
     */
    public static QsysWorkloadResult qsys_ldps_workload(double lambda, Distribution B, double[] alpha,
                                                        int N, double[] t, int ngrid) {
        // gating: reject anything outside the assumed model
        if (!Double.isFinite(lambda) || lambda <= 0) {
            throw new IllegalArgumentException(
                    "lambda must be a finite positive scalar, the rate of the Poisson arrival stream.");
        }
        if (B == null) {
            throw new IllegalArgumentException("B must be a Distribution giving the required service time.");
        }
        if (B.isDisabled() || B.isImmediate()) {
            throw new IllegalArgumentException(
                    "B must be an active service time distribution, not Disabled or Immediate.");
        }
        if (!B.isContinuous()) {
            throw new IllegalArgumentException("B must be a continuous distribution: Cohen (1979) Sect. 9 assumes "
                    + "an absolutely continuous required service time.");
        }
        double beta = B.getMean();
        if (!Double.isFinite(beta) || beta <= 0) {
            throw new IllegalArgumentException("B must have a finite positive mean.");
        }
        if (N < 1) {
            throw new IllegalArgumentException(
                    "N must be a finite positive integer, the blocking capacity of the service stage.");
        }
        if (alpha == null || alpha.length < N) {
            throw new IllegalArgumentException("alpha must supply the rate scaling for n=1.." + N + ", but only "
                    + (alpha == null ? 0 : alpha.length) + " entries were given.");
        }
        for (int k = 0; k < N; k++) {
            if (!Double.isFinite(alpha[k]) || alpha[k] <= 0) {
                throw new IllegalArgumentException("alpha(n) must be finite and strictly positive for n=1..N, since "
                        + "every request in a busy stage is served at a positive rate.");
            }
        }
        if (ngrid == 0) {
            ngrid = DEFAULT_NGRID;
        }
        if (ngrid < 2) {
            throw new IllegalArgumentException("ngrid must be an integer of at least 2.");
        }

        // stationary number in system, eq. (9.1). Accumulated in logs so that
        // large rho or large N do not overflow before normalization.
        double rho = lambda * beta;
        double[] logw = new double[N + 1];
        logw[0] = 0.0;
        for (int k = 1; k <= N; k++) {
            logw[k] = logw[k - 1] + Math.log(rho) - Math.log(alpha[k - 1]);
        }
        double mx = logw[0];
        for (int k = 1; k <= N; k++) {
            if (logw[k] > mx) {
                mx = logw[k];
            }
        }
        double[] p = new double[N + 1];
        double sum = 0.0;
        for (int k = 0; k <= N; k++) {
            p[k] = Math.exp(logw[k] - mx);
            sum += p[k];
        }
        for (int k = 0; k <= N; k++) {
            p[k] /= sum;
        }

        // see _kb/03-api-layer.md for rationale
        double m1e = beta * (1 + B.getSCV()) / 2;
        if (!Double.isFinite(m1e) || m1e <= 0) {
            throw new IllegalArgumentException("B must have a finite second moment: the equilibrium residual "
                    + "service time is otherwise undefined.");
        }
        boolean userGrid = t != null && t.length > 0;
        double tmax;
        if (userGrid) {
            tmax = 0.0;
            for (int j = 0; j < t.length; j++) {
                if (!Double.isFinite(t[j]) || t[j] < 0) {
                    throw new IllegalArgumentException("t must be a vector of finite non-negative times.");
                }
                if (t[j] > tmax) {
                    tmax = t[j];
                }
            }
            if (tmax <= 0) {
                tmax = m1e;
            }
        } else {
            int hmax = 0;
            for (int h = N; h >= 0; h--) {
                if (p[h] > 1e-12) {
                    hmax = h;
                    break;
                }
            }
            if (hmax < 1) {
                hmax = 1;
            }
            tmax = m1e * (hmax + 8 * Math.sqrt(hmax));
            tmax = Math.max(tmax, 8 * m1e);
        }

        double[] tg = new double[ngrid];
        for (int j = 0; j < ngrid; j++) {
            tg[j] = tmax * j / (ngrid - 1.0);
        }
        double dt = tg[1] - tg[0];

        // equilibrium residual service distribution, eq. (9.3): density (1-B(v))/beta
        double[] e = new double[ngrid];
        for (int j = 0; j < ngrid; j++) {
            e[j] = (1 - B.evalCDF(tg[j])) / beta;
        }

        // workload distribution, eq. (9.2). The h=0 term is degenerate at zero,
        // contributing the atom p_0 over the whole non-negative grid.
        double[] Fg = new double[ngrid];
        for (int j = 0; j < ngrid; j++) {
            Fg[j] = p[0];
        }
        double[] dens = new double[ngrid];
        System.arraycopy(e, 0, dens, 0, ngrid);
        for (int h = 1; h <= N; h++) {
            if (h > 1) {
                dens = convTrap(dens, e, dt);
            }
            double[] Psih = cumTrapz(tg, dens);
            for (int j = 0; j < ngrid; j++) {
                Fg[j] += p[h] * Psih[j];
            }
        }

        if (userGrid) {
            return new QsysWorkloadResult(interp1Linear(tg, Fg, t), t, p);
        }
        return new QsysWorkloadResult(Fg, tg, p);
    }

    /**
     * Convenience overload using the default grid and quadrature resolution.
     *
     * @param lambda rate of the Poisson arrival stream
     * @param B      required service time distribution
     * @param alpha  rate scaling alpha(n)=n*f(n) for n=1..N
     * @param N      blocking capacity
     * @return the workload CDF, its grid, and the stationary number in system
     */
    public static QsysWorkloadResult qsys_ldps_workload(double lambda, Distribution B, double[] alpha, int N) {
        return qsys_ldps_workload(lambda, B, alpha, N, null, 0);
    }

    /**
     * Convolution of two densities sampled on a uniform grid, using the
     * trapezoidal rule rather than the rectangle rule implied by a bare
     * convolution. The correction matters here because the equilibrium density
     * does not vanish at the origin: e(0)=1/beta. Writing t_i=i*dt,
     *
     *   (f*g)(t_i) = int_0^{t_i} f(v) g(t_i-v) dv
     *              ~ dt*[ sum_{j=0}^{i} f_j g_{i-j} - (f_0 g_i + f_i g_0)/2 ],
     *
     * i.e. the raw convolution less half of each endpoint. Without the
     * correction each convolution over-counts by dt*f_0*g_i, which accumulates
     * over h and drives the mixture CDF above one.
     */
    private static double[] convTrap(double[] f, double[] g, double dt) {
        int n = f.length;
        double[] c = new double[n];
        for (int i = 0; i < n; i++) {
            double acc = 0.0;
            for (int j = 0; j <= i; j++) {
                acc += f[j] * g[i - j];
            }
            c[i] = dt * (acc - (f[0] * g[i] + f[i] * g[0]) / 2);
        }
        return c;
    }

    /** Cumulative trapezoidal integral of y over the grid x. */
    private static double[] cumTrapz(double[] x, double[] y) {
        int n = x.length;
        double[] c = new double[n];
        c[0] = 0.0;
        for (int j = 1; j < n; j++) {
            c[j] = c[j - 1] + (x[j] - x[j - 1]) * (y[j] + y[j - 1]) / 2;
        }
        return c;
    }

    /** Linear interpolation of (x,y) onto the query points xq. */
    private static double[] interp1Linear(double[] x, double[] y, double[] xq) {
        double[] yq = new double[xq.length];
        int n = x.length;
        for (int i = 0; i < xq.length; i++) {
            double v = xq[i];
            if (v <= x[0]) {
                yq[i] = y[0];
                continue;
            }
            if (v >= x[n - 1]) {
                yq[i] = y[n - 1];
                continue;
            }
            int lo = 0;
            int hi = n - 1;
            while (hi - lo > 1) {
                int mid = (lo + hi) >>> 1;
                if (x[mid] <= v) {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            double w = (v - x[lo]) / (x[hi] - x[lo]);
            yq[i] = y[lo] * (1 - w) + y[hi] * w;
        }
        return yq;
    }
}
