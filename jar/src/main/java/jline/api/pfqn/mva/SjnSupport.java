/**
 * @file Shared internals of the shortest-job-next response time equations
 *
 * The two-moment Erlang-mixture reconstruction of the job size distribution, the fixed-grid
 * quadrature, the conditional waiting time recursion at one station and the utilization cap,
 * shared by Pfqn_mvasjn (population lattice) and Pfqn_amvasjn (Schweitzer fixed point) so the
 * exact and approximate routes cannot drift apart. Ported at parity from the MATLAB
 * matlab/src/api/pfqn/private/sjn_*.m helpers.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.mva;

import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.FastMath;

final class SjnSupport {
    private SjnSupport() {}

    /** Erlang mixture fitted to a mean and a squared coefficient of variation. */
    static final class Fit {
        final double[] w;
        final int[] k;
        final double[] mu;

        Fit(double[] w, int[] k, double[] mu) {
            this.w = w;
            this.k = k;
            this.mu = mu;
        }

        boolean isEmpty() {
            return w.length == 0;
        }
    }

    /**
     * Two-moment fit of a service time distribution: a branching Erlang (Erlang(k-1) and
     * Erlang(k) sharing a rate) below CV^2 = 1 and a balanced-means hyperexponential above it.
     * The mixture form is what makes theta(x) and the tail integrals closed form.
     */
    static Fit fit(double s, double cv2) {
        if (s <= 0) {
            return new Fit(new double[0], new int[0], new double[0]);
        }
        if (cv2 < 0) {
            throw new IllegalArgumentException("negative squared coefficient of variation");
        }
        if (FastMath.abs(cv2 - 1) < 1e-8) {
            return new Fit(new double[] {1.0}, new int[] {1}, new double[] {1.0 / s});
        }
        if (cv2 < 1) {
            int k = (int) FastMath.ceil(1.0 / cv2);
            double p = (k * cv2 - FastMath.sqrt(k * (1 + cv2) - k * k * cv2)) / (1 + cv2);
            double mu = (k - p) / s;
            return new Fit(new double[] {p, 1 - p}, new int[] {k - 1, k}, new double[] {mu, mu});
        }
        double p = 0.5 * (1 + FastMath.sqrt((cv2 - 1) / (cv2 + 1)));
        return new Fit(new double[] {p, 1 - p}, new int[] {1, 1},
                new double[] {2 * p / s, 2 * (1 - p) / s});
    }

    /** Density of the mixture on a grid. */
    static double[] pdf(Fit f, double[] x) {
        double[] y = new double[x.length];
        for (int j = 0; j < f.w.length; j++) {
            int k = f.k[j];
            double mu = f.mu[j];
            for (int i = 0; i < x.length; i++) {
                // realmin, the smallest normal, not the smallest denormal
                double xi = FastMath.max(x[i], Double.MIN_NORMAL);
                y[i] += f.w[j] * FastMath.exp(k * FastMath.log(mu) + (k - 1) * FastMath.log(xi)
                        - mu * x[i] - Gamma.logGamma(k));
            }
        }
        return y;
    }

    /** The primitive int_0^x t f(t) dt, in closed form. */
    static double[] theta(Fit f, double[] x) {
        double[] y = new double[x.length];
        for (int j = 0; j < f.w.length; j++) {
            int k = f.k[j];
            double mu = f.mu[j];
            for (int i = 0; i < x.length; i++) {
                y[i] += f.w[j] * (k / mu) * Gamma.regularizedGammaP(k + 1, mu * x[i]);
            }
        }
        return y;
    }

    /** The tail mass int_x^inf f(t) dt, in closed form. */
    static double ccdf(Fit f, double x) {
        double y = 0;
        for (int j = 0; j < f.w.length; j++) {
            y += f.w[j] * Gamma.regularizedGammaQ(f.k[j], f.mu[j] * x);
        }
        return y;
    }

    /**
     * The tail integral int_Lx^inf t^order exp(-c (t-Lx)) f(t) dt, evaluated in logarithms so
     * that exp(c Lx) cannot overflow against an underflowing incomplete gamma.
     */
    static double tailmom(Fit f, double Lx, double c, int order) {
        double y = 0;
        for (int j = 0; j < f.w.length; j++) {
            int k = f.k[j];
            double mu = f.mu[j];
            double rate = mu + c;
            double g = Gamma.regularizedGammaQ(k + order, rate * Lx);
            if (g <= 0) {
                continue;
            }
            double lg = c * Lx + k * FastMath.log(mu / rate) + FastMath.log(g);
            if (order == 1) {
                lg += FastMath.log(k / rate);
            }
            y += f.w[j] * FastMath.exp(lg);
        }
        return y;
    }

    /** Composite Simpson over an even number of subdivisions. */
    static double simpson(double[] y, double dx) {
        int n = y.length;
        double odd = 0;
        double even = 0;
        for (int i = 1; i < n - 1; i += 2) {
            odd += y[i];
        }
        for (int i = 2; i < n - 1; i += 2) {
            even += y[i];
        }
        return dx / 3 * (y[0] + y[n - 1] + 4 * odd + 2 * even);
    }

    /**
     * Cumulative Simpson: full panels at the odd nodes and a half panel at the even ones, so the
     * primitive is available at every grid node. Quadrature at arbitrary abscissae could not
     * provide it, the profile being needed again at the next population step.
     */
    static double[] cumsimpson(double[] y, double dx) {
        int n = y.length;
        double[] I = new double[n];
        for (int i = 2; i < n; i += 2) {
            I[i] = I[i - 2] + dx / 3 * (y[i - 2] + 4 * y[i - 1] + y[i]);
        }
        for (int i = 1; i < n; i += 2) {
            if (i + 1 < n) {
                I[i] = I[i - 1] + dx / 12 * (5 * y[i - 1] + 8 * y[i] - y[i + 1]);
            } else {
                I[i] = I[i - 1] + dx / 12 * (-y[i - 2] + 8 * y[i - 1] + 5 * y[i]);
            }
        }
        return I;
    }

    /** Job size grid of one SJN station and the population-independent integrals over it. */
    static final class Grid {
        final double Lx;
        final double dx;
        final double[] x;
        final double[][] f;      // (ngrid x R)
        final double[][] theta;  // (ngrid x R)
        final double[] tail0;
        final double[] tail1;
        final Fit[] fit;

        Grid(double Lx, double dx, double[] x, double[][] f, double[][] theta,
                double[] tail0, double[] tail1, Fit[] fit) {
            this.Lx = Lx;
            this.dx = dx;
            this.x = x;
            this.f = f;
            this.theta = theta;
            this.tail0 = tail0;
            this.tail1 = tail1;
            this.fit = fit;
        }
    }

    /**
     * Build the grid of one station. It spans [0, Lfactor * max_r s_r] because the conditional
     * waiting time has flattened out well before that point, its remainder being carried by the
     * analytic tail of the recursion rather than by quadrature.
     */
    static Grid setup(double[] S, double[] scv, int ns, double Lfactor) {
        int R = S.length;
        double smax = 0;
        for (int r = 0; r < R; r++) {
            smax = FastMath.max(smax, S[r]);
        }
        if (smax <= 0) {
            throw new IllegalArgumentException("the station has zero service demand in every class");
        }
        double Lx = Lfactor * smax;
        int ngrid = ns + 1;
        double[] x = new double[ngrid];
        double dx = Lx / ns;
        for (int i = 0; i < ngrid; i++) {
            x[i] = i * dx;
        }
        double[][] f = new double[ngrid][R];
        double[][] th = new double[ngrid][R];
        double[] tail0 = new double[R];
        double[] tail1 = new double[R];
        Fit[] fits = new Fit[R];
        for (int r = 0; r < R; r++) {
            fits[r] = fit(S[r], scv[r]);
            if (fits[r].isEmpty()) {
                continue;
            }
            double[] fr = pdf(fits[r], x);
            double[] tr = theta(fits[r], x);
            for (int i = 0; i < ngrid; i++) {
                f[i][r] = fr[i];
                th[i][r] = tr[i];
            }
            tail0[r] = ccdf(fits[r], Lx);
            tail1[r] = S[r] - tr[ngrid - 1];
        }
        return new Grid(Lx, dx, x, f, th, tail0, tail1, fits);
    }

    /** State of one SJN station at the reference population. */
    static final class State {
        final double[] lam;
        final double[] U;
        final double[] Q;
        final double[][] W;    // (ngrid x R)
        final double[][] phi;  // (ngrid x R)
        final double[] phiinf;

        State(double[] lam, double[] U, double[] Q, double[][] W, double[][] phi, double[] phiinf) {
            this.lam = lam;
            this.U = U;
            this.Q = Q;
            this.W = W;
            this.phi = phi;
            this.phiinf = phiinf;
        }
    }

    /** Outcome of one evaluation of the conditional waiting time equation. */
    static final class StationResult {
        final double C;
        final double[] W;
        final double[] phi;
        final double phiinf;
        final double[] tail;

        StationResult(double C, double[] W, double[] phi, double phiinf, double[] tail) {
            this.C = C;
            this.W = W;
            this.phi = phi;
            this.phiinf = phiinf;
            this.tail = tail;
        }
    }

    /**
     * One evaluation of the SJN conditional waiting time equation for a tagged customer of class
     * r at station m.
     *
     * <p>The quantity lam_k W_k(x) f_k(x) is the density, in the job size x, of the queued
     * class-k customers, so deflating it by beta_k is what turns the same equation into either
     * the exact recursion (beta = 1, the state already being the one at n - e_r) or the
     * Schweitzer closure (beta_r = (N_r-1)/N_r, the state being the one at N).</p>
     */
    static StationResult station(int m, int r, Grid G, double[] S, double[] scv, double[] V,
            State st, double[] beta, boolean useprio, int[] prio) {
        int R = S.length;
        int ngrid = G.x.length;
        double[] lamb = new double[R];
        double[] Ub = new double[R];
        double[] Qb = new double[R];
        double RL = 0;
        for (int k = 0; k < R; k++) {
            lamb[k] = beta[k] * st.lam[k];
            Ub[k] = beta[k] * st.U[k];
            Qb[k] = beta[k] * st.Q[k];
            RL += (1 + scv[k]) * S[k] * Ub[k] / 2;
        }
        double[] num = new double[ngrid];
        double[] den = new double[ngrid];
        double numinf;
        double deninf;
        if (useprio) {
            double base = RL;
            double uhi = 0;
            for (int k = 0; k < R; k++) {
                if (prio[k] < prio[r]) {
                    base += S[k] * (Qb[k] - Ub[k]);
                    uhi += Ub[k];
                }
            }
            for (int i = 0; i < ngrid; i++) {
                num[i] = base + lamb[r] * st.phi[i][r];
                den[i] = 1 - uhi - lamb[r] * G.theta[i][r];
            }
            numinf = base + lamb[r] * st.phiinf[r];
            deninf = 1 - uhi - lamb[r] * S[r];
        } else {
            double usum = 0;
            double phiinfsum = 0;
            for (int k = 0; k < R; k++) {
                usum += lamb[k] * S[k];
                phiinfsum += lamb[k] * st.phiinf[k];
            }
            for (int i = 0; i < ngrid; i++) {
                double n = RL;
                double d = 1;
                for (int k = 0; k < R; k++) {
                    n += lamb[k] * st.phi[i][k];
                    d -= lamb[k] * G.theta[i][k];
                }
                num[i] = n;
                den[i] = d;
            }
            numinf = RL + phiinfsum;
            deninf = 1 - usum;
        }
        for (int i = 0; i < ngrid; i++) {
            if (den[i] <= 0) {
                throw new RuntimeException(sjnSingularMessage(m));
            }
        }
        if (deninf <= 0) {
            throw new RuntimeException(sjnSingularMessage(m));
        }
        double[] W = new double[ngrid];
        for (int i = 0; i < ngrid; i++) {
            W[i] = num[i] / den[i];
        }
        double Winf = numinf / deninf;
        double slope;
        if (useprio) {
            slope = G.Lx * lamb[r] * G.f[ngrid - 1][r] * (st.W[ngrid - 1][r] + W[ngrid - 1])
                    / den[ngrid - 1];
        } else {
            double acc = 0;
            for (int k = 0; k < R; k++) {
                acc += lamb[k] * G.f[ngrid - 1][k] * (st.W[ngrid - 1][k] + W[ngrid - 1]);
            }
            slope = G.Lx * acc / den[ngrid - 1];
        }
        double a = Winf;
        double b = Winf - W[ngrid - 1];
        double c;
        if (b <= 0) {
            b = 0;
            c = 0;
        } else if (slope < 0) {
            throw new RuntimeException("the conditional waiting time at SJN station " + m
                    + " decreases in the job size, which the discipline forbids: the recursion has"
                    + " become numerically unstable.");
        } else {
            c = slope / b;
        }
        double[] integrand = new double[ngrid];
        for (int i = 0; i < ngrid; i++) {
            integrand[i] = W[i] * G.x[i] * G.f[i][r];
        }
        double[] phi = cumsimpson(integrand, G.dx);
        double phiinf = phi[ngrid - 1] + a * G.tail1[r] - b * tailmom(G.fit[r], G.Lx, c, 1);
        for (int i = 0; i < ngrid; i++) {
            integrand[i] = W[i] * G.f[i][r];
        }
        double Wbar = simpson(integrand, G.dx) + a * G.tail0[r] - b * tailmom(G.fit[r], G.Lx, c, 0);
        double C = V[r] * (S[r] + Wbar);
        return new StationResult(C, W, phi, phiinf, new double[] {a, b, c});
    }

    private static String sjnSingularMessage(int m) {
        return "the SJN recursion at station " + m + " has no solution: the work brought by jobs no"
                + " longer than the tagged one saturates the server, at which point long jobs starve"
                + " and the arrival theorem no longer holds. Reduce the load at that station or"
                + " model it with SolverCTMC or SolverLDES.";
    }

    /** Outcome of the utilization cap. */
    static final class CapResult {
        final double[][] C;
        final double[] X;
        final double[] kappa;
        final boolean bound;

        CapResult(double[][] C, double[] X, double[] kappa, boolean bound) {
            this.C = C;
            this.X = X;
            this.kappa = kappa;
            this.bound = bound;
        }
    }

    /**
     * Enforce U &lt;= umax at every SJN station by inflating its waiting time.
     *
     * <p>The response time equation is an open-system one and has no solution once the fraction
     * of the server taken by jobs no longer than x reaches one. A closed network never reaches it
     * in reality, but the approximation can, because it underestimates the residence time at a
     * congested SJN station and the throughput then exceeds the station capacity. What is imposed
     * is the utilization law sum_r X_r L_mr &lt;= umax, an exact property of the network.</p>
     *
     * <p>The constraint acts on the waiting time, i.e. on the excess C - L, and never on the
     * throughput directly, so that X (Z + sum_m C) = N still holds exactly and no jobs are lost.
     * The same factor scales the station's conditional waiting time profile.</p>
     */
    static CapResult cap(double[][] C, double[][] L, double[] N, double[] Z, int[] sjnset,
            double umax) {
        int nsjn = sjnset.length;
        double[] kappa = new double[nsjn];
        for (int q = 0; q < nsjn; q++) {
            kappa[q] = 1;
        }
        boolean bound = false;
        double[] X = thru(C, N, Z);
        if (nsjn == 0) {
            return new CapResult(C, X, kappa, false);
        }
        if (umax >= 1) {
            throw new IllegalArgumentException("the utilization cap must be strictly below one,"
                    + " the response time equation is singular at one");
        }
        int R = N.length;
        for (int sweep = 0; sweep < 20; sweep++) {
            boolean viol = false;
            for (int q = 0; q < nsjn; q++) {
                int m = sjnset[q];
                double rho = 0;
                for (int r = 0; r < R; r++) {
                    rho += X[r] * L[m][r];
                }
                if (rho <= umax) {
                    continue;
                }
                viol = true;
                bound = true;
                double[] Wq = new double[R];
                for (int r = 0; r < R; r++) {
                    Wq[r] = C[m][r] - L[m][r];
                }
                double hi = 2;
                while (rhoAt(C, L, m, N, Z, Wq, hi) > umax) {
                    hi *= 2;
                    if (hi > 1e12) {
                        throw new RuntimeException("station " + m + " cannot be brought under the"
                                + " utilization cap by any waiting time: its service demands alone"
                                + " saturate it at this population.");
                    }
                }
                double lo = 1;
                for (int b = 0; b < 200; b++) {
                    double mid = (lo + hi) / 2;
                    if (rhoAt(C, L, m, N, Z, Wq, mid) > umax) {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }
                kappa[q] *= hi;
                for (int r = 0; r < R; r++) {
                    C[m][r] = L[m][r] + hi * Wq[r];
                }
                X = thru(C, N, Z);
            }
            if (!viol) {
                return new CapResult(C, X, kappa, bound);
            }
        }
        throw new RuntimeException("the utilization cap did not settle across the SJN stations");
    }

    /** Throughputs implied by the residence times, keeping Little's law exact. */
    static double[] thru(double[][] C, double[] N, double[] Z) {
        int R = N.length;
        double[] X = new double[R];
        for (int r = 0; r < R; r++) {
            if (N[r] <= 0) {
                continue;
            }
            double den = Z[r];
            for (int m = 0; m < C.length; m++) {
                den += C[m][r];
            }
            X[r] = N[r] / den;
        }
        return X;
    }

    private static double rhoAt(double[][] C, double[][] L, int m, double[] N, double[] Z,
            double[] Wq, double kappa) {
        int R = N.length;
        double[] saved = new double[R];
        for (int r = 0; r < R; r++) {
            saved[r] = C[m][r];
            C[m][r] = L[m][r] + kappa * Wq[r];
        }
        double[] X = thru(C, N, Z);
        double rho = 0;
        for (int r = 0; r < R; r++) {
            rho += X[r] * L[m][r];
        }
        for (int r = 0; r < R; r++) {
            C[m][r] = saved[r];
        }
        return rho;
    }
}
