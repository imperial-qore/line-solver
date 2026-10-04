/**
 * @file Subtask dispersion and the delays that minimise it
 *
 * Ports of matlab/src/api/fj/fj_dispersion.m and fj_delay_opt.m, Section 6.2 of
 * A. Thomasian, "Analysis of Fork/Join and Related Queueing Systems", ACM
 * Computing Surveys 47(2), Article 17, 2014 (Eqs. (24)-(25)).
 *
 * @since LINE 3.0
 */
package jline.api.fj;

public final class FJ_dispersion {
    private FJ_dispersion() {}

    /** [Edisp, Emax, Emin] of fj_dispersion. */
    public static final class FJDispersionResult {
        public final double Edisp;
        public final double Emax;
        public final double Emin;

        public FJDispersionResult(double Edisp, double Emax, double Emin) {
            this.Edisp = Edisp;
            this.Emax = Emax;
            this.Emin = Emin;
        }
    }

    /**
     * Mean subtask dispersion of a split-merge system whose branches are shifted
     * Erlangs, the equivalent used in the delay-scheduling construction:
     *
     * E[D_d] = integral_0^inf [ 1 - prod_i F_i(x-d_i) - prod_i (1-F_i(x-d_i)) ] dx.
     *
     * That integrand is non-negative and vanishes at both ends; the difference
     * of the two products printed in the survey is not the dispersion and can go
     * negative.
     */
    public static FJDispersionResult fj_dispersion(int[] shape, double[] rate, double[] d,
                                                   double tol, int npanels) {
        int N = shape.length;
        if (rate.length != N || d.length != N) {
            throw new IllegalArgumentException("shape, rate and d must have the same length.");
        }
        if (N < 1) {
            throw new IllegalArgumentException("At least one branch is required.");
        }
        for (int i = 0; i < N; i++) {
            if (shape[i] < 1) {
                throw new IllegalArgumentException("Erlang stage counts must be positive.");
            }
            if (!(rate[i] > 0)) {
                throw new IllegalArgumentException("Erlang stage rates must be positive.");
            }
            if (d[i] < 0) {
                throw new IllegalArgumentException("Delays must be non-negative.");
            }
        }
        if (npanels % 2 != 0) {
            npanels++;
        }
        double dmax = d[0], mmax = (double) shape[0] / rate[0];
        for (int i = 0; i < N; i++) {
            if (d[i] > dmax) {
                dmax = d[i];
            }
            double mi = (double) shape[i] / rate[i];
            if (mi > mmax) {
                mmax = mi;
            }
        }
        double U = dmax + 8 * mmax;
        for (int it = 0; it < 60; it++) {
            double prodF = 1;
            for (int i = 0; i < N; i++) {
                prodF *= erlangCdf(U - d[i], shape[i], rate[i]);
            }
            if (1 - prodF < tol) {
                break;
            }
            U *= 2;
        }
        double h = U / npanels;
        double accMax = 0, accMin = 0;
        for (int i = 0; i <= npanels; i++) {
            double x = h * i;
            double Fprod = 1, Sprod = 1;
            for (int j = 0; j < N; j++) {
                double Fj = erlangCdf(x - d[j], shape[j], rate[j]);
                Fprod *= Fj;
                Sprod *= (1 - Fj);
            }
            double w = (i == 0 || i == npanels) ? 1 : ((i % 2 == 1) ? 4 : 2);
            accMax += w * (1 - Fprod);
            accMin += w * Sprod;
        }
        double scale = h / 3;
        double Emax = scale * accMax;
        double Emin = scale * accMin;
        return new FJDispersionResult(Emax - Emin, Emax, Emin);
    }

    public static FJDispersionResult fj_dispersion(int[] shape, double[] rate) {
        return fj_dispersion(shape, rate, new double[shape.length], 1e-10, 4000);
    }

    /** Erlang-k distribution function, zero on the negative half line. */
    private static double erlangCdf(double t, int k, double mu) {
        if (!(t > 0)) {
            return 0;
        }
        double acc = 0, term = 1;
        for (int j = 0; j < k; j++) {
            if (j > 0) {
                term = term * (mu * t) / j;
            }
            acc += term;
        }
        return 1 - Math.exp(-mu * t) * acc;
    }

    /** [d, Edisp, Emax] of fj_delay_opt. */
    public static final class FJDelayOptResult {
        public final double[] d;
        public final double Edisp;
        public final double Emax;

        public FJDelayOptResult(double[] d, double Edisp, double Emax) {
            this.d = d;
            this.Edisp = Edisp;
            this.Emax = Emax;
        }
    }

    /**
     * The deterministic delays that minimise the mean dispersion. Holding back a
     * fast branch costs little at the last completion and buys a great deal at
     * the first, so the minimiser is generally interior and strictly positive on
     * every branch but the slowest.
     *
     * Cyclic coordinate descent with a golden section line search on each
     * coordinate: deterministic, derivative-free, and the same sequence of
     * evaluations in all four codebases. Adding a constant to every delay shifts
     * both order statistics equally, so the representative with min(d) = 0 is
     * returned.
     */
    public static FJDelayOptResult fj_delay_opt(int[] shape, double[] rate, int maxsweeps,
                                                double dtol, int npanels) {
        int N = shape.length;
        if (rate.length != N) {
            throw new IllegalArgumentException("shape and rate must have the same length.");
        }
        if (N < 1) {
            throw new IllegalArgumentException("At least one branch is required.");
        }
        double tol = 1e-10;
        double[] d = new double[N];
        if (N < 2) {
            FJDispersionResult r0 = fj_dispersion(shape, rate, d, tol, npanels);
            return new FJDelayOptResult(d, r0.Edisp, r0.Emax);
        }
        double mmax = 0, smax = 0;
        for (int i = 0; i < N; i++) {
            double mi = (double) shape[i] / rate[i];
            double si = Math.sqrt(shape[i]) / rate[i];
            if (mi > mmax) {
                mmax = mi;
            }
            if (si > smax) {
                smax = si;
            }
        }
        double ub = mmax + 8 * smax;
        double invphi = (Math.sqrt(5.0) - 1) / 2;
        double fcur = fj_dispersion(shape, rate, d, tol, npanels).Edisp;
        for (int sweep = 0; sweep < maxsweeps; sweep++) {
            double fprev = fcur;
            for (int i = 0; i < N; i++) {
                double a = 0, b = ub;
                double c = b - invphi * (b - a), dd = a + invphi * (b - a);
                double[] probe = d.clone();
                probe[i] = c;
                double fc = fj_dispersion(shape, rate, probe, tol, npanels).Edisp;
                probe[i] = dd;
                double fd = fj_dispersion(shape, rate, probe, tol, npanels).Edisp;
                for (int it = 0; it < 60; it++) {
                    if (fc < fd) {
                        b = dd;
                        dd = c;
                        fd = fc;
                        c = b - invphi * (b - a);
                        probe[i] = c;
                        fc = fj_dispersion(shape, rate, probe, tol, npanels).Edisp;
                    } else {
                        a = c;
                        c = dd;
                        fc = fd;
                        dd = a + invphi * (b - a);
                        probe[i] = dd;
                        fd = fj_dispersion(shape, rate, probe, tol, npanels).Edisp;
                    }
                    if ((b - a) <= dtol * Math.max(1.0, ub)) {
                        break;
                    }
                }
                d[i] = (fc < fd) ? c : dd;
            }
            double dmin = d[0];
            for (int i = 1; i < N; i++) {
                if (d[i] < dmin) {
                    dmin = d[i];
                }
            }
            for (int i = 0; i < N; i++) {
                d[i] -= dmin;
            }
            fcur = fj_dispersion(shape, rate, d, tol, npanels).Edisp;
            if (Math.abs(fprev - fcur) <= dtol * Math.max(1.0, Math.abs(fprev))) {
                break;
            }
        }
        FJDispersionResult r = fj_dispersion(shape, rate, d, tol, npanels);
        return new FJDelayOptResult(d, r.Edisp, r.Emax);
    }

    public static FJDelayOptResult fj_delay_opt(int[] shape, double[] rate) {
        return fj_delay_opt(shape, rate, 40, 1e-8, 2000);
    }
}
