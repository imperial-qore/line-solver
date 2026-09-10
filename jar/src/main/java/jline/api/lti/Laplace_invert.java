/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.lti;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.function.UnaryOperator;

import org.apache.commons.math3.complex.Complex;

import jline.lib.lti.iltcme;
import jline.lib.lti.talbot;

/**
 * Numerical inversion of a Laplace transform.
 *
 * <p>Euler, Talbot, Gaver-Stehfest and CME are delegated to {@link iltcme},
 * whose Abate-Whitt evaluator and coefficient tables are already in the tree.
 * WEEKS (the Laguerre series) is implemented here and had no counterpart in any
 * codebase: {@code jline.lib.lti.laguerre} is Gauss-Laguerre QUADRATURE, which
 * is a different thing, and the native Python {@code lib_lti_laguerre} silently
 * delegated to Talbot.
 *
 * <p>Reference for Weeks: W. Weeks, J. ACM 13, 1966; J. Abate, G. Choudhury and
 * W. Whitt, INFORMS J. Computing 8(4), 1996; P. G. Harrison and
 * W. J. Knottenbelt, "Passage Time Distributions in Large Markov Chains", 2002,
 * Sec. 4.1-4.3, whose Fig. 1 is the automatic scaling search below.
 */
public final class Laplace_invert {

    private Laplace_invert() {
    }

    /** Default number of trapezoid pairs in the Laguerre quadrature. */
    public static final int WEEKS_DEFAULT_P0 = 200;

    /**
     * Laguerre coefficients q_n, n = 0..2*p0-1, of the damped and scaled
     * function f_{sigma,b}(t) = exp(-sigma t) f(t/b), whose generating function
     * is
     *
     * <pre>
     *     Q_{sigma,b}(z) = b/(1-z) * L( b(1+z)/(2(1-z)) + b*sigma ).
     * </pre>
     *
     * <p>NOTE ON THE PAPER. Eq. 10 as printed carries the factor (1-z) rather
     * than 1/(1-z). The scaled form above, printed later in the same section,
     * carries 1/(1-z) and is the correct one: with l_n(t) = exp(-t/2) L_n(t)
     * the transform of l_n is (s-1/2)^n/(s+1/2)^{n+1}, so L(s) = Q(z)/(s+1/2)
     * with z = (s-1/2)/(s+1/2) and s+1/2 = 1/(1-z). Implementing the printed
     * (1-z) is wrong at every t (163 per cent at t = 0.1 on Exp(2)).
     *
     * <p>Sec. 4.3 fixes the trapezoid count at 2*p0 and the radius at
     * r = 0.1^(4/p0) for every n, so the quadrature is one discrete Fourier
     * transform of Q sampled on the circle and the transform is evaluated 2*p0
     * times IN TOTAL rather than per coefficient. The DFT is evaluated
     * directly: at 2*p0 = 400 points that is 160k complex multiplies, which is
     * not worth a dependency.
     */
    public static double[] laplace_weeks_coeffs(UnaryOperator<Complex> F, double sigma, double b,
                                                int p0) {
        if (!(b > 0.0)) {
            throw new RuntimeException("laplace_weeks_coeffs: b must be positive");
        }
        if (p0 <= 0) {
            throw new RuntimeException("laplace_weeks_coeffs: p0 must be positive");
        }
        int N = 2 * p0;
        double r = Math.pow(0.1, 4.0 / p0);
        double twopi = 2.0 * Math.PI;

        Complex[] Q = new Complex[N];
        for (int j = 0; j < N; j++) {
            double u = twopi * j / N;
            Complex z = new Complex(r * Math.cos(u), r * Math.sin(u));
            Complex one = new Complex(1.0, 0.0);
            Complex s = one.add(z).multiply(b).divide(one.subtract(z).multiply(2.0))
                    .add(new Complex(b * sigma, 0.0));
            Q[j] = one.subtract(z).reciprocal().multiply(b).multiply(F.apply(s));
        }

        double[] q = new double[N];
        double rpow = 1.0;
        for (int n = 0; n < N; n++) {
            double re = 0.0;
            double im = 0.0;
            for (int j = 0; j < N; j++) {
                double u = -twopi * ((double) n) * ((double) j) / N;
                double c = Math.cos(u);
                double sn = Math.sin(u);
                re += Q[j].getReal() * c - Q[j].getImaginary() * sn;
                im += Q[j].getReal() * sn + Q[j].getImaginary() * c;
            }
            q[n] = re / N / rpow;
            rpow *= r;
        }
        return q;
    }

    /**
     * The automatic (sigma, b) search of Fig. 1: accept the first pair at which
     * the coefficients have decayed by term p0, doubling sigma from 0.001 and
     * stepping b by 4 whenever sigma passes 0.2.
     *
     * <p>REFUSES BY NAME when the box is exhausted. Raising b further is
     * counterproductive and excessive damping is unstable in finite precision,
     * and a density with a discontinuity in itself or its derivatives has no
     * usable Laguerre representation at all (Sec. 4.2). Returning the last
     * iterate would report noise as an answer; Euler handles those cases.
     */
    public static WeeksParams laplace_weeks_scaling(UnaryOperator<Complex> F, int p0, double tol) {
        double sigma = 0.0;
        double b = 1.0;
        while (true) {
            double[] q = laplace_weeks_coeffs(F, sigma, b, p0);
            if (Math.abs(q[p0]) <= tol && Math.abs(q[p0 + 1]) <= tol) {
                return new WeeksParams(sigma, b, q);
            }
            sigma = (sigma == 0.0) ? 0.001 : 2.0 * sigma;
            if (sigma > 0.2) {
                b = b + 4.0;
                if (b > 10.0) {
                    throw new RuntimeException("laplace_weeks_scaling: no suitable scaling "
                            + "parameters were found for the Laguerre inversion: the transform's "
                            + "density is not smooth enough for a Laguerre series. Use the euler "
                            + "method instead.");
                }
                sigma = 0.0;
            }
        }
    }

    /** The search at the default p0 = 200 and tolerance 1e-10. */
    public static WeeksParams laplace_weeks_scaling(UnaryOperator<Complex> F) {
        return laplace_weeks_scaling(F, WEEKS_DEFAULT_P0, 1e-10);
    }

    /**
     * Truncate at the FIRST index where the coefficients have decayed, never
     * the last. The quadrature divides by r^n with r &lt; 1, so past the genuine
     * decay the entries are rounding noise amplified by r^-n: at n = 2*p0 that
     * factor is 1e8, and scanning for the last entry above a threshold sums
     * 1e-8 of pure noise (worst error on Exp(2) 2.7e-09 instead of 1.9e-14).
     */
    private static int weeksNterms(double[] q) {
        int p0 = q.length / 2;
        for (int n = 1; n + 1 < p0; n++) {
            if (Math.abs(q[n]) <= 1e-13 && Math.abs(q[n + 1]) <= 1e-13) {
                return n;
            }
        }
        return p0;
    }

    /** l_n(t) = exp(-t/2) L_n(t) by the stable recursion of Sec. 4.1. */
    private static double[] laguerreFunctions(double t, int N) {
        double[] l = new double[N];
        if (N == 0) {
            return l;
        }
        l[0] = Math.exp(-t / 2.0);
        if (N > 1) {
            l[1] = (1.0 - t) * l[0];
        }
        for (int n = 2; n < N; n++) {
            l[n] = ((2.0 * n - 1.0 - t) / n) * l[n - 1] - ((n - 1.0) / n) * l[n - 2];
        }
        return l;
    }

    /**
     * Invert by the Laguerre series f(t) = sum_n q_n l_n(t), recovered as
     * exp(sigma*b*t) f_{sigma,b}(b*t).
     *
     * <p>Unlike Euler and Talbot the coefficients do not depend on t, so ONE
     * parameter set serves an arbitrary number of time points: the transform is
     * evaluated 2*p0 times in total, not 2*p0 times per t. That is the property
     * this method is here for, so build the {@link WeeksParams} once and reuse
     * it on a grid.
     */
    public static double laplace_invert_weeks(WeeksParams w, double t) {
        if (!(t > 0.0)) {
            return 0.0;
        }
        int n = weeksNterms(w.q);
        double[] l = laguerreFunctions(w.b * t, n);
        double acc = 0.0;
        for (int i = 0; i < n; i++) {
            acc += w.q[i] * l[i];
        }
        return Math.exp(w.sigma * w.b * t) * acc;
    }

    /** Convenience: build the parameters, then invert at one point. */
    public static double laplace_invert_weeks(UnaryOperator<Complex> F, double t) {
        return laplace_invert_weeks(laplace_weeks_scaling(F), t);
    }

    /**
     * Talbot contour nodes alpha_i, i = 1..n.
     *
     * <p>The nodes and weights themselves come from {@link jline.lib.lti.talbot},
     * which already carried them; only the inverter below was missing, so the
     * class was unreachable. The Euler and Gaver-Stehfest nodes have no
     * counterpart on this surface on purpose: they are built inline by
     * {@code iltcme.abateWhittWeights}, the one implementation the JAR carries
     * for those two. Talbot deforms the contour rather than shifting the real
     * axis, so it is not an Abate-Whitt scheme and that evaluator cannot host it.
     */
    public static Complex[] talbot_get_alpha(int n) {
        return talbot.getalpha(n).toArray(new Complex[0]);
    }

    /** Talbot contour weights omega_i, i = 1..n, for the nodes of ALPHA. */
    public static Complex[] talbot_get_omega(int n, Complex[] alpha) {
        ArrayList<Complex> a = (alpha == null) ? talbot.getalpha(n)
                : new ArrayList<Complex>(Arrays.asList(alpha));
        return talbot.getomega(n, a).toArray(new Complex[0]);
    }

    /**
     * Invert F at t by Talbot's deformed contour. N defaults to 32.
     *
     * <p>Talbot samples F off the real axis, so F must accept a genuinely complex
     * argument; unlike Gaver-Stehfest it cannot be fed a real-only transform.
     */
    public static double laplace_invert_talbot(UnaryOperator<Complex> F, double t, int n) {
        if (!(t > 0.0)) {
            throw new RuntimeException("laplace_invert_talbot: the time point must be positive");
        }
        int nn = (n > 0) ? n : 32;
        Complex[] alpha = talbot_get_alpha(nn);
        Complex[] omega = talbot_get_omega(nn, alpha);
        double res = 0.0;
        for (int i = 0; i < nn; i++) {
            res += omega[i].multiply(F.apply(alpha[i].divide(t))).getReal();
        }
        return res / t;
    }

    /** Invert F at t by the Euler (Abate-Whitt) method. N defaults to 41. */
    public static double laplace_invert_euler(UnaryOperator<Complex> F, double t, int n) {
        return laplace_invert(F, t, "euler", n);
    }

    /** Invert F at t by Gaver-Stehfest. N defaults to 12. */
    public static double laplace_invert_gaver_stehfest(UnaryOperator<Complex> F, double t, int n) {
        return laplace_invert(F, t, "gaver-stehfest", n);
    }

    /** Invert F at t by the Concentrated Matrix Exponential method. N defaults to 25. */
    public static double laplace_invert_cme(UnaryOperator<Complex> F, double t, int n) {
        return laplace_invert(F, t, "cme", n);
    }

    /**
     * Invert F at t by the named method: "euler", "talbot", "gaver-stehfest",
     * "cme" or "weeks".
     *
     * @param n number of terms; 0 takes the method's own default (41 Euler, 32
     *          Talbot, 12 Gaver-Stehfest, 25 CME, 200 the Weeks p0)
     */
    public static double laplace_invert(UnaryOperator<Complex> F, double t, String method, int n) {
        if (!(t > 0.0)) {
            throw new RuntimeException("laplace_invert: the time point must be positive");
        }
        String m = method == null ? "euler" : method.toLowerCase();
        if ("weeks".equals(m) || "laguerre".equals(m)) {
            return laplace_invert_weeks(laplace_weeks_scaling(F, n > 0 ? n : WEEKS_DEFAULT_P0,
                    1e-10), t);
        }
        if ("talbot".equals(m)) {
            return laplace_invert_talbot(F, t, n > 0 ? n : 32);
        }
        double[] T = new double[]{t};
        if ("euler".equals(m)) {
            return iltcme.ilt(F, T, n > 0 ? n : 41, "euler")[0];
        }
        if ("gaver".equals(m) || "gaver-stehfest".equals(m) || "gaver_stehfest".equals(m)) {
            return iltcme.ilt(F, T, n > 0 ? n : 12, "gaver")[0];
        }
        if ("cme".equals(m)) {
            return iltcme.ilt(F, T, n > 0 ? n : 25, "cme")[0];
        }
        throw new RuntimeException("laplace_invert: unknown method '" + method
                + "', expected euler, talbot, gaver-stehfest, cme or weeks");
    }

    /**
     * The DENSITY on a grid: the inversion clamped at zero.
     *
     * <p>A density cannot be negative, and a numerical inversion can undershoot
     * near the origin or in a tail.
     */
    public static double[] laplace_invert_pdf(UnaryOperator<Complex> F, double[] t, String method,
                                              int n) {
        double[] out = new double[t.length];
        String m = method == null ? "euler" : method.toLowerCase();
        if ("weeks".equals(m) || "laguerre".equals(m)) {
            // One expansion serves the whole grid; this is the point of Weeks.
            WeeksParams w = laplace_weeks_scaling(F, n > 0 ? n : WEEKS_DEFAULT_P0, 1e-10);
            for (int i = 0; i < t.length; i++) {
                out[i] = Math.max(0.0, laplace_invert_weeks(w, t[i]));
            }
            return out;
        }
        for (int i = 0; i < t.length; i++) {
            if (!(t[i] > 0.0)) {
                continue;
            }
            out[i] = Math.max(0.0, laplace_invert(F, t[i], method, n));
        }
        return out;
    }

    /**
     * The DISTRIBUTION on a grid, from the transform of the DENSITY.
     *
     * <p>F(s)/s is the transform of the CDF, so that is what is inverted --
     * passing the CDF's own transform here would invert it twice. The result is
     * clamped into [0,1] and made monotone by a running maximum, because a
     * numerical inversion is pointwise and nothing in it enforces either
     * property; a non-monotone "CDF" then yields negative probabilities
     * downstream.
     */
    public static double[] laplace_invert_cdf(UnaryOperator<Complex> F, double[] t, String method,
                                              int n) {
        UnaryOperator<Complex> Fc = new UnaryOperator<Complex>() {
            @Override
            public Complex apply(Complex s) {
                if (s.abs() < 1e-15) {
                    return new Complex(1.0, 0.0);
                }
                return F.apply(s).divide(s);
            }
        };
        double[] out = new double[t.length];
        String m = method == null ? "euler" : method.toLowerCase();
        if ("weeks".equals(m) || "laguerre".equals(m)) {
            WeeksParams w = laplace_weeks_scaling(Fc, n > 0 ? n : WEEKS_DEFAULT_P0, 1e-10);
            for (int i = 0; i < t.length; i++) {
                out[i] = Math.min(1.0, Math.max(0.0, laplace_invert_weeks(w, t[i])));
            }
        } else {
            for (int i = 0; i < t.length; i++) {
                if (!(t[i] > 0.0)) {
                    continue;
                }
                out[i] = Math.min(1.0, Math.max(0.0, laplace_invert(Fc, t[i], method, n)));
            }
        }
        for (int i = 1; i < out.length; i++) {
            out[i] = Math.max(out[i], out[i - 1]);
        }
        return out;
    }
}
