/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.qsys;

/**
 * Conditional downstream waiting time in an M/M/1 tandem.
 *
 * <p>The conditional mean waiting time of customer n+1 at the downstream station
 * of a two-station single-server tandem queue, given that customer n waited
 * {@code Wk} upstream and {@code Wk1} downstream.
 *
 * <p>The point of the tandem recursion is that the interarrival time at the
 * downstream station is the interdeparture time upstream, not an independent
 * draw. With {@code A ~ Exp(lambda)} the interarrival time upstream and
 * {@code S1, S1'} the service times upstream of customers n and n+1, that
 * interdeparture time is
 *
 * <pre>
 *   D = max(A - Wk - S1, 0) + S1',
 * </pre>
 *
 * <p>an idle period followed by the next service, and the downstream Lindley step
 * is {@code W2_{n+1} = (Wk1 + S2 - D)^+} with {@code S2 ~ Exp(mu2)} independent of
 * D.
 *
 * <p>Because A is exponential, {@code max(A - Wk - S1, 0)} is zero with
 * probability {@code 1-q} and {@code Exp(lambda)} with probability
 * {@code q = e^{-lambda Wk} mu1/(lambda+mu1)}, the probability the upstream server
 * goes idle, so D is either {@code Exp(mu1)} or the sum of {@code Exp(mu1)} and
 * {@code Exp(lambda)}. Averaging the downstream step over both cases needs only
 * two elementary transforms of
 *
 * <pre>
 *   g(d) = E[(Wk1 + S2 - d)^+] = Wk1 - d + 1/mu2      for d &lt;= Wk1,
 *                              = e^{-mu2 (d-Wk1)}/mu2  for d &gt; Wk1,
 * </pre>
 *
 * <p>namely {@code J(c) = int_0^inf e^{-cu} g(u) du} and
 * {@code Jw(c) = int_0^inf u e^{-cu} g(u) du}, both closed form, giving
 *
 * <pre>
 *   E[W2_{n+1} | Wk, Wk1] = (1-q) mu1 J(mu1) + q C,
 *   C = lambda mu1 (J(mu1) - J(lambda))/(lambda-mu1)  if lambda != mu1,
 *     = mu1^2 Jw(mu1)                                 if lambda == mu1.
 * </pre>
 *
 * <p>As {@code Wk} grows the upstream server never idles, q vanishes, and the mean
 * tends to {@code mu1 J(mu1) = E[g(S1')]}, as it must.
 *
 * <p>Two caveats, both inherited from the reference and both quantified here.
 *
 * <p>First, this is exact for the step taken in isolation, that is when the
 * conditioning pair is independent of the four primitives that drive the step. In
 * a running tandem it is not: the downstream wait {@code Wk1} was itself
 * determined by an interdeparture time containing {@code S1}, so conditioning on
 * {@code (Wk, Wk1)} is not conditioning on a Markov state of the tandem. Measured
 * against a 4e6-customer simulation of the real tandem at lambda = 0.8,
 * mu1 = mu2 = 1, the formula is within 0.4% to 1.3% away from the empty state and
 * 7% at {@code Wk = Wk1 = 0}, where the entanglement is strongest. Treat it as
 * exact for one isolated step and as a good approximation in a running tandem.
 *
 * <p>Second, this closed form was derived rather than transcribed from the
 * reference's theorem 4, because that theorem rests on its proposition 1, which
 * omits a service-time difference and so does not describe a tandem queue; see
 * {@link Qsys_tandem_lindley}. The two differ: at lambda = 0.8, mu1 = mu2 = 1 and
 * {@code Wk = Wk1 = 0} the published route gives 0.3016 against 0.3457 here, the
 * latter matching simulation of the step to 6e-4 relative error.
 *
 * <p>As in the reference, the upstream interarrival time is taken to be
 * {@code Exp(lambda)}, which by Burke's theorem is also the stationary
 * interdeparture law, so the same formula is applied at any pair of consecutive
 * stations of a longer M/M/1 tandem, with the caveat above compounding.
 *
 * <p>Port of MATLAB qsys_mm1_tandem_lindley.m.
 *
 * <p>Reference: S. Palomo, J. Pender, "Learning the Tandem Network Lindley
 * Recursion", Proc. Winter Simulation Conference, 2021, proposition 1 and
 * theorem 4, corrected as described above.
 *
 * @since LINE 3.1.0
 */
public final class Qsys_mm1_tandem_lindley {
    private Qsys_mm1_tandem_lindley() {}

    /**
     * Conditional downstream mean at a single pair of current waiting times.
     *
     * @param lambda external arrival rate upstream, positive
     * @param mu1    upstream service rate, positive
     * @param mu2    downstream service rate, positive
     * @param Wk     current upstream waiting time, finite and nonnegative
     * @param Wk1    current downstream waiting time, finite and nonnegative
     * @return the conditional means
     */
    public static QsysTandemLindleyResult qsys_mm1_tandem_lindley(double lambda, double mu1,
                                                                 double mu2, double Wk,
                                                                 double Wk1) {
        return qsys_mm1_tandem_lindley(lambda, mu1, mu2, new double[] {Wk}, new double[] {Wk1});
    }

    /**
     * Conditional downstream mean at several pairs of current waiting times.
     *
     * @param lambda external arrival rate upstream, positive
     * @param mu1    upstream service rate, positive
     * @param mu2    downstream service rate, positive
     * @param Wk     current upstream waiting times, finite and nonnegative
     * @param Wk1    current downstream waiting times, same length as Wk
     * @return the conditional means
     */
    public static QsysTandemLindleyResult qsys_mm1_tandem_lindley(double lambda, double mu1,
                                                                 double mu2, double[] Wk,
                                                                 double[] Wk1) {
        if (!(lambda > 0.0) || !Double.isFinite(lambda)) {
            throw new IllegalArgumentException("lambda=" + lambda + " must be positive and finite");
        }
        if (!(mu1 > 0.0) || !Double.isFinite(mu1)) {
            throw new IllegalArgumentException("mu1=" + mu1 + " must be positive and finite");
        }
        if (!(mu2 > 0.0) || !Double.isFinite(mu2)) {
            throw new IllegalArgumentException("mu2=" + mu2 + " must be positive and finite");
        }
        if (Wk == null || Wk1 == null || Wk.length == 0 || Wk.length != Wk1.length) {
            throw new IllegalArgumentException(
                    "Wk and Wk1 must be nonempty and of equal length");
        }
        for (int i = 0; i < Wk.length; i++) {
            if (!Double.isFinite(Wk[i]) || Wk[i] < 0.0
                    || !Double.isFinite(Wk1[i]) || Wk1[i] < 0.0) {
                throw new IllegalArgumentException(
                        "Wk and Wk1 must hold finite nonnegative values");
            }
        }

        int nw = Wk.length;
        double[] mean = new double[nw];
        double[] interdep = new double[nw];
        double[] idle = new double[nw];
        boolean equalRates = Math.abs(lambda - mu1) <= 1e-9 * Math.max(lambda, mu1);

        for (int i = 0; i < nw; i++) {
            double x = Wk[i];
            double y = Wk1[i];
            double q = Math.exp(-lambda * x) * mu1 / (lambda + mu1);
            double base = mu1 * j(mu1, y, mu2);
            double conv;
            if (!equalRates) {
                conv = lambda * mu1 / (lambda - mu1) * (j(mu1, y, mu2) - j(lambda, y, mu2));
            } else {
                // the two rates coincide, the interdeparture time is Erlang(2,mu1)
                conv = mu1 * mu1 * jw(mu1, y, mu2);
            }
            mean[i] = (1.0 - q) * base + q * conv;
            interdep[i] = 1.0 / mu1 + q / lambda;
            idle[i] = q;
        }

        return new QsysTandemLindleyResult(mean, interdep, idle);
    }

    /** J(c) = int_0^inf e^{-c u} E[(y + S2 - u)^+] du with S2 ~ Exp(mu2). */
    private static double j(double c, double y, double mu2) {
        double e = Math.exp(-c * y);
        return (y + 1.0 / mu2) * (1.0 - e) / c
                - (1.0 - e * (1.0 + c * y)) / (c * c)
                + e / (mu2 * (c + mu2));
    }

    /** Jw(c) = int_0^inf u e^{-c u} E[(y + S2 - u)^+] du with S2 ~ Exp(mu2). */
    private static double jw(double c, double y, double mu2) {
        double e = Math.exp(-c * y);
        double d = c + mu2;
        return (y + 1.0 / mu2) * (1.0 - e * (1.0 + c * y)) / (c * c)
                - (2.0 - e * (2.0 + 2.0 * c * y + c * c * y * y)) / (c * c * c)
                + e * (y / d + 1.0 / (d * d)) / mu2;
    }
}
