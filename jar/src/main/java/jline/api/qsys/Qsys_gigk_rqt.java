/**
 * @file Robust Queueing Theory (RQT) worst-case system time of a G/G/k queue
 *
 * Worst-case analysis of an FCFS queue whose arrival and service processes are
 * described by polyhedral uncertainty sets rather than by distributions, per
 * C. Bandi, D. Bertsimas, N. Youssef (2015), "Robust Queueing Theory",
 * Operations Research 63(3), 676-700.
 *
 * @since LINE 3.0
 */
package jline.api.qsys;

import java.util.function.DoubleUnaryOperator;

import org.apache.commons.math3.util.FastMath;

public final class Qsys_gigk_rqt {
    private Qsys_gigk_rqt() {}

    /**
     * Robust Queueing Theory worst-case system time of a G/G/k FCFS queue.
     *
     * The arrival and service processes are not described by distributions but
     * by the polyhedral uncertainty sets
     * U^a = { T : (sum_{i=k+1}^n T_i - (n-k)/lambda)/(n-k)^(1/alpha_a) &gt;= -Gamma_a },
     * U^s = { X : (sum_{i=k}^n X_i - (n-k+1)/mu)/(n-k+1)^(1/alpha_s) &lt;= Gamma_s },
     * whose shape follows the (generalized) central limit theorem: alpha=2 is the
     * finite-variance regime and alpha in (1,2) the heavy-tailed one.
     *
     * <p>The first returned value is the closed-form bound of Theorem 3 (Theorem 8
     * when the two tail coefficients differ, with alphabar = min(alpha_a,alpha_s)),
     * W &lt;= (ab-1)/ab^(ab/(ab-1)) lambda^(1/(ab-1)) (Gamma_a+Gamma_s/k^(1/ab))^(ab/(ab-1))
     * / (1-rho)^(1/(ab-1)) + k/lambda, which for k=1 reduces to Theorem 2 and, at
     * alphabar=2, to the Kingman-like form (lambda/4)(Gamma_a+Gamma_s)^2/(1-rho)
     * + 1/lambda. The third is the exact worst case over the uncertainty sets,
     * eq. (45), the supremum over the integer x &gt;= 1 of
     * x/mu + Gamma_s x^(1/alpha_s) - k(x-1)/lambda + Gamma_a (k(x-1))^(1/alpha_a).
     * The arrival deviation ADDS to the worst case, since the adversary shortens
     * the interarrival times.
     *
     * <p>W is a SYSTEM time (waiting plus service), and its additive term is
     * k/lambda rather than the mean service time 1/mu.
     *
     * @param lambda  arrival rate
     * @param mu      service rate of each server
     * @param Gamma_a variability parameter of the arrival uncertainty set
     * @param Gamma_s variability parameter of the service uncertainty set
     * @param k       number of servers
     * @param alpha_a arrival tail coefficient in (1,2]
     * @param alpha_s service tail coefficient in (1,2]
     * @return array {W, rhohat, Sworst}: closed-form bound on the system time,
     *         modified utilization, and exact worst-case system time
     */
    public static double[] qsys_gigk_rqt(double lambda, double mu, double Gamma_a, double Gamma_s,
                                         int k, double alpha_a, double alpha_s) {
        if (alpha_a <= 1 || alpha_a > 2 || alpha_s <= 1 || alpha_s > 2) {
            throw new RuntimeException("RQT tail coefficients must lie in (1,2].");
        }
        double rho = lambda / (k * mu);
        if (lambda <= 0) {
            return new double[]{1.0 / mu, 0.0, 1.0 / mu};
        }
        if (rho >= 1) {
            double inf = Double.POSITIVE_INFINITY;
            return new double[]{inf, 1.0, inf};
        }

        // Theorem 8 collapses to Theorem 3 when the two tails agree
        double ab = FastMath.min(alpha_a, alpha_s);
        double beta = Gamma_a + Gamma_s / FastMath.pow(k, 1.0 / ab);
        double W;
        if (beta <= 0) {
            // a nonpositive effective variability leaves only the deterministic term
            W = k / lambda;
        } else {
            W = (ab - 1) / FastMath.pow(ab, ab / (ab - 1)) * FastMath.pow(lambda, 1.0 / (ab - 1))
                    * FastMath.pow(beta, ab / (ab - 1)) / FastMath.pow(1 - rho, 1.0 / (ab - 1))
                    + k / lambda;
        }
        double rhohat = W * lambda / (1 + W * lambda);

        final double fmu = mu, fga = Gamma_a, fgs = Gamma_s, faa = alpha_a, fas = alpha_s, flam = lambda;
        final int fk = k;
        DoubleUnaryOperator objf = new DoubleUnaryOperator() {
            @Override
            public double applyAsDouble(double x) {
                double y = FastMath.max(x - 1.0, 0.0);
                return x / fmu + fgs * FastMath.pow(x, 1.0 / fas) - fk * y / flam
                        + fga * FastMath.pow(fk * y, 1.0 / faa);
            }
        };
        // the continuous maximizer of the bounding problem, eq. (16), sizes the scan
        double xstar = beta > 0 ? FastMath.pow(lambda * beta / (ab * (1 - rho)), ab / (ab - 1)) : 1.0;
        double xhi = FastMath.max(4.0, FastMath.ceil(4.0 * xstar));
        int n = 400;
        double loExp = 0.0, hiExp = FastMath.log10(xhi);
        double[] xs = new double[n];
        int m = 0;
        for (int i = 0; i < n; i++) {
            double x = FastMath.max(1.0, FastMath.round(FastMath.pow(10.0, loExp + (hiExp - loExp) * i / (n - 1.0))));
            if (m == 0 || x != xs[m - 1]) {
                xs[m++] = x;   // unique(round(logspace(...))), as in MATLAB
            }
        }
        double Sworst = Double.NEGATIVE_INFINITY;
        int imax = 0;
        for (int i = 0; i < m; i++) {
            double v = objf.applyAsDouble(xs[i]);
            if (v > Sworst) {
                Sworst = v;
                imax = i;
            }
        }
        // refine on the continuous relaxation, then round back onto the lattice
        double xlo = xs[FastMath.max(0, imax - 1)];
        double xup = xs[FastMath.min(m - 1, imax + 1)];
        if (xup > xlo) {
            double xc = goldenMax(objf, xlo, xup, 1e-8);
            double[] cand = {FastMath.floor(xc), FastMath.ceil(xc)};
            for (int i = 0; i < cand.length; i++) {
                if (cand[i] >= 1.0) {
                    Sworst = FastMath.max(Sworst, objf.applyAsDouble(cand[i]));
                }
            }
        }
        return new double[]{W, rhohat, Sworst};
    }

    // Golden-section maximization of a unimodal (in practice) function on [a,b].
    private static double goldenMax(DoubleUnaryOperator f, double a, double b, double tol) {
        final double gr = (FastMath.sqrt(5.0) - 1.0) / 2.0;
        double c = b - gr * (b - a);
        double d = a + gr * (b - a);
        double fc = f.applyAsDouble(c);
        double fd = f.applyAsDouble(d);
        int maxit = 500;
        for (int it = 0; it < maxit && FastMath.abs(b - a) > tol; it++) {
            if (fc > fd) {
                b = d;
                d = c;
                fd = fc;
                c = b - gr * (b - a);
                fc = f.applyAsDouble(c);
            } else {
                a = c;
                c = d;
                fc = fd;
                d = a + gr * (b - a);
                fd = f.applyAsDouble(d);
            }
        }
        return (a + b) / 2.0;
    }
}
