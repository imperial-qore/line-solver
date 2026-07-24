/**
 * @file Robust Queueing (RQ) approximation for a single G/GI/1 queue
 *
 * Mean steady-state workload of a single G/GI/1 queue partially characterized
 * by its arrival rate, index of dispersion for counts (IDC) and the first two
 * moments of the service time, per W. Whitt and W. You (2018), "A Robust
 * Queueing Network Analyzer Based on Indices of Dispersion", eqs. (13),(16)-(18).
 *
 * @since LINE 3.0
 */
package jline.api.qsys;

import java.util.function.DoubleUnaryOperator;

import org.apache.commons.math3.util.FastMath;

public final class Qsys_gig1_rq {
    private Qsys_gig1_rq() {}

    /**
     * Robust Queueing (RQ) approximation for a single G/GI/1 queue.
     *
     * Implements the mean steady-state workload
     *   Z* = sup_{x&gt;=0} { -(1-rho) x + sqrt( 2 rho x (I_a(x) + c2_s) / mu ) }
     * and the derived steady-state performance measures.
     *
     * @param rho   traffic intensity lambda/mu (0&lt;rho&lt;1)
     * @param mu    service rate
     * @param cs2   service SCV c2_s
     * @param IaFun arrival IDC handle, IaFun(x) -&gt; I_a(x) at time argument x&gt;0
     * @return array {Z, W, Q, X}: mean workload E[Z], waiting time E[W], queue
     *         length E[Q] (waiting + in service), and number in system E[X]
     */
    public static double[] qsys_gig1_rq(double rho, double mu, double cs2, DoubleUnaryOperator IaFun) {
        if (rho <= 0) {
            return new double[]{0.0, 0.0, 0.0, 0.0};
        }
        if (rho >= 1) {
            double inf = Double.POSITIVE_INFINITY;
            return new double[]{inf, inf, inf, inf};
        }
        double lambda = rho * mu;

        final double frho = rho, fmu = mu, fcs2 = cs2;
        final DoubleUnaryOperator fIa = IaFun;
        // objective f(x) = -(1-rho)x + sqrt( 2 rho x (I_a(x)+c2_s)/mu )
        DoubleUnaryOperator objf = new DoubleUnaryOperator() {
            @Override
            public double applyAsDouble(double x) {
                if (x <= 0) {
                    return 0.0;
                }
                double ia = fIa.applyAsDouble(x);
                return -(1.0 - frho) * x + FastMath.sqrt(FastMath.max(0.0, 2.0 * frho * x * (ia + fcs2) / fmu));
            }
        };

        // Coarse log-spaced scan for bracketing followed by golden-section
        // refinement on the best bracket (mirrors the MATLAB scan+fminbnd).
        int n = 200;
        double loExp = -6.0, hiExp = 8.0;
        double[] xs = new double[n];
        double[] fv = new double[n];
        int imax = 0;
        for (int i = 0; i < n; i++) {
            double ex = loExp + (hiExp - loExp) * i / (n - 1);
            xs[i] = FastMath.pow(10.0, ex);
            fv[i] = objf.applyAsDouble(xs[i]);
            if (fv[i] > fv[imax]) {
                imax = i;
            }
        }
        double lo = xs[FastMath.max(0, imax - 1)];
        double hi = xs[FastMath.min(n - 1, imax + 1)];
        double xopt = goldenMax(objf, lo, hi, 1e-10);
        double Z = FastMath.max(fv[imax], objf.applyAsDouble(xopt));
        Z = FastMath.max(Z, 0.0);

        // derived measures (eqs. 16-18)
        double W = FastMath.max(0.0, Z / rho - (cs2 + 1.0) / (2.0 * mu));
        double Q = lambda * W;      // E[Q] waiting (Little's law on waiting time)
        double X = Q + rho;         // E[X] number in system including one in service
        return new double[]{Z, W, Q, X};
    }

    // Golden-section maximization of a unimodal (in practice) function on [a,b].
    private static double goldenMax(DoubleUnaryOperator f, double a, double b, double tol) {
        final double gr = (FastMath.sqrt(5.0) - 1.0) / 2.0;   // 0.618...
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
