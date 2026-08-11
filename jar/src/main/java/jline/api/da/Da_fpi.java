/**
 * @file Generic decomposition-aggregation fixed-point driver
 *
 * @since LINE 3.0
 */
package jline.api.da;

/**
 * Generic damped successive-substitution driver for decomposition-aggregation
 * (DA) fixed-point iterations. Each call to the sweep performs one DA pass:
 * solve the isolated submodels given the current coupling iterate x, exchange
 * flows or rates, and return the updated iterate together with the comparison
 * baseline for the convergence test.
 *
 * Mirrors matlab/src/api/da/da_fpi.m.
 */
public final class Da_fpi {
    private Da_fpi() {}

    /**
     * One DA sweep from iterate x at sweep count it. Returns the updated
     * iterate xnew and the baseline xref for the convergence test; return
     * xref = x for a standard successive-substitution test, or a mid-sweep
     * checkpoint when the method compares against a renormalized iterate.
     */
    public interface Sweep<T> {
        SweepResult<T> sweep(T x, int it);
    }

    /** Pair returned by a sweep: updated iterate and comparison baseline. */
    public static final class SweepResult<T> {
        public final T xnew;
        public final T xref;

        public SweepResult(T xnew, T xref) {
            this.xnew = xnew;
            this.xref = xref;
        }
    }

    /** Convergence measure between the updated iterate and the baseline. */
    public interface Norm<T> {
        double eval(T xnew, T xref);
    }

    /** Final iterate, sweeps executed, and convergence flag. */
    public static final class Result<T> {
        public final T x;
        public final int it;
        public final boolean converged;

        Result(T x, int it, boolean converged) {
            this.x = x;
            this.it = it;
            this.converged = converged;
        }
    }

    /** Driver options. nanstop replicates legacy while-loop drivers whose
     * "continue while delta &gt; tol" test exits on a NaN measure; when false
     * (default) a NaN measure keeps iterating, as in legacy "break if delta
     * &lt; tol" drivers. Convergence is not tested before miniter sweeps. */
    public static final class Options<T> {
        public final int iterMax;
        public final double iterTol;
        public final Norm<T> norm;
        public boolean nanstop = false;
        public int miniter = 1;

        public Options(int iterMax, double iterTol, Norm<T> norm) {
            this.iterMax = iterMax;
            this.iterTol = iterTol;
            this.norm = norm;
        }
    }

    /** Drive the fixed point by successive substitution. */
    public static <T> Result<T> run(Sweep<T> sweep, T x0, Options<T> options) {
        T x = x0;
        int it = 0;
        boolean converged = false;
        for (it = 1; it <= options.iterMax; it++) {
            SweepResult<T> r = sweep.sweep(x, it);
            double delta = options.norm.eval(r.xnew, r.xref);
            x = r.xnew;
            if (it >= options.miniter) {
                if (delta < options.iterTol) {
                    converged = true;
                    break;
                } else if (options.nanstop && Double.isNaN(delta)) {
                    break;
                }
            }
        }
        if (it > options.iterMax) {
            it = options.iterMax;
        }
        return new Result<T>(x, it, converged);
    }

    /** Max absolute elementwise difference norm for double[] iterates. */
    public static Norm<double[]> maxAbsDiff() {
        return new Norm<double[]>() {
            @Override
            public double eval(double[] xnew, double[] xref) {
                double d = 0.0;
                boolean seen = false;
                for (int i = 0; i < xnew.length; i++) {
                    double a = Math.abs(xnew[i] - xref[i]);
                    if (!Double.isNaN(a)) {
                        if (!seen || a > d) {
                            d = a;
                            seen = true;
                        }
                    }
                }
                return seen ? d : Double.NaN;
            }
        };
    }
}
