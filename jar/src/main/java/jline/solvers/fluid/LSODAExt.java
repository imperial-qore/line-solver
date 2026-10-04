package jline.solvers.fluid;

import odesolver.LSODA;

import static jline.io.InputOutput.line_warning;
import static jline.io.InputOutput.mfilename;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

/**
 * Extended LSODA solver that supports configurable maximum internal steps.
 *
 * The upstream lsoda-java library hardcodes mxstep=1000000 inside the lsoda()
 * method, ignoring any reflection-based field setting. This subclass overrides
 * integrate() to call the public lsoda() method directly with the desired
 * mxstep value.
 */
public class LSODAExt extends LSODA {

    /**
     * Measurement scaffold for the currently-invisible failure class: an lsoda()
     * call that sets a negative istate and returns NORMALLY, rather than throwing.
     * Such a call has its (un-integrated) state copied out as if it had converged.
     * <p>
     * This counts and reports rather than throwing, ON PURPOSE. It is not yet
     * established that this failure mode is reachable at all: the failures we can
     * see today arrive as RuntimeExceptions and are handled by the retry in
     * {@link jline.solvers.fluid.analyzers.ClosingAndStateDepMethodsAnalyzer}.
     * If the count is zero, every reachable failure goes through the exception
     * path and a throw here is free. If it is not zero, some body of currently
     * "passing" results is quietly wrong, and WHICH models those are is a larger
     * question than the missing check -- so that must be measured and reported
     * before a throw converts them into failures.
     */
    private static final java.util.concurrent.atomic.AtomicLong INTEGRATE_CALLS =
            new java.util.concurrent.atomic.AtomicLong();
    private static final java.util.concurrent.atomic.AtomicLong NEGATIVE_ISTATE_CALLS =
            new java.util.concurrent.atomic.AtomicLong();
    /** First-seen istate code, so the report names a code without flooding stdout. */
    private static final java.util.concurrent.atomic.AtomicLong FIRST_NEGATIVE_ISTATE =
            new java.util.concurrent.atomic.AtomicLong(0);

    /** Total integrate() calls since the last reset. */
    public static long getIntegrateCallCount() {
        return INTEGRATE_CALLS.get();
    }

    /** Calls that returned normally with a negative istate, i.e. silent failures. */
    public static long getNegativeIstateCount() {
        return NEGATIVE_ISTATE_CALLS.get();
    }

    /** The first negative istate code observed, or 0 if none. */
    public static long getFirstNegativeIstate() {
        return FIRST_NEGATIVE_ISTATE.get();
    }

    public static void resetIstateCounters() {
        INTEGRATE_CALLS.set(0);
        NEGATIVE_ISTATE_CALLS.set(0);
        FIRST_NEGATIVE_ISTATE.set(0);
    }

    private final int maxSteps;
    private final double relativeTol;
    private final double absoluteTol;
    private final double hmaxVal;
    private final double hminVal;

    /**
     * Create an LSODAExt solver with configurable max steps.
     *
     * @param minStep   minimum step size (absolute value)
     * @param maxStep   maximum step size (absolute value)
     * @param rtol      relative tolerance
     * @param atol      absolute tolerance
     * @param meth      method flag (12 = nonstiff/stiff auto)
     * @param miter     iteration method (5 = diagonal Jacobian)
     * @param maxSteps  maximum number of internal steps (replaces hardcoded 1000000)
     */
    public LSODAExt(double minStep, double maxStep, double rtol, double atol,
                    int meth, int miter, int maxSteps) {
        super(minStep, maxStep, rtol, atol, meth, miter);
        this.maxSteps = maxSteps;
        this.relativeTol = rtol;
        this.absoluteTol = atol;
        // lsoda() takes the max STEP and inverts it itself; 0 means unbounded
        this.hmaxVal = (maxStep > 0 && !Double.isInfinite(maxStep)) ? maxStep : 0.0;
        this.hminVal = minStep;
    }

    @Override
    public double integrate(FirstOrderDifferentialEquations equations, double t0,
                            double[] y0, double t, double[] yOut) {
        this.ode = equations;
        int neq = equations.getDimension();

        double[] rtolArr = new double[]{0, this.relativeTol};
        double[] atolArr = new double[]{0, this.absoluteTol};

        // see _kb/06-solver-catalog.md (JAR-only implementation notes: LSODAExt mxstep reflection hack)
        this.lsoda(neq, y0, t0, t, 1, rtolArr, atolArr, 1, 1, 1,
                0, this.maxSteps, 0, 0, 0,
                0.0, 0.0, this.hmaxVal, this.hminVal);

        // lsoda() is VOID: it publishes its status on the inherited public field
        // istate, not as a return value. This read is the whole point of the
        // measurement -- a negative istate returned NORMALLY (rather than thrown)
        // means the state copied out below never integrated, and no caller can
        // currently tell that from a converged solve.
        INTEGRATE_CALLS.incrementAndGet();
        if (this.istate < 0) {
            long seen = NEGATIVE_ISTATE_CALLS.incrementAndGet();
            FIRST_NEGATIVE_ISTATE.compareAndSet(0, this.istate);
            if (seen == 1) {
                // Reported once only: this sits inside the layer ODE loop and a
                // per-call message would bury the run. The counters carry the rest.
                line_warning(mfilename(new Object(){}),
                        "lsoda returned NORMALLY with istate=" + this.istate
                        + " over t in [" + t0 + ", " + t + "], neq=" + neq
                        + ": the state copied out did not integrate. Further"
                        + " occurrences are counted, not printed"
                        + " (LSODAExt.getNegativeIstateCount()).");
            }
        }

        // Copy result to output array
        System.arraycopy(this.y, 1, yOut, 0, neq);
        return 0.0;
    }
}
