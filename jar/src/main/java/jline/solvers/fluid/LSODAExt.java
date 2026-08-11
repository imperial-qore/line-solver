package jline.solvers.fluid;

import odesolver.LSODA;
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

    private final int maxSteps;
    private final double relativeTol;
    private final double absoluteTol;
    private final double hmaxInv;
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
        // hmaxInv = 1/maxStep (inverse of max step size), hmin = minStep
        this.hmaxInv = (maxStep > 0) ? 1.0 / maxStep : 0.0;
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
                0.0, 0.0, this.hmaxInv, this.hminVal);

        // Copy result to output array
        System.arraycopy(this.y, 1, yOut, 0, neq);
        return 0.0;
    }
}
