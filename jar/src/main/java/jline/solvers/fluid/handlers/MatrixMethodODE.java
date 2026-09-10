/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */


package jline.solvers.fluid.handlers;

import jline.lang.NetworkStruct;
import jline.solvers.fluid.analyzers.FluidStateRateMultiplier;
import jline.GlobalConstants;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixEquation;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.util.FastMath;

import java.util.List;

import static org.apache.commons.math3.util.FastMath.min;

public class MatrixMethodODE implements FirstOrderDifferentialEquations {

    private final Matrix W;
    private final Matrix SQ;
    private final Matrix S;
    private final Matrix Qa;
    private final Matrix ALambda;
    private final int numDimensions;
    private final boolean[] isSourceState;
    private final FluidStateRateMultiplier stateMult;
    private Matrix pQa;
    private boolean[] isInfState;
    /**
     * When set, min(E[n], c) is replaced by E[min(n, c)] under the station's
     * equilibrium geometric marginal. Set ONLY by the degeneracy repair in
     * MatrixMethodAnalyzer: min() is FLAT above the server count, so a network of
     * saturated stations has a CONTINUUM of fixed points and the integrator
     * returns whichever one it stopped at. See BUGS.md and _kb/06-solver-catalog.md.
     */
    private boolean varClosure = false;

    public MatrixMethodODE(
            Matrix W, Matrix SQ, Matrix S, Matrix Qa, Matrix ALambda, int numDimensions,
            boolean[] isSourceState) {
        this(W, SQ, S, Qa, ALambda, numDimensions, isSourceState, (FluidStateRateMultiplier) null);
    }

    /**
     * @param stateMult time-varying per-state rate multiplier, or null for the
     *                  autonomous (constant-rate) ODE
     */
    public MatrixMethodODE(
            Matrix W, Matrix SQ, Matrix S, Matrix Qa, Matrix ALambda, int numDimensions,
            boolean[] isSourceState, FluidStateRateMultiplier stateMult) {
        this.W = W.copy();
        this.SQ = SQ.copy();
        this.S = S.copy();
        this.Qa = Qa.copy();
        this.ALambda = ALambda.copy();
        this.numDimensions = numDimensions;
        this.isSourceState = isSourceState;
        this.stateMult = stateMult;
        this.pQa = new Matrix(0, 0);
        this.isInfState = new boolean[numDimensions];
    }

    public MatrixMethodODE(
            Matrix W,
            Matrix SQ,
            Matrix S,
            Matrix Qa,
            Matrix ALambda,
            int numDimensions,
            boolean[] isSourceState,
            NetworkStruct sn,
            List<Double> pStarValues) {
        this(W, SQ, S, Qa, ALambda, numDimensions, isSourceState, sn, pStarValues, null);
    }

    /**
     * @param stateMult time-varying per-state rate multiplier, or null for the
     *                  autonomous (constant-rate) ODE
     */
    public MatrixMethodODE(
            Matrix W,
            Matrix SQ,
            Matrix S,
            Matrix Qa,
            Matrix ALambda,
            int numDimensions,
            boolean[] isSourceState,
            NetworkStruct sn,
            List<Double> pStarValues,
            FluidStateRateMultiplier stateMult) {

        this(W, SQ, S, Qa, ALambda, numDimensions, isSourceState, stateMult);

        // Read the exponent off Qa, which is already filtered by keep and by the
        // immediate-state elimination; walking the FULL (station,class,phase)
        // space instead misaligns the exponents whenever a state was dropped
        int nStates = this.Qa.getNumCols();
        this.pQa = new Matrix(nStates, 1);
        this.isInfState = new boolean[nStates];
        for (int i = 0; i < nStates; i++) {
            int ist = (int) this.Qa.get(0, i);
            int idx = pStarValues.size() == 1 ? 0 : FastMath.min(ist, pStarValues.size() - 1);
            pQa.set(i, 0, pStarValues.get(idx));
            // An INF station has k = infinity, so its share is 1 with no min()
            // to smooth; S holds the population there, which the p-norm would
            // otherwise read as a k = N queue
            this.isInfState[i] = Double.isInfinite(sn.nservers.get(ist, 0));
        }
    }

    @Override
    public void computeDerivatives(double t, double[] x, double[] dxdt)
            throws MaxCountExceededException, DimensionMismatchException {

        Matrix xDMS = new Matrix(x.length, 1);
        for (int i = 0; i < x.length; i++) {
            xDMS.set(i, 0, x[i]);
        }

        Matrix theta = theta(xDMS);
        applyStateMultiplier(theta, t);

        MatrixEquation computeDerivatives = new MatrixEquation();
        computeDerivatives.alias(W, "W", theta, "theta", ALambda, "ALambda");
        computeDerivatives.process("dxdt = W' * theta + ALambda");
        Matrix dxdtTmp = computeDerivatives.lookupSimple("dxdt");

        for (int i = 0; i < dxdt.length; i++) {
            dxdt[i] = dxdtTmp.get(i);
        }
    }

    /**
     * The mass in service, i.e. theta(x) of Ruuskanen et al., PEVA 151 (2021):
     * eq. (12) without smoothing, eq. (26)-(27) with an exponent set.
     *
     * <p>The analyzer reads the metrics off this same theta: eq. (23) takes the
     * utilization from the share the ODE integrated, so U and T must not revert
     * to the hard min() after a smoothed solve. The state multiplier is NOT
     * applied here, since the throughput scales by it while the utilization
     * does not.</p>
     */
    /**
     * Turns on the variance-carrying saturation term and declares which states
     * belong to INF stations, which carry no min() to close at all.
     */
    public void setVarClosure(boolean[] infStates) {
        this.varClosure = true;
        if (infStates != null) {
            this.isInfState = infStates.clone();
        }
    }

    public Matrix theta(Matrix x) {
        MatrixEquation calculateSumXQa = new MatrixEquation();
        calculateSumXQa.alias(x, "x", SQ, "SQ", GlobalConstants.FineTol, "distribZero");
        calculateSumXQa.process("sumXQa = distribZero + SQ * x");
        Matrix sumXQa = calculateSumXQa.lookupSimple("sumXQa");

        int nStates = x.getNumRows();
        Matrix SQa = new Matrix(nStates, 1);
        for (int i = 0; i < nStates; i++) {
            SQa.set(i, 0, S.get((int) Qa.get(0, i), 0));
        }

        boolean smoothed = this.pQa.getNumRows() > 0;
        Matrix theta = new Matrix(nStates, 1);
        for (int i = 0; i < nStates; i++) {
            if (isSourceState != null && i < isSourceState.length && isSourceState[i]) {
                theta.set(i, 0, 0.0);
                continue;
            }
            double xVal = x.get(i, 0);
            double sumVal = sumXQa.get(i, 0);
            double sVal = SQa.get(i, 0);
            if (smoothed && !(isInfState != null && i < isInfState.length && isInfState[i])) {
                double pVal = pQa.get(i, 0);
                double ghatVal = 1.0 / FastMath.pow(1 + FastMath.pow(sumVal / sVal, pVal), 1.0 / pVal);
                theta.set(i, 0, Double.isNaN(ghatVal) ? 0.0 : xVal * ghatVal);
            } else if (smoothed) {
                theta.set(i, 0, xVal);
            } else if (varClosure) {
                // E[min(n, c)] UNDER A GEOMETRIC MARGINAL, not min(E[n], c):
                //   n ~ Geometric(mean m) => E[min(n,c)] = sum_{k=1..c} p^k
                //                          = m * (1 - p^c),  p = m/(1+m).
                // Strictly increasing in m (slope 1/(1+m)^2 at c = 1, still 1e-2
                // at m = 9, a restoring force the integrator can follow inside
                // its horizon) and with the same asymptotes, -> c as m -> inf and
                // -> m as m -> 0. An INF station has a server per job, so there
                // is no min() to close and S holds the whole population there.
                double eMin;
                if (isInfState != null && i < isInfState.length && isInfState[i]) {
                    eMin = sumVal;
                } else {
                    double p = sumVal > 0 ? sumVal / (1.0 + sumVal) : 0.0;
                    eMin = sumVal * (1.0 - FastMath.pow(p, FastMath.max(sVal, 0.0)));
                    if (Double.isNaN(eMin) || Double.isInfinite(eMin)) {
                        eMin = 0.0;
                    }
                    eMin = min(eMin, sumVal);
                }
                theta.set(i, 0, xVal / sumVal * eMin);
            } else {
                theta.set(i, 0, xVal / sumVal * min(sumVal, sVal));
            }
        }
        return theta;
    }

    /**
     * Scales the effective service population in place by the time-varying
     * per-state rate multiplier. The drift out of a state is linear in that
     * (station,class) service rate, so scaling theta scales every rate out of
     * the state by exactly m(t), which is what a time-varying rate means.
     */
    private void applyStateMultiplier(Matrix theta, double t) {
        if (stateMult == null) {
            return;
        }
        double[] m = stateMult.multAt(t);
        int n = Math.min(theta.getNumRows(), m.length);
        for (int i = 0; i < n; i++) {
            theta.set(i, 0, theta.get(i, 0) * m[i]);
        }
    }

    @Override
    public int getDimension() {
        return numDimensions;
    }
}
