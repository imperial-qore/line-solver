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

        this.pQa = new Matrix(SQ.getNumRows(), 1);
        int row = 0;
        for (int i = 0; i < sn.nstations; i++) {
            double pStarValue = pStarValues.get(i);
            for (int j = 0; j < sn.nclasses; j++) {
                int nPhases = (int) sn.phases.get(i, j);
                if (nPhases == 0) nPhases = 1;
                for (int k = 0; k < nPhases; k++) {
                    if (row < pQa.getNumRows()) {
                        pQa.set(row, 0, pStarValue);
                    }
                    row++;
                }
            }
        }
    }

    @Override
    public void computeDerivatives(double t, double[] x, double[] dxdt)
            throws MaxCountExceededException, DimensionMismatchException {

        Matrix xDMS = new Matrix(x.length, 1);
        for (int i = 0; i < x.length; i++) {
            xDMS.set(i, 0, x[i]);
        }

        MatrixEquation calculateSumXQa = new MatrixEquation();
        calculateSumXQa.alias(xDMS, "x", SQ, "SQ", GlobalConstants.FineTol, "distribZero");
        calculateSumXQa.process("sumXQa = distribZero + SQ * x");
        Matrix sumXQa = calculateSumXQa.lookupSimple("sumXQa");

        int QaCols = this.Qa.getNumCols();
        Matrix SQa = new Matrix(QaCols, 1);
        for (int i = 0; i < QaCols; i++) {
            SQa.set(i, 0, S.get((int) Qa.get(0, i), 0));
        }

        Matrix dxdtTmp;
        if (this.pQa.getNumRows() == 0) {
            dxdtTmp = computeDerivativesWithoutSmoothing(xDMS, sumXQa, SQa, t);
        } else {
            dxdtTmp = computeDerivativesUsingPNormSmoothing(xDMS, sumXQa, SQa, t);
        }

        for (int i = 0; i < dxdt.length; i++) {
            dxdt[i] = dxdtTmp.get(i);
        }
    }

    private Matrix computeDerivativesUsingPNormSmoothing(
            Matrix x, Matrix sumXQa, Matrix SQa, double t) {

        Matrix ghat = Matrix.createLike(new Matrix(x));
        for (int i = 0; i < x.getNumRows(); i++) {
            double xVal = sumXQa.get(i, 0);
            double cVal = SQa.get(i, 0);
            double pVal = pQa.get(i, 0);
            double ghatVal = 1.0 / FastMath.pow(1 + FastMath.pow(xVal / cVal, pVal), 1.0 / pVal);
            if (Double.isNaN(ghatVal)) {
                ghat.set(i, 0, 0);
            } else {
                ghat.set(i, 0, ghatVal);
            }
        }

        Matrix thetaEff = new Matrix(x.getNumRows(), 1);
        for (int i = 0; i < x.getNumRows(); i++) {
            if (isSourceState != null && i < isSourceState.length && isSourceState[i]) {
                thetaEff.set(i, 0, 0.0);
            } else {
                thetaEff.set(i, 0, x.get(i, 0) * ghat.get(i, 0));
            }
        }

        applyStateMultiplier(thetaEff, t);

        MatrixEquation computeDerivatives = new MatrixEquation();
        computeDerivatives.alias(W, "W", thetaEff, "theta", ALambda, "ALambda");
        computeDerivatives.process("dxdt = W' * theta + ALambda");
        return computeDerivatives.lookupSimple("dxdt");
    }

    private Matrix computeDerivativesWithoutSmoothing(
            Matrix x, Matrix sumXQa, Matrix SQa, double t) {

        int nStates = x.getNumRows();
        Matrix theta = new Matrix(nStates, 1);
        for (int i = 0; i < nStates; i++) {
            if (isSourceState != null && i < isSourceState.length && isSourceState[i]) {
                theta.set(i, 0, 0.0);
            } else {
                double xVal = x.get(i, 0);
                double sumVal = sumXQa.get(i, 0);
                double sVal = SQa.get(i, 0);
                theta.set(i, 0, xVal / sumVal * min(sumVal, sVal));
            }
        }

        applyStateMultiplier(theta, t);

        MatrixEquation computeDerivatives = new MatrixEquation();
        computeDerivatives.alias(W, "W", theta, "theta", ALambda, "ALambda");
        computeDerivatives.process("dxdt = W' * theta + ALambda");
        return computeDerivatives.lookupSimple("dxdt");
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
