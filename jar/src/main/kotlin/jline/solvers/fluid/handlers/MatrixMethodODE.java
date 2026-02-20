/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */


package jline.solvers.fluid.handlers;

import jline.lang.NetworkStruct;
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
    private Matrix pQa;

    public MatrixMethodODE(
            Matrix W, Matrix SQ, Matrix S, Matrix Qa, Matrix ALambda, int numDimensions,
            boolean[] isSourceState) {
        this.W = W.copy();
        this.SQ = SQ.copy();
        this.S = S.copy();
        this.Qa = Qa.copy();
        this.ALambda = ALambda.copy();
        this.numDimensions = numDimensions;
        this.isSourceState = isSourceState;
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

        this(W, SQ, S, Qa, ALambda, numDimensions, isSourceState);

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
            dxdtTmp = computeDerivativesWithoutSmoothing(xDMS, sumXQa, SQa);
        } else {
            dxdtTmp = computeDerivativesUsingPNormSmoothing(xDMS, sumXQa, SQa);
        }

        for (int i = 0; i < dxdt.length; i++) {
            dxdt[i] = dxdtTmp.get(i);
        }
    }

    private Matrix computeDerivativesUsingPNormSmoothing(
            Matrix x, Matrix sumXQa, Matrix SQa) {

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

        MatrixEquation computeDerivatives = new MatrixEquation();
        computeDerivatives.alias(W, "W", thetaEff, "theta", ALambda, "ALambda");
        computeDerivatives.process("dxdt = W' * theta + ALambda");
        return computeDerivatives.lookupSimple("dxdt");
    }

    private Matrix computeDerivativesWithoutSmoothing(
            Matrix x, Matrix sumXQa, Matrix SQa) {

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

        MatrixEquation computeDerivatives = new MatrixEquation();
        computeDerivatives.alias(W, "W", theta, "theta", ALambda, "ALambda");
        computeDerivatives.process("dxdt = W' * theta + ALambda");
        return computeDerivatives.lookupSimple("dxdt");
    }

    @Override
    public int getDimension() {
        return numDimensions;
    }
}
