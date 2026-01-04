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
    private Matrix pQa;

    public MatrixMethodODE(
            Matrix W, Matrix SQ, Matrix S, Matrix Qa, Matrix ALambda, int numDimensions) {
        this.W = W.copy();
        this.SQ = SQ.copy();
        this.S = S.copy();
        this.Qa = Qa.copy();
        this.ALambda = ALambda.copy();
        this.numDimensions = numDimensions;
        this.pQa = new Matrix(0, 0);
    }

    public MatrixMethodODE(
            Matrix W,
            Matrix SQ,
            Matrix S,
            Matrix Qa,
            Matrix ALambda,
            int numDimensions,
            NetworkStruct sn,
            List<Double> pStarValues) {

        this(W, SQ, S, Qa, ALambda, numDimensions);

        this.pQa = new Matrix(SQ.getNumRows(), 1);
        int row = 0;
        for (int i = 0; i < sn.nstations; i++) {
            double pStarValue = pStarValues.get(i);
            for (int j = 0; j < sn.nclasses; j++) {
                int nPhases = (int) sn.phases.get(i, j);
                for (int k = 0; k < nPhases; k++) {
                    pQa.set(row, 0, pStarValue);
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
        if (this.pQa.getNumRows() == 0) { // If no pStar values have been specified
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

        // ghat = smoothed processor-share constraint approximation, per Ruuskanen et. al
        Matrix ghat = Matrix.createLike(new Matrix(x));
        for (int i = 0; i < x.getNumRows(); i++) {
            // x, c and p as per Ruuskanen's Julia implementation
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

        MatrixEquation computeDerivatives = new MatrixEquation();
        computeDerivatives.alias(W, "W", x, "x", ghat, "ghat", ALambda, "ALambda");
        computeDerivatives.process("dxdt = W' * (x .* ghat) + ALambda");
        return computeDerivatives.lookupSimple("dxdt");
    }

    private Matrix computeDerivativesWithoutSmoothing(
            Matrix x, Matrix sumXQa, Matrix SQa) {

        int sumXQaRows = sumXQa.getNumRows();
        int sumXQaCols = sumXQa.getNumCols();
        Matrix minOfSumXQaAndSQa = new Matrix(sumXQaRows, sumXQaCols);
        for (int i = 0; i < sumXQaRows; i++) {
            for (int j = 0; j < sumXQaCols; j++) {
                minOfSumXQaAndSQa.set(i, j, min(sumXQa.get(i, j), SQa.get(i, j)));
            }
        }

        MatrixEquation computeDerivatives = new MatrixEquation();
        computeDerivatives.alias(
                W, "W", x, "x", sumXQa, "sumXQa", minOfSumXQaAndSQa, "minOfSumXQaAndSQa", ALambda, "ALambda");
        computeDerivatives.process("dxdt = W' * (x ./ sumXQa .* minOfSumXQaAndSQa) + ALambda");
        return computeDerivatives.lookupSimple("dxdt");
    }

    @Override
    public int getDimension() {
        return numDimensions;
    }
}
