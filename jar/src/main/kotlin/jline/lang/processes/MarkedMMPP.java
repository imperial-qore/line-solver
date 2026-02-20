/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import jline.util.Maths;
import jline.util.RandomManager;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import org.apache.commons.math3.util.FastMath;

import java.io.Serializable;
import java.util.Random;

import static jline.api.mam.Map_momentKt.map_moment;
import static jline.api.mam.Mmap_timereverseKt.mmap_timereverse;
import static jline.io.InputOutputKt.line_error;
import static jline.io.InputOutputKt.mfilename;

/**
 * A Marked Markov-Modulated Poisson Process (M3PP)
 */
public class MarkedMMPP extends Marked implements Serializable {

    public MarkedMMPP() {
        super("MarkedMMPP", 0);
    }

    public MarkedMMPP(MatrixCell mmap) {
        super("MMAP", mmap.size());
        setProcess(mmap);
        Matrix D0 = mmap.get(0);
        Matrix D1 = mmap.get(1);
        int nPhases = D0.getNumCols();
        this.setParam(1, "D0", D0);
        this.setParam(2, "D1", D1);
        if (!D1.isDiag()) {
            line_error(mfilename(new Object() {
            }), "Non-diagonal D1 matrix in Marked MMPP");
        }
        for (int i = 2; i < mmap.size(); i++) {
            this.setParam(i + 1, "D1" + i, mmap.get(i));
            if (!mmap.get(i).isDiag()) {
                line_error(mfilename(new Object() {
                }), "Non-diagonal D1K matrix in Marked MMPP");
            }
        }

        this.nPhases = nPhases;
    }

    /**
     * Evaluates the cumulative distribution function (CDF) at time t.
     *
     * @param t the time value
     * @return the CDF value at time t
     */
    public double evalCDF(double t) {
        if (t <= 0) {
            return 0.0;
        }

        try {
            // Use aggregated MAP approach for CDF calculation
            Matrix D0 = D(0);
            MatrixCell proc = getProcess();

            // Create aggregated D1
            Matrix D1_agg = new Matrix(D0.getNumRows(), D0.getNumCols());
            for (int k = 1; k < proc.size(); k++) {
                D1_agg.addEq(1.0, proc.get(k));
            }

            // CDF calculation using MAP formula
            Matrix pie = jline.api.mam.Map_pieKt.map_pie(D0, D1_agg);
            Matrix D0t = Matrix.scaleMult(D0, t);
            Matrix expT = Maths.matrixExp(D0t);
            Matrix ones = Matrix.ones(D0.getNumRows(), 1);
            Matrix result = pie.mult(expT).mult(ones);
            return 1.0 - result.get(0, 0);
        } catch (Exception e) {
            // Fallback: exponential approximation
            double mean = getMean();
            return 1.0 - Math.exp(-t / mean);
        }
    }

    /**
     * Evaluates the Laplace-Stieltjes Transform (LST) at parameter s.
     *
     * @param s the transform parameter
     * @return the LST value at parameter s
     */
    public double evalLST(double s) {
        if (s < 0) {
            throw new IllegalArgumentException("LST parameter s must be non-negative");
        }

        try {
            Matrix D0 = D(0);
            MatrixCell proc = getProcess();

            // Create aggregated D1
            Matrix D1_agg = new Matrix(D0.getNumRows(), D0.getNumCols());
            for (int k = 1; k < proc.size(); k++) {
                D1_agg.addEq(1.0, proc.get(k));
            }

            // LST calculation using MAP formula
            Matrix pie = jline.api.mam.Map_pieKt.map_pie(D0, D1_agg);
            Matrix I = Matrix.eye(D0.getNumRows());
            Matrix negT = Matrix.scaleMult(D0, -1.0);
            Matrix sI = Matrix.scaleMult(I, s);
            Matrix resolvent = negT.add(1.0, sI).inv();
            Matrix ones = Matrix.ones(D1_agg.getNumCols(), 1);
            Matrix t = D1_agg.mult(ones);
            Matrix result = pie.mult(resolvent).mult(t);
            return result.get(0, 0);
        } catch (Exception e) {
            // Fallback: exponential approximation
            double mean = getMean();
            return 1.0 / (1.0 + s * mean);
        }
    }

    /**
     * Returns the process representation as a MatrixCell.
     *
     * @return the process matrices
     */
    @Override
    public MatrixCell getProcess() {
        return this.process;
    }

    public double getRate() {
        return 1.0 / this.getMean();
    }

    public double getSCV() {
        double mean = this.getMean();
        return this.getVar() / mean / mean;
    }

    public double getSkewness() {
        double E1 = map_moment(D(0), D(1), 1);
        double E2 = map_moment(D(0), D(1), 2);
        double E3 = map_moment(D(0), D(1), 3);
        double skew = E3 - 3 * E2 * E1 + 2 * E1 * E1 * E1;
        double scv = (E2 - E1 * E1) / E1 / E1;
        skew = skew / FastMath.pow(Math.sqrt(scv) * E1, 3);
        return skew;
    }

    public double getVar() {
        double E1 = map_moment(D(0), D(1), 1);
        double E2 = map_moment(D(0), D(1), 2);
        return E2 - E1 * E1;
    }

    /**
     * Normalizes the MarkedMMPP so that D0+sum(D_i) rows form a proper infinitesimal generator.
     * For MMPP, all D_i matrices should be diagonal with non-negative elements.
     */
    public void normalize() {
        Matrix D0 = D(0).copy();
        MatrixCell proc = getProcess();

        // Ensure non-negative off-diagonal elements in D0 and diagonal elements in D_i
        for (int i = 0; i < nPhases; i++) {
            for (int j = 0; j < nPhases; j++) {
                if (i != j && D0.get(i, j) < 0) {
                    D0.set(i, j, 0.0);
                }
                // For MMPP, D_k matrices should be diagonal with non-negative elements
                for (int k = 1; k < proc.size(); k++) {
                    Matrix Dk = proc.get(k);
                    if (i == j && Dk.get(i, j) < 0) {
                        Dk.set(i, j, 0.0); // Diagonal elements should be non-negative
                    } else if (i != j) {
                        Dk.set(i, j, 0.0); // Off-diagonal elements should be zero for MMPP
                    }
                }
            }
        }

        // Adjust diagonal elements of D0 to make row sums zero
        for (int i = 0; i < nPhases; i++) {
            double rowSum = 0.0;
            // Sum off-diagonal elements of D0
            for (int j = 0; j < nPhases; j++) {
                if (i != j) {
                    rowSum += D0.get(i, j);
                }
            }
            // Sum diagonal elements in D_k matrices (k > 0) for row i
            for (int k = 1; k < proc.size(); k++) {
                Matrix Dk = proc.get(k);
                rowSum += Dk.get(i, i); // Only diagonal elements for MMPP
            }
            // Set diagonal element to make row sum zero
            D0.set(i, i, -rowSum);
        }

        // Update D0 in the process
        proc.set(0, D0);
        setProcess(proc);
    }

    @Override
    public double[] sample(int n) {
        return this.sample(n, RandomManager.getThreadRandomAsRandom());
    }

    /**
     * Generates samples from the MarkedMMPP using the specified random generator.
     *
     * @param n      the number of samples to generate
     * @param random the random number generator to use
     * @return array of inter-arrival times
     */
    @Override
    public double[] sample(int n, Random random) {
        try {
            // For MMPP, we can use MAP sampling with aggregated diagonal matrices
            Matrix D0 = D(0);
            MatrixCell proc = getProcess();

            // Create aggregated D1 by summing all diagonal arrival matrices
            Matrix D1_agg = new Matrix(D0.getNumRows(), D0.getNumCols());
            for (int k = 1; k < proc.size(); k++) {
                D1_agg.addEq(1.0, proc.get(k));
            }

            // Use MAP sampling with aggregated arrival matrix
            return jline.api.mam.Map_sampleKt.map_sample(D0, D1_agg, n, random);
        } catch (Exception e) {
            // Fallback: generate exponential samples
            double[] samples = new double[n];
            double mean = getMean();
            for (int i = 0; i < n; i++) {
                samples[i] = -mean * Math.log(random.nextDouble());
            }
            return samples;
        }
    }

    public MarkedMMPP toTimeReversed() {
        MatrixCell reversed = mmap_timereverse(this.process);
        return new MarkedMMPP(reversed);
    }
}
