/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import jline.api.mam.Dmap_pie;
import jline.lib.butools.dmap.RandomDMAP;
import jline.lib.butools.dmap.SamplesFromDMAP;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import java.io.Serializable;
import java.util.Random;

/**
 * A Discrete Markovian Arrival Process (DMAP).
 * Unlike MAP (continuous-time), D0+D1 is a stochastic matrix (row sums = 1).
 */
public class DMAP extends MarkovModulated implements Serializable {

    public DMAP(MatrixCell map) {
        this(map.get(0), map.get(1));
    }

    public DMAP() {
        super("DMAP", 2);
    }

    public DMAP(Matrix D0, Matrix D1) {
        super("DMAP", 2);
        this.nPhases = D0.getNumCols();
        this.setParam(1, "D0", D0);
        this.setParam(2, "D1", D1);
        MatrixCell rep = new MatrixCell();
        rep.set(0, D0);
        rep.set(1, D1);
        this.setProcess(rep);
    }

    public static DMAP rand() {
        return DMAP.rand(2);
    }

    public static DMAP rand(int order) {
        jline.util.Pair<Matrix, Matrix> pair = RandomDMAP.randomDMAP(order, 10.0, 0, 1000, 1e-7, new Random());
        return new DMAP(pair.getFirst(), pair.getSecond());
    }

    @Override
    public Matrix D(int i) {
        return this.process.get(i);
    }

    private Matrix getStationaryVector() {
        return Dmap_pie.dmap_pie(D(0), D(1));
    }

    /**
     * Mean inter-arrival time: pi * (I - D0)^{-1} * e
     */
    @Override
    public double getMean() {
        int n = D(0).getNumRows();
        Matrix pi = getStationaryVector();
        Matrix I = Matrix.eye(n);
        Matrix ImD0inv = I.sub(D(0)).inv();
        Matrix ones = Matrix.ones(n, 1);
        return pi.mult(ImD0inv).mult(ones).get(0, 0);
    }

    /**
     * Variance: 2*pi*(I-D0)^{-2}*e - mean - mean^2
     */
    @Override
    public double getVar() {
        int n = D(0).getNumRows();
        Matrix pi = getStationaryVector();
        Matrix I = Matrix.eye(n);
        Matrix ImD0inv = I.sub(D(0)).inv();
        Matrix ones = Matrix.ones(n, 1);
        double mean = getMean();
        double term = pi.mult(ImD0inv.mult(ImD0inv)).mult(ones).get(0, 0);
        return 2.0 * term - mean - mean * mean;
    }

    @Override
    public double getSCV() {
        double mean = getMean();
        return getVar() / (mean * mean);
    }

    /**
     * Row sums of (I - D0).
     */
    @Override
    public Matrix getMu() {
        Matrix D0 = D(0);
        int n = D0.getNumRows();
        Matrix res = new Matrix(n, 1, n);
        for (int i = 0; i < n; i++) {
            double rowSum = 0.0;
            for (int j = 0; j < n; j++) {
                if (i == j) {
                    rowSum += 1.0 - D0.get(i, j);
                } else {
                    rowSum -= D0.get(i, j);
                }
            }
            res.set(i, 0, rowSum);
        }
        return res;
    }

    @Override
    public double[] sample(int n, Random random) {
        Random rng = (random != null) ? random : new Random();
        int[] intSamples = SamplesFromDMAP.samplesFromDMAP(D(0), D(1), n, null, 1e-14, rng);
        double[] result = new double[intSamples.length];
        for (int i = 0; i < intSamples.length; i++) {
            result[i] = intSamples[i];
        }
        return result;
    }

    @Override
    public double[] sample(int n) {
        return sample(n, new Random());
    }

    @Override
    public double getRate() {
        return 1.0 / getMean();
    }

    /**
     * CDF(t) = 1 - pi * D0^floor(t) * e
     */
    @Override
    public double evalCDF(double t) {
        if (t < 1.0) return 0.0;
        int k = (int) Math.floor(t);
        int n = D(0).getNumRows();
        Matrix pi = getStationaryVector();
        Matrix D0k = Matrix.eye(n);
        Matrix D0 = D(0);
        for (int i = 0; i < k; i++) {
            D0k = D0k.mult(D0);
        }
        Matrix ones = Matrix.ones(n, 1);
        return 1.0 - pi.mult(D0k).mult(ones).get(0, 0);
    }

    /**
     * Embedded chain P = (I - D0)^{-1} * D1.
     */
    public Matrix getTransitionMatrix() {
        int n = D(0).getNumRows();
        Matrix I = Matrix.eye(n);
        return I.sub(D(0)).inv().mult(D(1));
    }

    @Override
    public long getNumberOfPhases() {
        return ((Matrix) this.getParam(1).getValue()).getNumCols();
    }

    public void normalize() {
        Matrix D0 = D(0).copy();
        Matrix D1 = D(1).copy();
        int n = D0.getNumRows();
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (D0.get(i, j) < 0) D0.set(i, j, 0.0);
                if (D1.get(i, j) < 0) D1.set(i, j, 0.0);
            }
        }
        for (int i = 0; i < n; i++) {
            double rowSum = 0.0;
            for (int j = 0; j < n; j++) {
                rowSum += D0.get(i, j) + D1.get(i, j);
            }
            if (rowSum > 0) {
                for (int j = 0; j < n; j++) {
                    D0.set(i, j, D0.get(i, j) / rowSum);
                    D1.set(i, j, D1.get(i, j) / rowSum);
                }
            }
        }
        this.setParam(1, "D0", D0);
        this.setParam(2, "D1", D1);
        MatrixCell res = new MatrixCell();
        res.set(0, D0);
        res.set(1, D1);
        this.setProcess(res);
    }

    @Override
    public MatrixCell getProcess() {
        MatrixCell res = new MatrixCell();
        res.set(0, this.D(0));
        res.set(1, this.D(1));
        return res;
    }

    @Override
    public String toString() {
        Matrix D0 = D(0);
        return String.format("DMAP(%dx%d, mean=%.6f, scv=%.6f)",
                D0.getNumRows(), D0.getNumCols(), getMean(), getSCV());
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (obj == null || getClass() != obj.getClass()) return false;
        DMAP other = (DMAP) obj;
        Matrix D0 = D(0), D1 = D(1);
        Matrix oD0 = other.D(0), oD1 = other.D(1);
        if (D0.getNumRows() != oD0.getNumRows()) return false;
        for (int i = 0; i < D0.getNumRows(); i++)
            for (int j = 0; j < D0.getNumCols(); j++) {
                if (Math.abs(D0.get(i, j) - oD0.get(i, j)) > 1e-10) return false;
                if (Math.abs(D1.get(i, j) - oD1.get(i, j)) > 1e-10) return false;
            }
        return true;
    }

    @Override
    public int hashCode() {
        Matrix D0 = D(0), D1 = D(1);
        int result = 1;
        for (int i = 0; i < D0.getNumRows(); i++)
            for (int j = 0; j < D0.getNumCols(); j++) {
                result = 31 * result + Double.hashCode(D0.get(i, j));
                result = 31 * result + Double.hashCode(D1.get(i, j));
            }
        return result;
    }
}
