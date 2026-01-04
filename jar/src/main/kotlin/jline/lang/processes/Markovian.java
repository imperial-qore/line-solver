/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import static jline.GlobalConstants.Inf;

import jline.api.mam.*;
import jline.GlobalConstants;
import jline.util.Pair;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import org.apache.commons.math3.util.FastMath;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * An abstract class for a Markovian distribution
 */
public class Markovian extends ContinuousDistribution implements Serializable {

    protected MatrixCell process; // Associated renewal process: <0, D0>, <1, D1>, ...
    protected int nPhases;

    /**
     * Creates a new Markovian distribution with the specified name and parameter count.
     * 
     * @param name the name of this distribution type
     * @param numParam the number of parameters
     */
    public Markovian(String name, int numParam) {
        super(name, numParam, new Pair<Double, Double>(0.0, Inf));
        nPhases = 0;
    }

    /**
     * Gets the i-th matrix of the Markovian arrival process representation.
     * 
     * @param i the matrix index (0 for D0, 1 for D1, etc.)
     * @return the matrix at index i
     */
    public Matrix D(int i) {
        return this.process.get(i);
    }

    /**
     * Evaluates the cumulative distribution function at the given point.
     * 
     * @param t the point at which to evaluate the CDF
     * @return the CDF value at point t
     */
    public double evalCDF(double t) {
        Matrix D0 = D(0);
        Matrix D1 = D(1);
        return Map_cdfKt.map_cdf(D0, D1, Matrix.singleton(t)).toDouble();
    }

    /**
     * Evaluates the cumulative distribution function at multiple points.
     * 
     * @param tset array of points at which to evaluate the CDF
     * @return array of CDF values
     */
    public double[] evalCDF(double[] tset) {
        Matrix D0 = D(0);
        Matrix D1 = D(1);
        return Map_cdfKt.map_cdf(D0, D1, new Matrix(tset)).toArray1D();
    }

    @Override
    public double evalLST(double s) {
        Matrix pie = Map_pieKt.map_pie(D(0), D(1));
        Matrix A = D(0);
        Matrix e = Matrix.ones(nPhases, 1);
        Matrix sI = Matrix.eye(nPhases).scale(s);
        return pie.mult((sI.sub(A)).inv()).mult(A.scale(-1)).mult(e).value();
    }

    /**
     * Evaluates the mean count at time t.
     * 
     * @param t the time point
     * @return the mean count at time t
     */
    public double evalMeanT(double t) {
        return t / getMean();
    }

    /**
     * Evaluates the variance count at time t.
     * 
     * @param t the time point
     * @return the variance count at time t
     */
    public double evalVarT(double t) {
        Matrix D0 = D(0);
        Matrix D1 = D(1);
        return Map_varcountKt.map_varcount(D0, D1, t);
    }

    /**
     * Gets the autocorrelation function at the specified lags.
     * 
     * @param lags matrix of lag values
     * @return matrix of autocorrelation values
     */
    public Matrix getACF(Matrix lags) {
        return Map_acfKt.map_acf(D(0), D(1), lags);
    }

    /**
     * Gets the embedded Markov chain transition matrix.
     * 
     * @return the embedded chain matrix
     */
    public Matrix getEmbedded() {
        return Map_embeddedKt.map_embedded(D(0), D(1));
    }

    /**
     * Gets the stationary probability vector of the embedded Markov chain.
     * 
     * @return the embedded chain stationary probabilities
     */
    public Matrix getEmbeddedProb() {
        return Map_pieKt.map_pie(D(0), D(1));
    }

    /**
     * Gets the index of dispersion for counts (IDC).
     * 
     * @return the IDC value
     */
    public double getIDC() {
        return Map_idcKt.map_idc(D(0), D(1));
    }

    /**
     * Gets the index of dispersion for intervals (IDI).
     * For renewal processes, IDI = IDC.
     * 
     * @return the IDI value
     */
    public double getIDI() {
        return getIDC();
    }

    /**
     * Gets the initial probability vector.
     * 
     * @return the initial probabilities
     */
    public Matrix getInitProb() {
        return getEmbeddedProb();
    }

    /**
     * Gets the mean of this Markovian distribution.
     * 
     * @return the mean value, or NaN if the process contains NaN values
     */
    public double getMean() {
        MatrixCell rep = getProcess();
        if (rep.get(0).hasNaN()) {
            return Double.NaN;
        } else {
            Matrix D0 = D(0);
            Matrix D1 = D(1);
            return Map_meanKt.map_mean(D0, D1);
        }
    }

    /**
     * Sets the mean of this Markovian distribution by scaling the process.
     * 
     * @param newMean the new mean value
     */
    public void setMean(double newMean) {
        Matrix D0 = D(0);
        Matrix D1 = D(1);
        if (newMean == 0.0) {
            this.process = Map_scaleKt.map_scale(D0, D1, GlobalConstants.Zero);
        } else {
            this.process = Map_scaleKt.map_scale(D0, D1, newMean);
        }
    }

    /**
     * Gets the first three moments of this distribution.
     * 
     * @return list containing the first, second, and third moments
     */
    public List<Double> getMoments() {
        List<Double> moments = new ArrayList<>();
        for (int i = 1; i <= 3; i++) {
            moments.add(Map_momentKt.map_moment(D(0), D(1), i));
        }
        return moments;
    }

    /**
     * Gets the diagonal rate matrix containing the negative diagonal elements of D0.
     * 
     * @return column vector of rates
     */
    public Matrix getMu() {
        Matrix aph_1 = D(0);
        int size = FastMath.min(aph_1.getNumCols(), aph_1.getNumRows());
        Matrix res = new Matrix(size, 1, size);
        for (int i = 0; i < size; i++) {
            res.set(i, 0, -aph_1.get(i, i));
        }
        return res;
//    	Matrix aph_1 = D(1);
//    	int size = FastMath.min(aph_1.getNumCols(), aph_1.getNumRows());
//    	Matrix res = new Matrix(size, 1, size);
//    	for(int i = 0; i < size; i++) {
//			for(int j = 0; j < size; j++) {
//				res.set(i, 0, res.get(i)+aph_1.get(i, j));
//			}
//    	}
//    	return res;
    }

    /**
     * Gets the number of phases in this Markovian distribution.
     * 
     * @return the number of phases
     */
    public long getNumberOfPhases() {
        MatrixCell PH = this.getProcess();
        return PH.get(1).getNumCols();
    }

    /**
     * Gets the exit probability vector (phi).
     * 
     * @return column vector of exit probabilities for each phase
     */
    public Matrix getPhi() {
        MatrixCell aph = this.getProcess();
        Matrix ones = new Matrix(aph.get(0).getNumRows(), 1, aph.get(0).getNumRows());
        Matrix res = new Matrix(aph.get(1).getNumRows(), 1);
        Matrix mu = getMu();

        ones.fill(1.0);
        aph.get(1).mult(ones, res);
        res.divideRows(mu.getNonZeroValues(), 0);
        return res;
    }

    /**
     * Gets the matrix representation of this Markovian process.
     * 
     * @return MatrixCell containing D0, D1, ... matrices
     */
    public MatrixCell getProcess() {
        return process;
    }

    /**
     * Sets the matrix representation of this Markovian process.
     * 
     * @param D MatrixCell containing D0, D1, ... matrices
     */
    public void setProcess(MatrixCell D) {
        process = D;
    }

    @Override
    public double getRate() {
        return 1.0 / getMean();
    }

    /**
     * Sets the rate of this Markovian distribution by scaling the process.
     * 
     * @param newRate the new rate value
     */
    public void setRate(double newRate) {
        Matrix D0 = D(0);
        Matrix D1 = D(1);
        if (newRate == 0.0) {
            this.process = Map_scaleKt.map_scale(D0, D1, 1 / GlobalConstants.Zero);
        } else {
            this.process = Map_scaleKt.map_scale(D0, D1, 1 / newRate);
        }
    }

    public double getSCV() {
        MatrixCell rep = getProcess();
        if (rep.get(0).hasNaN()) {
            return Double.NaN;
        } else {
            Matrix D0 = D(0);
            Matrix D1 = D(1);
            return Map_scvKt.map_scv(D0, D1);
        }
    }

    public double getSkewness() {
        return Map_skewKt.map_skew(D(0), D(1));
    }

    public Matrix getSubgenerator() {
        return D(0).add(D(1));
    }

    @Override
    public double getVar() {
        return Map_varKt.map_var(D(0), D(1));
    }

    public double getVariance() {
        MatrixCell rep = getProcess();
        if (rep.get(0).hasNaN()) {
            return Double.NaN;
        } else {
            Matrix D0 = D(0);
            Matrix D1 = D(1);
            return Map_varKt.map_var(D0, D1);
        }
    }

    public double[] sample(int n) {
        return this.sample(n, null);
    }

    public double[] sample(int n, Random random) {
        return Map_sampleKt.map_sample(D(0), D(1), n, random);
    }

    // =================== KOTLIN-STYLE PROPERTY ALIASES ===================
    
    /**
     * Kotlin-style property alias for getMean()
     */
    public double mean() {
        return getMean();
    }
    
    /**
     * Kotlin-style property alias for getRate()
     */
    public double rate() {
        return getRate();
    }
    
    /**
     * Kotlin-style property alias for getSCV()
     */
    public double scv() {
        return getSCV();
    }
    
    /**
     * Kotlin-style property alias for getSkewness()
     */
    public double skewness() {
        return getSkewness();
    }
    
    /**
     * Kotlin-style property alias for getVar()
     */
    public double var() {
        return getVar();
    }
    
    /**
     * Kotlin-style property alias for getVariance()
     */
    public double variance() {
        return getVariance();
    }
    
    /**
     * Kotlin-style property alias for getACF()
     */
    public Matrix acf(int maxLag) {
        // Create a matrix with lags from 0 to maxLag
        int[] lags = new int[maxLag + 1];
        for (int i = 0; i <= maxLag; i++) {
            lags[i] = i;
        }
        return getACF(new Matrix(lags));
    }
    
    /**
     * Kotlin-style property alias for getEmbedded()
     */
    public Matrix embedded() {
        return getEmbedded();
    }
    
    /**
     * Kotlin-style property alias for getEmbeddedProb()
     */
    public Matrix embeddedProb() {
        return getEmbeddedProb();
    }
    
    /**
     * Kotlin-style property alias for getIDC()
     */
    public double idc() {
        return getIDC();
    }
    
    /**
     * Kotlin-style property alias for getIDI()
     */
    public double idi() {
        return getIDI();
    }
    
    /**
     * Kotlin-style property alias for getInitProb()
     */
    public Matrix initProb() {
        return getInitProb();
    }
    
    /**
     * Kotlin-style property alias for getMoments()
     */
    public List<Double> moments() {
        return getMoments();
    }
    
    /**
     * Kotlin-style property alias for getMu()
     */
    public Matrix mu() {
        return getMu();
    }
    
    /**
     * Kotlin-style property alias for getNumberOfPhases()
     */
    public long numberOfPhases() {
        return getNumberOfPhases();
    }
    
    /**
     * Kotlin-style property alias for getNumberOfPhases()
     */
    public long numPhases() {
        return getNumberOfPhases();
    }
    
    /**
     * Kotlin-style property alias for getPhi()
     */
    public Matrix phi() {
        return getPhi();
    }
    
    /**
     * Kotlin-style property alias for getProcess()
     */
    public MatrixCell process() {
        return getProcess();
    }
    
    /**
     * Kotlin-style property alias for getSubgenerator()
     */
    public Matrix subgenerator() {
        return getSubgenerator();
    }

}
