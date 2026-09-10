package jline.lib.m3a;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import de.xypron.jcobyla.Calcfc;
import de.xypron.jcobyla.Cobyla;

import jline.api.mam.Map_acf;
import jline.api.mam.Map_idc;
import jline.api.mam.Map_moment;
import jline.api.mam.Map_scv;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * Utility functions for M3A (Markovian Arrival Process with 3-moment Approximation) compression.
 */
public final class M3aUtils {
    private M3aUtils() {}

    /**
     * Computes the autocorrelation function of an MMAP up to the specified lag.
     */
    public static Double[] computeAutocorrelation(MatrixCell MMAP, int maxLag) {
        Matrix D0 = MMAP.get(0);
        Matrix D1 = MMAP.get(1);

        Matrix lags = new Matrix(1, maxLag, maxLag);
        for (int i = 0; i < maxLag; i++) {
            lags.set(0, i, (double) (i + 1));
        }

        Matrix acfMatrix = Map_acf.map_acf(D0, D1, lags);
        Double[] result = new Double[maxLag];
        for (int i = 0; i < maxLag; i++) {
            result[i] = Double.valueOf(acfMatrix.get(i));
        }
        return result;
    }

    /**
     * Computes the index of dispersion for counts (IDC) of an MMAP.
     */
    public static double computeIDC(MatrixCell MMAP, double timeWindow) {
        Matrix D0 = MMAP.get(0);
        Matrix D1 = MMAP.get(1);
        return Map_idc.map_idc(D0, D1);
    }

    /**
     * Computes the coefficient of variation of an MMAP.
     */
    public static double computeCoeffVar(MatrixCell MMAP) {
        Matrix D0 = MMAP.get(0);
        Matrix D1 = MMAP.get(1);
        return Math.sqrt(Map_scv.map_scv(D0, D1));
    }

    /**
     * Computes the first n moments of an MMAP.
     */
    public static Double[] computeMoments(MatrixCell MMAP, int n) {
        Matrix D0 = MMAP.get(0);
        Matrix D1 = MMAP.get(1);

        Double[] moments = new Double[n];
        for (int k = 1; k <= n; k++) {
            moments[k - 1] = Double.valueOf(Map_moment.map_moment(D0, D1, k));
        }
        return moments;
    }

    /**
     * Computes the spectral gap of an MMAP generator matrix.
     */
    public static double computeSpectralGap(MatrixCell MMAP) {
        Matrix Q = MMAP.get(0).add(1.0, MMAP.get(1));
        Matrix eigenvalues = Q.eigval().values;

        List<Double> eigenList = new ArrayList<Double>();
        for (int i = 0; i < eigenvalues.getNumRows(); i++) {
            eigenList.add(Double.valueOf(eigenvalues.get(i, 0)));
        }
        Collections.sort(eigenList, Collections.reverseOrder());

        if (eigenList.size() >= 2) {
            return eigenList.get(0) - eigenList.get(1);
        }
        return 0.0;
    }

    public static Double[] optimizeParameters(Double[] initialParams,
                                              final ObjectiveFunction objectiveFunction,
                                              final ConstraintFunction[] constraints) {
        return optimizeParameters(initialParams, objectiveFunction, constraints, 1e-6);
    }

    /**
     * Optimizes MMAP parameters using COBYLA optimization.
     */
    public static Double[] optimizeParameters(final Double[] initialParams,
                                              final ObjectiveFunction objectiveFunction,
                                              final ConstraintFunction[] constraints,
                                              double tolerance) {
        final double[] x = new double[initialParams.length];
        for (int i = 0; i < initialParams.length; i++) {
            x[i] = initialParams[i];
        }

        Calcfc calcfc = new Calcfc() {
            public double compute(int n, int m, double[] xVars, double[] con) {
                Double[] xArray = new Double[n];
                for (int i = 0; i < n; i++) xArray[i] = Double.valueOf(xVars[i]);

                for (int i = 0; i < constraints.length; i++) {
                    con[i] = constraints[i].evaluate(xArray);
                }
                return objectiveFunction.evaluate(xArray);
            }
        };

        Cobyla.findMinimum(calcfc, initialParams.length, constraints.length,
                x, 0.5, tolerance, 0, 1000);

        Double[] result = new Double[x.length];
        for (int i = 0; i < x.length; i++) {
            result[i] = Double.valueOf(x[i]);
        }
        return result;
    }

    public static double computeKLDivergence(MatrixCell MMAP1, MatrixCell MMAP2) {
        return computeKLDivergence(MMAP1, MMAP2, 10000, 23000);
    }

    public static double computeKLDivergence(MatrixCell MMAP1, MatrixCell MMAP2, int numSamples) {
        return computeKLDivergence(MMAP1, MMAP2, numSamples, 23000);
    }

    /**
     * Computes the Kullback-Leibler divergence between two MMAPs.
     */
    public static double computeKLDivergence(MatrixCell MMAP1, MatrixCell MMAP2,
                                             int numSamples, int seed) {
        Double[] samples1 = sampleInterArrivalTimes(MMAP1, numSamples, seed);
        Double[] samples2 = sampleInterArrivalTimes(MMAP2, numSamples, seed);

        double[] hist1 = computeHistogram(samples1, 50);
        double[] hist2 = computeHistogram(samples2, 50);

        double kl = 0.0;
        for (int i = 0; i < hist1.length; i++) {
            if (hist1[i] > 0 && hist2[i] > 0) {
                kl += hist1[i] * Math.log(hist1[i] / hist2[i]);
            }
        }
        return kl;
    }

    /**
     * Validates that a matrix represents a valid MMAP.
     */
    public static boolean validateMMAP(MatrixCell MMAP) {
        if (MMAP.size() < 2) return false;

        Matrix D0 = MMAP.get(0);
        Matrix D1 = MMAP.get(1);

        if (D0.getNumRows() != D0.getNumCols() || D1.getNumRows() != D1.getNumCols()) return false;
        if (D0.getNumRows() != D1.getNumRows()) return false;

        for (int i = 0; i < D0.getNumRows(); i++) {
            if (D0.get(i, i) > 0) return false;
        }

        for (int i = 0; i < D0.getNumRows(); i++) {
            for (int j = 0; j < D0.getNumCols(); j++) {
                if (i != j && D0.get(i, j) < 0) return false;
            }
        }

        for (int i = 0; i < D1.getNumRows(); i++) {
            for (int j = 0; j < D1.getNumCols(); j++) {
                if (D1.get(i, j) < 0) return false;
            }
        }

        Matrix Q = D0.add(1.0, D1);
        Matrix e = Matrix.ones(Q.getNumRows(), 1);
        Matrix rowSums = Q.mult(e);

        for (int i = 0; i < rowSums.getNumRows(); i++) {
            if (Math.abs(rowSums.get(i, 0)) > 1e-10) return false;
        }
        return true;
    }

    private static Double[] sampleInterArrivalTimes(MatrixCell MMAP, int numSamples, int seed) {
        Double[] samples = new Double[numSamples];

        java.util.Random random = new java.util.Random(seed);
        int state = 0;
        Matrix Q = MMAP.get(0).add(1.0, MMAP.get(1));
        double[] rates = new double[Q.getNumRows()];
        for (int i = 0; i < Q.getNumRows(); i++) {
            rates[i] = -Q.get(i, i);
        }

        for (int i = 0; i < numSamples; i++) {
            double u = random.nextDouble();
            samples[i] = Double.valueOf(-Math.log(u) / rates[state]);
            state = (state + 1) % Q.getNumRows();
        }
        return samples;
    }

    private static double[] computeHistogram(Double[] samples, int numBins) {
        double[] hist = new double[numBins];
        double minVal = samples[0];
        double maxVal = samples[0];
        for (Double s : samples) {
            if (s < minVal) minVal = s;
            if (s > maxVal) maxVal = s;
        }
        double binWidth = (maxVal - minVal) / numBins;

        for (Double sample : samples) {
            int binIndex = Math.min(numBins - 1, (int) ((sample - minVal) / binWidth));
            hist[binIndex] += 1.0;
        }

        double total = 0.0;
        for (double v : hist) total += v;
        if (total > 0) {
            for (int i = 0; i < hist.length; i++) {
                hist[i] /= total;
            }
        }
        return hist;
    }

    /**
     * Functional interface for objective functions used by COBYLA.
     */
    public interface ObjectiveFunction {
        double evaluate(Double[] x);
    }

    /**
     * Functional interface for constraint functions used by COBYLA.
     */
    public interface ConstraintFunction {
        double evaluate(Double[] x);
    }
}
