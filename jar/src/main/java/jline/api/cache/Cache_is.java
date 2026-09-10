/**
 * @file Cache Analysis via Importance Sampling
 *
 * @since LINE 3.0
 */
package jline.api.cache;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.util.FastMath;

import jline.io.InputOutput;
import jline.io.Ret;
import jline.util.matrix.Matrix;

public final class Cache_is {
    private Cache_is() {}

    /**
     * Estimate the normalizing constant of the cache steady state distribution
     * using Monte Carlo importance sampling.
     */
    public static Ret.cacheIs cache_is(Matrix gamma, Matrix m, int samples) {
        return cache_is(gamma, m, samples, null, null);
    }

    /**
     * Estimate the (cost-capped) normalizing constant with the feasibility
     * indicator I{S in O} of Casale-Gast, IEEE/ACM Trans. Networking 29(2),
     * 2021, Sec. IX-B.
     *
     * @param gamma Cache access factors (n x h).
     * @param m Cache capacity vector (1 x h).
     * @param samples Number of Monte Carlo samples.
     * @param sigma Item storage costs (sizes); null or empty for none.
     * @param k Per-list storage cost caps; null or empty for none.
     * @return the estimate and its logarithm.
     */
    public static Ret.cacheIs cache_is(Matrix gamma, Matrix m, int samples, Matrix sigma, Matrix k) {
        boolean capped = sigma != null && k != null && !sigma.isEmpty() && !k.isEmpty();
        // Remove items with zero gamma
        boolean[] rowsToKeep = new boolean[gamma.getNumRows()];
        boolean[] colsToKeep = new boolean[gamma.getNumCols()];
        for (int i = 0; i < gamma.getNumCols(); i++) {
            colsToKeep[i] = true;
        }
        int validRows = 0;
        for (int i = 0; i < gamma.getNumRows(); i++) {
            if (gamma.getRow(i).elementSum() > 0) {
                rowsToKeep[i] = true;
                validRows++;
            }
        }

        Matrix filteredGamma = (validRows < gamma.getNumRows()) ? gamma.getSlice(rowsToKeep, colsToKeep) : gamma;
        Matrix filteredSigma = sigma;
        if (capped && validRows < gamma.getNumRows()) {
            filteredSigma = new Matrix(1, validRows);
            int c = 0;
            for (int i = 0; i < gamma.getNumRows(); i++) {
                if (rowsToKeep[i]) {
                    filteredSigma.set(0, c, sigma.get(i));
                    c++;
                }
            }
        }

        int n = filteredGamma.getNumRows();
        int h = filteredGamma.getNumCols();
        int mt = (int) m.elementSum();

        if (n == 0 || mt == 0) {
            return new Ret.cacheIs(1.0, 0.0);
        }

        if (n < mt) {
            InputOutput.line_warning("cache_is", "Number of items (%d) less than cache capacity (%d).", n, mt);
            return new Ret.cacheIs(0.0, Double.NEGATIVE_INFINITY);
        }

        if (n == mt) {
            double E = Cache_erec.cache_erec(filteredGamma, m, filteredSigma, k).value();
            return new Ret.cacheIs(E, FastMath.log(E));
        }

        Matrix logGamma = new Matrix(n, h);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < h; j++) {
                logGamma.set(i, j, FastMath.log(filteredGamma.get(i, j) + 1e-300));
            }
        }

        double logMFact = 0.0;
        for (int j = 0; j < h; j++) {
            logMFact += factln((int) m.get(j));
        }

        double logCombinations = factln(n) - factln(mt) - factln(n - mt);

        double logMultinomial = factln(mt) - logMFact;

        double[] lZSamples = new double[samples];
        Random random = new Random();

        for (int s = 0; s < samples; s++) {
            int[] selected = sampleWithoutReplacement(n, mt, random);

            int[][] assignment = assignItemsToLevels(m, selected, random);

            double logStateProb = logMFact;
            boolean feasible = true;
            for (int j = 0; j < h && feasible; j++) {
                if (capped) {
                    double listCost = 0.0;
                    for (int item : assignment[j]) {
                        listCost += filteredSigma.get(item);
                    }
                    if (listCost > k.get(j)) {
                        feasible = false;
                        break;
                    }
                }
                for (int item : assignment[j]) {
                    logStateProb += logGamma.get(item, j);
                }
            }
            if (!feasible) {
                lZSamples[s] = Double.NEGATIVE_INFINITY; // I{S_v in O} = 0
                continue;
            }

            double logProposal = -logCombinations - logMultinomial;

            lZSamples[s] = logStateProb - logProposal;
        }

        double lE = logMeanExp(lZSamples);
        double E = FastMath.exp(lE);

        return new Ret.cacheIs(E, lE);
    }

    public static Ret.cacheIs cache_is(Matrix gamma, Matrix m) {
        return cache_is(gamma, m, 100000);
    }

    static int[] sampleWithoutReplacement(int n, int k, Random random) {
        int[] selected = new int[k];
        List<Integer> available = new ArrayList<Integer>();
        for (int i = 0; i < n; i++) available.add(i);

        for (int i = 0; i < k; i++) {
            int idx = random.nextInt(available.size());
            selected[i] = available.get(idx);
            available.remove(idx);
        }

        return selected;
    }

    static int[][] assignItemsToLevels(Matrix m, int[] selected, Random random) {
        int h = m.getNumElements();
        int[][] assignment = new int[h][];

        int[] shuffled = selected.clone();
        for (int i = shuffled.length - 1; i >= 1; i--) {
            int j = random.nextInt(i + 1);
            int temp = shuffled[i];
            shuffled[i] = shuffled[j];
            shuffled[j] = temp;
        }

        int idx = 0;
        for (int j = 0; j < h; j++) {
            int levelSize = (int) m.get(j);
            int[] arr = new int[levelSize];
            for (int k = 0; k < levelSize; k++) {
                arr[k] = shuffled[idx + k];
            }
            assignment[j] = arr;
            idx += levelSize;
        }

        return assignment;
    }

    static double factln(int n) {
        if (n < 0) return Double.NaN;
        if (n <= 1) return 0.0;
        if (n < 20) {
            double result = 0.0;
            for (int i = 2; i <= n; i++) {
                result += FastMath.log((double) i);
            }
            return result;
        }
        double nd = (double) n;
        return nd * FastMath.log(nd) - nd + 0.5 * FastMath.log(2 * FastMath.PI * nd);
    }

    static double logMeanExp(double[] x) {
        if (x.length == 0) return Double.NEGATIVE_INFINITY;
        double maxVal = x[0];
        for (int i = 1; i < x.length; i++) {
            if (x[i] > maxVal) maxVal = x[i];
        }
        if (maxVal == Double.NEGATIVE_INFINITY) return Double.NEGATIVE_INFINITY;

        double sum = 0.0;
        for (double xi : x) {
            sum += FastMath.exp(xi - maxVal);
        }

        return maxVal + FastMath.log(sum / x.length);
    }
}
