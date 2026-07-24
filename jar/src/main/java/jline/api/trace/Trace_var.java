/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.trace;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import jline.util.Pair;

/**
 * Single trace variance computation.
 *
 * <p>Computes sample variance and related statistics for empirical trace data.
 *
 * @since LINE 3.0
 */
public final class Trace_var {
    private Trace_var() {}

    /**
     * Computes the variance of the trace data.
     */
    public static double trace_var(double[] trace) {
        double e1 = 0.0;
        double e2 = 0.0;
        for (int i = 0; i < trace.length; i++) {
            e2 += trace[i] * trace[i];
            e1 += trace[i];
        }
        return e2 / trace.length - (e1 / trace.length) * (e1 / trace.length);
    }

    /**
     * Computes the squared coefficient of variation for trace S.
     */
    public static double trace_scv(double[] S) {
        double mean = Trace_mean.trace_mean(S);
        return trace_var(S) / (mean * mean);
    }

    /**
     * Computes the autocorrelation function for trace S at the specified lags.
     */
    public static double[] trace_acf(double[] S, int[] lags) {
        int maxLag = 1;
        for (int l : lags) {
            if (l > maxLag) maxLag = l;
        }

        if (maxLag > S.length - 2) {
            List<Integer> filtered = new ArrayList<Integer>();
            for (int l : lags) {
                if (l <= S.length - 2 && l > 0) filtered.add(l);
            }
            if (filtered.isEmpty()) return new double[0];
            int[] flt = new int[filtered.size()];
            for (int i = 0; i < flt.length; i++) flt[i] = filtered.get(i);
            return trace_acf(S, flt);
        }

        double mean = Trace_mean.trace_mean(S);
        double[] centered = new double[S.length];
        for (int i = 0; i < S.length; i++) centered[i] = S[i] - mean;

        double[] autocov = autocov(centered);
        double[] rho = new double[lags.length];
        for (int i = 0; i < lags.length; i++) {
            rho[i] = autocov[lags[i]] / autocov[0];
        }
        return rho;
    }

    public static double[] trace_acf(double[] S) {
        return trace_acf(S, new int[]{1});
    }

    private static double[] autocov(double[] X) {
        int M = X.length;
        double[] acv = new double[M - 1];
        for (int p = 0; p < M - 1; p++) {
            double sum = 0.0;
            int count = 0;
            for (int i = 0; i < M - p; i++) {
                sum += X[i] * X[i + p];
                count++;
            }
            acv[p] = sum / count;
        }
        return acv;
    }

    /**
     * Estimates the auto-correlation decay rate of a trace.
     *
     * @return [GAMMA, RHO0, RESIDUALS]
     */
    public static double[] trace_gamma(double[] T, int limit) {
        double M1 = Trace_mean.trace_mean(T);
        double M2 = 0.0;
        for (double v : T) M2 += v * v;
        M2 /= T.length;

        int maxLag = Math.min(limit, T.length - 1);
        int[] lag = new int[maxLag];
        for (int i = 0; i < maxLag; i++) lag[i] = i + 1;
        double[] rho = trace_acf(T, lag);

        double VAR = M2 - M1 * M1;
        double SCV = VAR / (M1 * M1);
        double RHO0 = 0.5 * (1.0 - 1.0 / SCV);

        double bestGamma = 0.99;
        double minResiduals = Double.MAX_VALUE;

        for (int gamma = 990; gamma <= 999; gamma++) {
            double g = gamma / 1000.0;
            double residuals = 0.0;
            for (int i = 0; i < lag.length; i++) {
                double expected = RHO0 * Math.pow(g, lag[i]);
                residuals += (rho[i] - expected) * (rho[i] - expected);
            }
            if (residuals < minResiduals) {
                minResiduals = residuals;
                bestGamma = g;
            }
        }

        return new double[]{bestGamma, RHO0, minResiduals};
    }

    public static double[] trace_gamma(double[] T) {
        return trace_gamma(T, 1000);
    }

    /**
     * Computes the counting process of S.
     */
    public static int[] trace_iat2counts(double[] S, double scale) {
        int n = S.length;
        double[] CS = new double[n + 1];
        CS[0] = 0.0;
        for (int i = 1; i <= n; i++) {
            CS[i] = CS[i - 1] + S[i - 1];
        }

        ArrayList<Integer> C = new ArrayList<Integer>();
        for (int i = 0; i < n - 1; i++) {
            int cur = i;
            while (cur + 1 < n && CS[cur + 1] - CS[i] <= scale) {
                cur++;
            }
            C.add(cur - i);
            if (cur == n - 1) break;
        }
        int[] result = new int[C.size()];
        for (int i = 0; i < result.length; i++) result[i] = C.get(i);
        return result;
    }

    /**
     * Computes the Index of Dispersion for Intervals of a trace.
     */
    public static Pair<double[], int[]> trace_idi(double[] S, int[] kset, String option, int n) {
        ArrayList<Double> IDIk = new ArrayList<Double>();
        ArrayList<Integer> support = new ArrayList<Integer>();

        for (int k : kset) {
            if (option == null) {
                int supportVal = S.length - k - 1;
                support.add(supportVal);
                double[] Sk = new double[S.length - k];
                for (int t = 0; t < S.length - k - 1; t++) {
                    double sum = 0.0;
                    for (int j = t; j < t + k; j++) sum += S[j];
                    Sk[t] = sum;
                }
                double variance = trace_var(Sk);
                double mean = Trace_mean.trace_mean(Sk);
                IDIk.add(k * variance / (mean * mean));
            } else if (option.equals("aggregate")) {
                int keff = k / n;
                int supportVal = S.length / keff;
                support.add(supportVal);
                double[] Sk = new double[S.length - keff];
                for (int t = 0; t < S.length - keff - 1; t++) {
                    double sum = 0.0;
                    for (int j = t; j < t + keff; j++) sum += S[j];
                    Sk[t] = sum;
                }
                double variance = trace_var(Sk);
                double mean = Trace_mean.trace_mean(Sk);
                IDIk.add(k * variance / (mean * mean));
            }
        }

        double[] idiArr = new double[IDIk.size()];
        for (int i = 0; i < idiArr.length; i++) idiArr[i] = IDIk.get(i);
        int[] supArr = new int[support.size()];
        for (int i = 0; i < supArr.length; i++) supArr[i] = support.get(i);
        return new Pair<double[], int[]>(idiArr, supArr);
    }

    public static Pair<double[], int[]> trace_idi(double[] S, int[] kset) {
        return trace_idi(S, kset, null, 1);
    }

    /**
     * Computes the Index of Dispersion for Counts.
     */
    public static double trace_idc(double[] S) {
        int[] kset = new int[]{Math.min(1000, S.length / 30)};
        Pair<double[], int[]> result = trace_idi(S, kset);
        return result.getLeft()[0];
    }

    /**
     * Computes the probability mass function of discrete data.
     */
    public static Pair<double[], int[]> trace_pmf(int[] X) {
        java.util.TreeSet<Integer> set = new java.util.TreeSet<Integer>();
        for (int v : X) set.add(v);
        int[] uniqueValues = new int[set.size()];
        int idx = 0;
        for (Integer v : set) uniqueValues[idx++] = v;

        double[] pmf = new double[uniqueValues.length];
        for (int i = 0; i < uniqueValues.length; i++) {
            int count = 0;
            for (int v : X) if (v == uniqueValues[i]) count++;
            pmf[i] = (double) count / X.length;
        }
        return new Pair<double[], int[]>(pmf, uniqueValues);
    }

    /**
     * Shuffles the trace data randomly.
     */
    public static double[] trace_shuffle(double[] S) {
        double[] result = S.clone();
        Random random = new Random();
        for (int i = result.length - 1; i >= 0; i--) {
            int j = random.nextInt(i + 1);
            double temp = result[i];
            result[i] = result[j];
            result[j] = temp;
        }
        return result;
    }

    /**
     * Computes the bicovariance of the trace.
     */
    public static Pair<double[], int[][]> trace_bicov(double[] S, int[] GRID) {
        ArrayList<int[]> BiCovLags = new ArrayList<int[]>();
        ArrayList<Double> BiCov = new ArrayList<Double>();

        for (int i : GRID) {
            for (int j : GRID) {
                BiCovLags.add(new int[]{1, i, j});
            }
        }

        for (int[] lags : BiCovLags) {
            double jointMoment = trace_joint(S, lags, new int[]{1, 1, 1});
            BiCov.add(jointMoment);
        }

        double[] bcov = new double[BiCov.size()];
        for (int i = 0; i < bcov.length; i++) bcov[i] = BiCov.get(i);
        int[][] lagsArr = BiCovLags.toArray(new int[0][]);
        return new Pair<double[], int[][]>(bcov, lagsArr);
    }

    /**
     * Computes the counts in each bin with specified timescale.
     */
    public static Pair<int[], int[]> trace_iat2bins(double[] S, double scale) {
        int n = S.length;
        double[] CS = new double[n + 1];
        CS[0] = 0.0;
        for (int i = 1; i <= n; i++) CS[i] = CS[i - 1] + S[i - 1];

        int bins = (int) Math.ceil((CS[n] - CS[0]) / scale);
        int[] C = new int[bins];
        ArrayList<Integer> bC = new ArrayList<Integer>();

        int cur = 0;
        int last = 0;

        for (int i = 0; i < bins; i++) {
            if (cur >= n - 1) break;
            while (cur < n - 1 && CS[cur + 1] <= (i + 1) * scale) {
                cur++;
            }
            C[i] = cur - last;
            for (int j = 0; j < cur - last; j++) bC.add(i);
            last = cur;
        }
        int[] bcArr = new int[bC.size()];
        for (int i = 0; i < bcArr.length; i++) bcArr[i] = bC.get(i);
        return new Pair<int[], int[]>(C, bcArr);
    }

    /**
     * Computes joint moments E[X^{k_1}_{i} X^{k_2}_{i+j} ...] for a trace.
     */
    public static double trace_joint(double[] S, int[] lag, int[] order) {
        int[] sortedLag = lag.clone();
        Arrays.sort(sortedLag);
        int K = sortedLag.length;
        int[] adjustedLag = new int[K];
        int baseLag = sortedLag[0];
        for (int i = 0; i < K; i++) adjustedLag[i] = sortedLag[i] - baseLag;

        int maxLag = 0;
        for (int v : adjustedLag) if (v > maxLag) maxLag = v;
        int validLength = S.length - maxLag;
        if (validLength <= 0) return 0.0;

        double sum = 0.0;
        for (int i = 0; i < validLength; i++) {
            double product = 1.0;
            for (int j = 0; j < order.length; j++) {
                int idx = i + adjustedLag[Math.min(j, adjustedLag.length - 1)];
                if (idx < S.length) {
                    product *= Math.pow(S[idx], order[j]);
                }
            }
            sum += product;
        }
        return sum / validLength;
    }

    /**
     * Computes comprehensive summary statistics for a trace.
     *
     * @return [MEAN,SCV,MAD,SKEW,KURT,Q25,Q50,Q75,P95,MIN,MAX,IQR,ACF1-4,IDC_SCV_RATIO]
     */
    public static double[] trace_summary(double[] m) {
        double mean = Trace_mean.trace_mean(m);
        double scv = trace_scv(m);
        double[] sortedM = m.clone();
        Arrays.sort(sortedM);

        double q25 = percentile(sortedM, 25.0);
        double q50 = percentile(sortedM, 50.0);
        double q75 = percentile(sortedM, 75.0);
        double p95 = percentile(sortedM, 95.0);
        double min = sortedM[0];
        double max = sortedM[sortedM.length - 1];
        double iqr = q75 - q25;

        double[] absDevs = new double[m.length];
        for (int i = 0; i < m.length; i++) absDevs[i] = Math.abs(m[i] - q50);
        Arrays.sort(absDevs);
        double mad = absDevs[m.length / 2];

        double variance = trace_var(m);
        double std = Math.sqrt(variance);
        double skew = 0.0;
        double kurt = 0.0;
        for (double value : m) {
            double z = (value - mean) / std;
            skew += z * z * z;
            kurt += z * z * z * z;
        }
        skew /= m.length;
        kurt = kurt / m.length - 3.0;

        double[] acf = trace_acf(m, new int[]{1, 2, 3, 4});

        double idc = trace_idc(m);
        double idcScvRatio = idc / scv;

        return new double[]{mean, scv, mad, skew, kurt, q25, q50, q75, p95, min, max, iqr,
                acf[0], acf[1], acf[2], acf[3], idcScvRatio};
    }

    private static double percentile(double[] sortedData, double p) {
        double index = (p / 100.0) * (sortedData.length - 1);
        int lower = (int) Math.floor(index);
        int upper = (int) Math.ceil(index);
        if (lower == upper) return sortedData[lower];
        double weight = index - lower;
        return sortedData[lower] * (1 - weight) + sortedData[upper] * weight;
    }
}
