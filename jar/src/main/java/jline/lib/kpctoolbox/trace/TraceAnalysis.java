/**
 * Trace analysis functions for the KPC-Toolbox.
 *
 * Ported from MATLAB: matlab/lib/kpctoolbox/trace/
 */
package jline.lib.kpctoolbox.trace;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.TreeSet;

import jline.util.Pair;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.descriptive.moment.Kurtosis;
import org.apache.commons.math3.stat.descriptive.moment.Skewness;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;
import org.apache.commons.math3.util.FastMath;

public final class TraceAnalysis {
    private TraceAnalysis() {}

    public static double trace_mean(double[] S) {
        return StatUtils.mean(S);
    }

    public static double trace_var(double[] S) {
        return StatUtils.variance(S);
    }

    public static double trace_scv(double[] S) {
        double mean = trace_mean(S);
        double variance = trace_var(S);
        return variance / (mean * mean);
    }

    public static double[] autocov(double[] X) {
        int M = X.length;
        if (M < 2) {
            throw new IllegalArgumentException("Trace too short");
        }
        double mean = StatUtils.mean(X);

        int padLen = nextPowerOf2(2 * M);
        double[] xPad = new double[padLen];
        for (int i = 0; i < M; i++) {
            xPad[i] = X[i] - mean;
        }

        FastFourierTransformer fft = new FastFourierTransformer(DftNormalization.STANDARD);
        Complex[] xHat = fft.transform(xPad, TransformType.FORWARD);

        Complex[] product = new Complex[padLen];
        for (int i = 0; i < padLen; i++) {
            double re = xHat[i].getReal();
            double im = xHat[i].getImaginary();
            product[i] = new Complex(re * re + im * im, 0.0);
        }

        Complex[] result = fft.transform(product, TransformType.INVERSE);

        double[] acv = new double[M - 1];
        for (int i = 0; i < M - 1; i++) {
            acv[i] = result[i].getReal() / (M - i);
        }
        return acv;
    }

    private static int nextPowerOf2(int n) {
        int power = 1;
        while (power < n) power *= 2;
        return power;
    }

    public static double[] trace_acf(double[] S, int[] lags) {
        if (lags.length == 0) return new double[0];

        int n = S.length;
        int maxLag = Integer.MIN_VALUE;
        for (int l : lags) if (l > maxLag) maxLag = l;
        if (maxLag == Integer.MIN_VALUE) maxLag = 1;

        if (maxLag > n - 2) {
            List<Integer> validLags = new ArrayList<Integer>();
            for (int l : lags) if (l <= n - 2) validLags.add(l);
            if (validLags.isEmpty()) return new double[0];
            int[] arr = new int[validLags.size()];
            for (int i = 0; i < arr.length; i++) arr[i] = validLags.get(i);
            return trace_acf(S, arr);
        }

        double[] acv = autocov(S);
        double[] rho = new double[lags.length];
        double variance = acv[0];
        if (variance == 0.0) return new double[lags.length];

        for (int idx = 0; idx < lags.length; idx++) {
            int lag = lags[idx];
            if (lag > 0 && lag < acv.length) {
                rho[idx] = acv[lag] / variance;
            } else if (lag == 0) {
                rho[idx] = 1.0;
            }
        }
        return rho;
    }

    public static double trace_acf(double[] S) {
        return trace_acf(S, new int[]{1})[0];
    }

    public static double trace_skew(double[] S) {
        List<Double> filtered = new ArrayList<Double>();
        for (double v : S) if (!Double.isNaN(v)) filtered.add(v);
        int n = filtered.size();
        if (n < 3) return Double.NaN;

        double[] arr = new double[n];
        for (int i = 0; i < n; i++) arr[i] = filtered.get(i);
        double mean = StatUtils.mean(arr);
        double s2 = 0.0, m3 = 0.0;
        for (double x : arr) {
            double res = x - mean;
            s2 += res * res;
            m3 += res * res * res;
        }
        s2 /= n;
        m3 /= n;
        double correction = FastMath.sqrt((n - 1.0) / n) * n / (n - 2.0);
        return m3 * correction / FastMath.pow(s2, 1.5);
    }

    public static double trace_joint(double[] S, int[] lags, int[] orders) {
        if (lags.length != orders.length) {
            throw new IllegalArgumentException("lags and orders must have same length");
        }
        int[] cumLags = new int[lags.length];
        int sum = 0;
        for (int i = 0; i < lags.length; i++) {
            sum += lags[i];
            cumLags[i] = sum;
        }
        int minLag = Integer.MAX_VALUE;
        for (int v : cumLags) if (v < minLag) minLag = v;
        if (minLag == Integer.MAX_VALUE) minLag = 0;
        for (int i = 0; i < cumLags.length; i++) cumLags[i] -= minLag;

        int maxLag = 0;
        for (int v : cumLags) if (v > maxLag) maxLag = v;

        int n = S.length;
        if (maxLag >= n) {
            throw new IllegalArgumentException("Lag exceeds trace length");
        }

        double total = 0.0;
        int validCount = n - maxLag;
        for (int t = 0; t < validCount; t++) {
            double product = 1.0;
            for (int i = 0; i < cumLags.length; i++) {
                product *= FastMath.pow(S[t + cumLags[i]], orders[i]);
            }
            total += product;
        }
        return total / validCount;
    }

    public static Pair<double[], int[][]> trace_bicov(double[] S, int[] grid) {
        List<int[]> lagPairs = new ArrayList<int[]>();
        for (int i : grid) {
            for (int j : grid) {
                lagPairs.add(new int[]{1, i, j});
            }
        }
        double[] bicov = new double[lagPairs.size()];
        for (int idx = 0; idx < lagPairs.size(); idx++) {
            bicov[idx] = trace_joint(S, lagPairs.get(idx), new int[]{1, 1, 1});
        }
        int[][] arr = lagPairs.toArray(new int[0][]);
        return new Pair<double[], int[][]>(bicov, arr);
    }

    public static double[] trace_idi(double[] S, int[] kset) {
        double[] idiValues = new double[kset.length];
        for (int idx = 0; idx < kset.length; idx++) {
            int k = kset[idx];
            if (k >= S.length) {
                idiValues[idx] = Double.NaN;
                continue;
            }
            int numSums = S.length - k;
            if (numSums <= 0) {
                idiValues[idx] = Double.NaN;
                continue;
            }
            double[] Sk = new double[numSums];
            for (int t = 0; t < numSums; t++) {
                double sum = 0.0;
                for (int i = 0; i < k; i++) sum += S[t + i];
                Sk[t] = sum;
            }
            double meanSk = StatUtils.mean(Sk);
            double varSk = StatUtils.variance(Sk);
            if (meanSk != 0.0) {
                idiValues[idx] = k * varSk / (meanSk * meanSk);
            } else {
                idiValues[idx] = Double.NaN;
            }
        }
        return idiValues;
    }

    public static double trace_idi(double[] S) {
        int k = Math.min(1000, S.length / 30);
        if (k <= 0) return Double.NaN;
        return trace_idi(S, new int[]{k})[0];
    }

    public static double trace_idc(double[] S) {
        return trace_idi(S);
    }

    public static double trace_gamma(double[] T) {
        return trace_gamma(T, 1000);
    }

    public static double trace_gamma(double[] T, int limit) {
        double M1 = trace_mean(T);
        double[] sq = new double[T.length];
        for (int i = 0; i < T.length; i++) sq[i] = T[i] * T[i];
        double M2 = StatUtils.mean(sq);

        int maxLag = Math.min(limit, T.length - 2);
        if (maxLag <= 0) return Double.NaN;

        int[] lags = new int[maxLag];
        for (int i = 0; i < maxLag; i++) lags[i] = i + 1;
        double[] rho = trace_acf(T, lags);

        double variance = M2 - M1 * M1;
        double scv = variance / (M1 * M1);
        double rho0 = 0.5 * (1 - 1 / scv);

        double sumLogRatio = 0.0;
        int count = 0;
        for (int i = 0; i < lags.length; i++) {
            double rhoK = rho[i];
            int k = lags[i];
            if (rhoK > 0 && rhoK < 1) {
                if (rho0 > 0 && rhoK / rho0 > 0) {
                    sumLogRatio += FastMath.log(rhoK / rho0) / k;
                    count++;
                }
            }
        }

        if (count > 0) {
            double v = FastMath.exp(sumLogRatio / count);
            if (v < 0.0) v = 0.0;
            if (v > 1.0) v = 1.0;
            return v;
        } else {
            return 0.5;
        }
    }

    public static double[] trace_shuffle(double[] S) {
        double[] result = S.clone();
        Random random = new Random();
        for (int i = result.length - 1; i >= 1; i--) {
            int j = random.nextInt(i + 1);
            double temp = result[i];
            result[i] = result[j];
            result[j] = temp;
        }
        return result;
    }

    public static int[] trace_iat2counts(double[] S, double scale) {
        int n = S.length;
        if (n <= 1) return new int[0];

        double[] CS = new double[n];
        CS[0] = S[0];
        for (int i = 1; i < n; i++) {
            CS[i] = CS[i - 1] + S[i];
        }

        List<Integer> counts = new ArrayList<Integer>();
        for (int i = 0; i < n - 1; i++) {
            int cur = i;
            while (cur + 1 < n && CS[cur + 1] - CS[i] <= scale) {
                cur++;
            }
            counts.add(cur - i);
            if (cur >= n - 1) break;
        }
        int[] result = new int[counts.size()];
        for (int i = 0; i < result.length; i++) result[i] = counts.get(i);
        return result;
    }

    public static Pair<int[], int[]> trace_iat2bins(double[] S, double scale) {
        int n = S.length;
        if (n == 0) return new Pair<int[], int[]>(new int[0], new int[0]);

        double[] CS = new double[n];
        CS[0] = S[0];
        for (int i = 1; i < n; i++) CS[i] = CS[i - 1] + S[i];

        double totalTime = CS[n - 1] - CS[0];
        int numBins = (int) Math.ceil(totalTime / scale);
        if (numBins <= 0) return new Pair<int[], int[]>(new int[0], new int[0]);

        int[] counts = new int[numBins];
        List<Integer> binMembership = new ArrayList<Integer>();

        int cur = 0;
        int last = 0;
        for (int binIdx = 0; binIdx < numBins; binIdx++) {
            double binEnd = (binIdx + 1) * scale;
            while (cur + 1 < n && CS[cur + 1] <= binEnd) cur++;
            counts[binIdx] = cur - last;
            for (int j = 0; j < cur - last; j++) {
                binMembership.add(binIdx + 1);
            }
            last = cur;
        }
        int[] bm = new int[binMembership.size()];
        for (int i = 0; i < bm.length; i++) bm[i] = binMembership.get(i);
        return new Pair<int[], int[]>(counts, bm);
    }

    public static Pair<double[], double[]> trace_pmf(double[] X) {
        double[] sorted = X.clone();
        Arrays.sort(sorted);
        TreeSet<Double> uniqueSet = new TreeSet<Double>();
        for (double v : sorted) uniqueSet.add(v);
        Double[] unique = uniqueSet.toArray(new Double[0]);
        double n = X.length;
        double[] pmf = new double[unique.length];
        double[] uniqueArr = new double[unique.length];
        for (int idx = 0; idx < unique.length; idx++) {
            double value = unique[idx];
            int c = 0;
            for (double v : X) if (v == value) c++;
            pmf[idx] = c / n;
            uniqueArr[idx] = value;
        }
        return new Pair<double[], double[]>(pmf, uniqueArr);
    }

    /**
     * Summary statistics for a trace.
     */
    public static class TraceSummary {
        public final double mean;
        public final double scv;
        public final double mad;
        public final double skewness;
        public final double kurtosis;
        public final double[] quartiles;
        public final double percentile95;
        public final double min;
        public final double max;
        public final double iqr;
        public final double[] acf;
        public final double idc;
        public final int length;

        public TraceSummary(double mean, double scv, double mad, double skewness, double kurtosis,
                            double[] quartiles, double percentile95, double min, double max,
                            double iqr, double[] acf, double idc, int length) {
            this.mean = mean;
            this.scv = scv;
            this.mad = mad;
            this.skewness = skewness;
            this.kurtosis = kurtosis;
            this.quartiles = quartiles;
            this.percentile95 = percentile95;
            this.min = min;
            this.max = max;
            this.iqr = iqr;
            this.acf = acf;
            this.idc = idc;
            this.length = length;
        }

        @Override
        public boolean equals(Object other) {
            if (this == other) return true;
            if (other == null || getClass() != other.getClass()) return false;
            TraceSummary o = (TraceSummary) other;
            return Double.compare(mean, o.mean) == 0 && Double.compare(scv, o.scv) == 0;
        }

        @Override
        public int hashCode() {
            int result = Double.valueOf(mean).hashCode();
            result = 31 * result + Double.valueOf(scv).hashCode();
            return result;
        }
    }

    public static TraceSummary trace_summary(double[] m) {
        double mean = trace_mean(m);
        double scv = trace_scv(m);
        double skewness = new Skewness().evaluate(m);
        double kurtosis = new Kurtosis().evaluate(m);

        Percentile percentile = new Percentile();
        percentile.setData(m);
        double q25 = percentile.evaluate(25.0);
        double q50 = percentile.evaluate(50.0);
        double q75 = percentile.evaluate(75.0);
        double p95 = percentile.evaluate(95.0);
        double iqr = q75 - q25;

        double median = q50;
        double[] absDeviations = new double[m.length];
        for (int i = 0; i < m.length; i++) absDeviations[i] = FastMath.abs(m[i] - median);
        double mad = new Percentile().evaluate(absDeviations, 50.0);

        int[] lags = new int[10];
        for (int i = 0; i < 10; i++) lags[i] = i + 1;
        double[] acf = trace_acf(m, lags);
        double idc = trace_idc(m);

        double minVal = Double.POSITIVE_INFINITY, maxVal = Double.NEGATIVE_INFINITY;
        for (double v : m) {
            if (v < minVal) minVal = v;
            if (v > maxVal) maxVal = v;
        }
        if (m.length == 0) {
            minVal = Double.NaN;
            maxVal = Double.NaN;
        }

        return new TraceSummary(mean, scv, mad, skewness, kurtosis,
                new double[]{q25, q50, q75}, p95, minVal, maxVal,
                iqr, acf, idc, m.length);
    }

    public static double[] mtrace_mean(List<double[]> traces) {
        double[] result = new double[traces.size()];
        for (int i = 0; i < traces.size(); i++) {
            result[i] = trace_mean(traces.get(i));
        }
        return result;
    }
}
