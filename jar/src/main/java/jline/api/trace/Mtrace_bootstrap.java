/**
 * @file Multi-class trace bootstrap resampling
 *
 * Implements bootstrap resampling methods for multi-class empirical trace data.
 *
 * @since LINE 3.0
 */
package jline.api.trace;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

public final class Mtrace_bootstrap {
    private Mtrace_bootstrap() {}

    public static BootstrapResults mtrace_bootstrap(double[] interArrivalTimes, int[] classLabels) {
        return mtrace_bootstrap(interArrivalTimes, classLabels, 1000, null, null);
    }

    public static BootstrapResults mtrace_bootstrap(double[] interArrivalTimes, int[] classLabels,
                                                    int numBootstraps) {
        return mtrace_bootstrap(interArrivalTimes, classLabels, numBootstraps, null, null);
    }

    public static BootstrapResults mtrace_bootstrap(double[] interArrivalTimes, int[] classLabels,
                                                    int numBootstraps, Integer blockSize) {
        return mtrace_bootstrap(interArrivalTimes, classLabels, numBootstraps, blockSize, null);
    }

    /**
     * Performs bootstrap resampling on a multi-class trace.
     */
    public static BootstrapResults mtrace_bootstrap(double[] interArrivalTimes, int[] classLabels,
                                                    int numBootstraps, Integer blockSize, Long seed) {
        if (interArrivalTimes.length != classLabels.length) {
            throw new IllegalArgumentException("Inter-arrival times and class labels must have same length");
        }

        Random random = (seed != null) ? new Random(seed) : new Random();
        int actualBlockSize = (blockSize != null) ? blockSize : determineOptimalBlockSize(interArrivalTimes);

        TraceStatistics originalStats = computeTraceStatistics(interArrivalTimes, classLabels);

        List<TraceStatistics> bootstrapStats = new ArrayList<TraceStatistics>();

        for (int b = 0; b < numBootstraps; b++) {
            double[] bootIAT;
            int[] bootLabels;
            if (actualBlockSize > 1) {
                Object[] sample = blockBootstrapSample(interArrivalTimes, classLabels, actualBlockSize, random);
                bootIAT = (double[]) sample[0];
                bootLabels = (int[]) sample[1];
            } else {
                Object[] sample = iidBootstrapSample(interArrivalTimes, classLabels, random);
                bootIAT = (double[]) sample[0];
                bootLabels = (int[]) sample[1];
            }

            TraceStatistics bootStats = computeTraceStatistics(bootIAT, bootLabels);
            bootstrapStats.add(bootStats);
        }

        Map<String, double[]> confidenceIntervals = computeConfidenceIntervals(bootstrapStats, 0.95);

        return new BootstrapResults(originalStats, bootstrapStats, confidenceIntervals, actualBlockSize);
    }

    private static Object[] blockBootstrapSample(double[] interArrivalTimes, int[] classLabels,
                                                 int blockSize, Random random) {
        int n = interArrivalTimes.length;
        int numBlocks = (n + blockSize - 1) / blockSize;

        List<Double> bootIAT = new ArrayList<Double>();
        List<Integer> bootLabels = new ArrayList<Integer>();

        for (int k = 0; k < numBlocks; k++) {
            int startIdx = random.nextInt(n - blockSize + 1);
            int endIdx = Math.min(startIdx + blockSize, n);

            for (int i = startIdx; i < endIdx; i++) {
                bootIAT.add(interArrivalTimes[i]);
                bootLabels.add(classLabels[i]);
            }
        }

        int targetSize = Math.min(bootIAT.size(), n);
        double[] iatArr = new double[targetSize];
        int[] labelsArr = new int[targetSize];
        for (int i = 0; i < targetSize; i++) {
            iatArr[i] = bootIAT.get(i);
            labelsArr[i] = bootLabels.get(i);
        }
        return new Object[] { iatArr, labelsArr };
    }

    private static Object[] iidBootstrapSample(double[] interArrivalTimes, int[] classLabels, Random random) {
        int n = interArrivalTimes.length;
        double[] bootIAT = new double[n];
        int[] bootLabels = new int[n];
        for (int i = 0; i < n; i++) {
            int idx = random.nextInt(n);
            bootIAT[i] = interArrivalTimes[idx];
            bootLabels[i] = classLabels[idx];
        }
        return new Object[] { bootIAT, bootLabels };
    }

    private static int determineOptimalBlockSize(double[] interArrivalTimes) {
        int n = interArrivalTimes.length;
        int maxLag = Math.min(50, n / 4);
        double[] autocorrs = new double[maxLag];

        double mean = 0.0;
        for (double v : interArrivalTimes) mean += v;
        mean /= n;

        double variance = 0.0;
        for (double v : interArrivalTimes) variance += (v - mean) * (v - mean);
        variance /= n;

        for (int lag = 1; lag < maxLag; lag++) {
            double sumProduct = 0.0;
            int count = 0;
            for (int i = 0; i < n - lag; i++) {
                sumProduct += (interArrivalTimes[i] - mean) * (interArrivalTimes[i + lag] - mean);
                count++;
            }
            autocorrs[lag] = (count > 0 && variance > 0) ? sumProduct / (count * variance) : 0.0;
        }

        int blockSize = 1;
        for (int lag = 1; lag < maxLag; lag++) {
            if (Math.abs(autocorrs[lag]) < 0.1) {
                blockSize = lag;
                break;
            }
        }

        return Math.max(1, Math.min(blockSize, n / 10));
    }

    private static TraceStatistics computeTraceStatistics(double[] interArrivalTimes, int[] classLabels) {
        int n = interArrivalTimes.length;
        int numClasses = 0;
        for (int v : classLabels) {
            if (v + 1 > numClasses) numClasses = v + 1;
        }

        double mean = 0.0;
        for (double v : interArrivalTimes) mean += v;
        mean /= n;

        double variance = 0.0;
        for (double v : interArrivalTimes) variance += (v - mean) * (v - mean);
        variance /= n;

        double arrivalRate = 1.0 / mean;
        double scv = variance / (mean * mean);

        double thirdMoment = 0.0;
        for (double v : interArrivalTimes) thirdMoment += Math.pow(v - mean, 3.0);
        thirdMoment /= n;

        double skewness = (variance > 0) ? thirdMoment / Math.pow(variance, 1.5) : 0.0;

        int[] classCounts = new int[numClasses];
        for (int label : classLabels) {
            if (label >= 0 && label < numClasses) {
                classCounts[label]++;
            }
        }
        double[] classProportions = new double[numClasses];
        for (int i = 0; i < numClasses; i++) {
            classProportions[i] = (double) classCounts[i] / n;
        }

        int maxLag = Math.min(10, n / 4);
        double[] lagCorrelations = new double[maxLag];

        for (int lag = 1; lag < maxLag; lag++) {
            double sumProduct = 0.0;
            int count = 0;
            for (int i = 0; i < n - lag; i++) {
                sumProduct += (interArrivalTimes[i] - mean) * (interArrivalTimes[i + lag] - mean);
                count++;
            }
            lagCorrelations[lag] = (count > 0 && variance > 0) ? sumProduct / (count * variance) : 0.0;
        }

        return new TraceStatistics(arrivalRate, scv, skewness, classProportions, lagCorrelations);
    }

    private static Map<String, double[]> computeConfidenceIntervals(List<TraceStatistics> bootstrapStats,
                                                                    double confidence) {
        double alpha = 1.0 - confidence;
        double lowerQuantile = alpha / 2.0;
        double upperQuantile = 1.0 - alpha / 2.0;

        Map<String, double[]> intervals = new HashMap<String, double[]>();

        List<Double> rates = new ArrayList<Double>();
        for (TraceStatistics s : bootstrapStats) rates.add(s.getArrivalRate());
        Collections.sort(rates);
        intervals.put("arrivalRate", new double[] { percentile(rates, lowerQuantile), percentile(rates, upperQuantile) });

        List<Double> scvs = new ArrayList<Double>();
        for (TraceStatistics s : bootstrapStats) scvs.add(s.getScv());
        Collections.sort(scvs);
        intervals.put("scv", new double[] { percentile(scvs, lowerQuantile), percentile(scvs, upperQuantile) });

        List<Double> skews = new ArrayList<Double>();
        for (TraceStatistics s : bootstrapStats) skews.add(s.getSkewness());
        Collections.sort(skews);
        intervals.put("skewness", new double[] { percentile(skews, lowerQuantile), percentile(skews, upperQuantile) });

        return intervals;
    }

    private static double percentile(List<Double> sortedData, double p) {
        if (sortedData.isEmpty()) return 0.0;
        int n = sortedData.size();
        double index = p * (n - 1);
        int lower = (int) index;
        int upper = Math.min(lower + 1, n - 1);
        double fraction = index - lower;
        return sortedData.get(lower) + fraction * (sortedData.get(upper) - sortedData.get(lower));
    }

    /**
     * Trace statistics container.
     */
    public static final class TraceStatistics {
        private final double arrivalRate;
        private final double scv;
        private final double skewness;
        private final double[] classProportions;
        private final double[] lagCorrelations;

        public TraceStatistics(double arrivalRate, double scv, double skewness,
                               double[] classProportions, double[] lagCorrelations) {
            this.arrivalRate = arrivalRate;
            this.scv = scv;
            this.skewness = skewness;
            this.classProportions = classProportions;
            this.lagCorrelations = lagCorrelations;
        }

        public double getArrivalRate() { return arrivalRate; }
        public double getScv() { return scv; }
        public double getSkewness() { return skewness; }
        public double[] getClassProportions() { return classProportions; }
        public double[] getLagCorrelations() { return lagCorrelations; }
    }

    /**
     * Bootstrap results container.
     */
    public static final class BootstrapResults {
        private final TraceStatistics original;
        private final List<TraceStatistics> bootstrapSamples;
        private final Map<String, double[]> confidenceIntervals;
        private final int blockSize;

        public BootstrapResults(TraceStatistics original, List<TraceStatistics> bootstrapSamples,
                                Map<String, double[]> confidenceIntervals, int blockSize) {
            this.original = original;
            this.bootstrapSamples = bootstrapSamples;
            this.confidenceIntervals = confidenceIntervals;
            this.blockSize = blockSize;
        }

        public TraceStatistics getOriginal() { return original; }
        public List<TraceStatistics> getBootstrapSamples() { return bootstrapSamples; }
        public Map<String, double[]> getConfidenceIntervals() { return confidenceIntervals; }
        public int getBlockSize() { return blockSize; }
    }
}
