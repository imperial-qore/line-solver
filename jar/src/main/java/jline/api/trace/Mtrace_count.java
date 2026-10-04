/**
 * @file Multi-class trace count statistics
 *
 * @since LINE 3.0
 */
package jline.api.trace;

public final class Mtrace_count {
    private Mtrace_count() {}

    /**
     * Computes count statistics from a multi-class trace over specified time windows.
     */
    public static CountStatistics mtrace_count(double[] interArrivalTimes, int[] classLabels,
                                               double windowSize, int classIndex) {
        if (interArrivalTimes.length != classLabels.length) {
            throw new IllegalArgumentException("Inter-arrival times and class labels must have same length");
        }
        if (!(windowSize > 0)) {
            throw new IllegalArgumentException("Window size must be positive");
        }

        double totalTime = 0.0;
        for (double t : interArrivalTimes) totalTime += t;
        int numWindows = Math.max(1, (int) (totalTime / windowSize));
        int maxLabel = 0;
        for (int l : classLabels) {
            if (l > maxLabel) maxLabel = l;
        }
        int numClasses = maxLabel + 1;

        int[][] counts;
        if (classIndex >= 0) {
            counts = new int[][] { generateCountProcess(interArrivalTimes, classLabels, windowSize, classIndex) };
        } else {
            counts = new int[numClasses][];
            for (int c = 0; c < numClasses; c++) {
                counts[c] = generateCountProcess(interArrivalTimes, classLabels, windowSize, c);
            }
        }

        SingleClassCountStats[] countStats = new SingleClassCountStats[counts.length];
        for (int i = 0; i < counts.length; i++) {
            countStats[i] = computeCountProcessStatistics(counts[i], windowSize);
        }

        return new CountStatistics(windowSize, numWindows, countStats, totalTime);
    }

    public static CountStatistics mtrace_count(double[] interArrivalTimes, int[] classLabels, double windowSize) {
        return mtrace_count(interArrivalTimes, classLabels, windowSize, -1);
    }

    private static int[] generateCountProcess(double[] interArrivalTimes, int[] classLabels,
                                              double windowSize, int targetClass) {
        double totalTime = 0.0;
        for (double t : interArrivalTimes) totalTime += t;
        int numWindows = Math.max(1, (int) (totalTime / windowSize));
        int[] counts = new int[numWindows];

        double currentTime = 0.0;
        int windowIndex = 0;

        for (int i = 0; i < interArrivalTimes.length; i++) {
            currentTime += interArrivalTimes[i];

            int targetWindow = Math.min((int) (currentTime / windowSize), numWindows - 1);

            if (classLabels[i] == targetClass) {
                for (int w = windowIndex; w <= targetWindow; w++) {
                    if (w < numWindows) {
                        counts[w]++;
                    }
                }
            }

            windowIndex = targetWindow;
        }

        return counts;
    }

    private static SingleClassCountStats computeCountProcessStatistics(int[] counts, double windowSize) {
        if (counts.length == 0) {
            return new SingleClassCountStats(0.0, 0.0, 0.0, 0.0, counts);
        }

        double sum = 0.0;
        for (int c : counts) sum += c;
        double mean = sum / counts.length;

        double variance;
        if (counts.length > 1) {
            double sqSum = 0.0;
            for (int c : counts) sqSum += (c - mean) * (c - mean);
            variance = sqSum / counts.length;
        } else {
            variance = 0.0;
        }

        double idc = (mean > 0) ? variance / mean : 0.0;

        double skewness;
        if (variance > 0 && counts.length > 2) {
            double thirdSum = 0.0;
            for (int c : counts) thirdSum += Math.pow(c - mean, 3.0);
            double thirdMoment = thirdSum / counts.length;
            skewness = thirdMoment / Math.pow(variance, 1.5);
        } else {
            skewness = 0.0;
        }

        return new SingleClassCountStats(mean, variance, idc, skewness, counts);
    }

    /**
     * Compute multi-scale count statistics
     */
    public static CountStatistics[] mtrace_count_multiscale(double[] interArrivalTimes,
                                                            int[] classLabels, double[] windowSizes) {
        CountStatistics[] result = new CountStatistics[windowSizes.length];
        for (int i = 0; i < windowSizes.length; i++) {
            result[i] = mtrace_count(interArrivalTimes, classLabels, windowSizes[i]);
        }
        return result;
    }

    /**
     * Convert inter-arrival times to count process
     */
    public static int[] mtrace_iat2counts(double[] interArrivalTimes, double timeScale) {
        double totalTime = 0.0;
        for (double t : interArrivalTimes) totalTime += t;
        int numWindows = Math.max(1, (int) (totalTime / timeScale));
        int[] counts = new int[numWindows];

        double currentTime = 0.0;
        int windowIndex = 0;

        for (double iat : interArrivalTimes) {
            currentTime += iat;
            int targetWindow = Math.min((int) (currentTime / timeScale), numWindows - 1);

            for (int w = windowIndex; w <= targetWindow; w++) {
                if (w < numWindows) {
                    counts[w]++;
                    break;
                }
            }

            windowIndex = targetWindow;
        }

        return counts;
    }
}
