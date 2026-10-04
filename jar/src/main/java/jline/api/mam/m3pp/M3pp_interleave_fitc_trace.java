/**
 * @file M3PP interleaved fitting from empirical trace data
 *
 * Implements M3PP interleaving and fitting directly from empirical trace data
 * with multiple arrival classes.
 *
 * @since LINE 3.0
 */
package jline.api.mam.m3pp;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import jline.util.Pair;
import jline.util.matrix.MatrixCell;

public final class M3pp_interleave_fitc_trace {
    private M3pp_interleave_fitc_trace() {}

    public static Pair<MatrixCell, List<MatrixCell>> m3pp_interleave_fitc_trace(
            double[] interArrivalTimes, int[] classLabels) {
        return m3pp_interleave_fitc_trace(interArrivalTimes, classLabels, null, null, null);
    }

    public static Pair<MatrixCell, List<MatrixCell>> m3pp_interleave_fitc_trace(
            double[] interArrivalTimes, int[] classLabels, Double t, Double tinf) {
        return m3pp_interleave_fitc_trace(interArrivalTimes, classLabels, t, tinf, null);
    }

    /**
     * Interleaves k M3PP to fit a multi-class trace with m classes.
     */
    public static Pair<MatrixCell, List<MatrixCell>> m3pp_interleave_fitc_trace(
            double[] interArrivalTimes, int[] classLabels,
            Double t, Double tinf, boolean[][] mapping) {
        if (interArrivalTimes.length != classLabels.length) {
            throw new IllegalArgumentException(
                    "Inter-arrival times and class labels must have the same length");
        }

        TreeSet<Integer> uniqueLabelsSet = new TreeSet<Integer>();
        for (int label : classLabels) uniqueLabelsSet.add(label);
        Integer[] uniqueLabels = uniqueLabelsSet.toArray(new Integer[0]);
        int m = uniqueLabels.length;

        double meanInterArrival = average(interArrivalTimes);
        double totalTime = sum(interArrivalTimes);
        double computedT = (t != null) ? t : (10.0 * meanInterArrival);
        double computedTinf = (tinf != null) ? tinf : Math.max(10.0 * computedT, totalTime / 100.0);

        double[][] countsT = computeMulticlassCountsAtScale(interArrivalTimes, classLabels, computedT);
        double[][] countsTinf = computeMulticlassCountsAtScale(interArrivalTimes, classLabels, computedTinf);

        double arrivalRate = 1.0 / meanInterArrival;

        double[] classRates = new double[m];
        for (int i = 0; i < m; i++) {
            int classCount = 0;
            for (int label : classLabels) {
                if (label == uniqueLabels[i]) classCount++;
            }
            classRates[i] = ((double) classCount / classLabels.length) * arrivalRate;
        }

        boolean[][] finalMapping = (mapping != null) ? mapping : computeCorrelationBasedMapping(countsTinf, m, 0.75);
        int k = finalMapping[0].length;

        validateMapping(finalMapping, m);

        System.out.println("Fitting " + m + " classes with " + k + " M3PP(2,m_j) processes");

        boolean[][] filters = new boolean[k][m];
        for (int j = 0; j < k; j++) {
            for (int i = 0; i < m; i++) {
                filters[j][i] = finalMapping[i][j];
            }
        }

        double[] processRates = new double[k];
        for (int j = 0; j < k; j++) {
            double sum = 0.0;
            for (int i = 0; i < m; i++) {
                if (filters[j][i]) sum += classRates[i];
            }
            processRates[j] = sum;
        }

        double[][] classRatesPerProcess = new double[k][];
        for (int j = 0; j < k; j++) {
            List<Double> arr = new ArrayList<Double>();
            for (int i = 0; i < m; i++) {
                if (filters[j][i]) arr.add(classRates[i]);
            }
            classRatesPerProcess[j] = new double[arr.size()];
            for (int i = 0; i < arr.size(); i++) classRatesPerProcess[j][i] = arr.get(i);
        }

        double[] idcT = new double[k];
        double[] idcTinf = new double[k];

        for (int j = 0; j < k; j++) {
            double[] aggregateCountsT = new double[countsT.length];
            for (int wIdx = 0; wIdx < countsT.length; wIdx++) {
                double s = 0.0;
                for (int classIdx = 0; classIdx < m; classIdx++) {
                    if (filters[j][classIdx]) s += countsT[wIdx][classIdx];
                }
                aggregateCountsT[wIdx] = s;
            }

            double[] aggregateCountsTinf = new double[countsTinf.length];
            for (int wIdx = 0; wIdx < countsTinf.length; wIdx++) {
                double s = 0.0;
                for (int classIdx = 0; classIdx < m; classIdx++) {
                    if (filters[j][classIdx]) s += countsTinf[wIdx][classIdx];
                }
                aggregateCountsTinf[wIdx] = s;
            }

            double meanT = average(aggregateCountsT);
            double varT = 0.0;
            for (double v : aggregateCountsT) varT += (v - meanT) * (v - meanT);
            varT /= aggregateCountsT.length;
            idcT[j] = varT / (processRates[j] * computedT);

            double meanTinf = average(aggregateCountsTinf);
            double varTinf = 0.0;
            for (double v : aggregateCountsTinf) varTinf += (v - meanTinf) * (v - meanTinf);
            varTinf /= aggregateCountsTinf.length;
            idcTinf[j] = varTinf / (processRates[j] * computedTinf);
        }

        double[][] gtc = new double[k][];
        for (int j = 0; j < k; j++) {
            int numClassesInProcess = 0;
            for (boolean f : filters[j]) if (f) numClassesInProcess++;
            double[] result = new double[numClassesInProcess];

            int resultIdx = 0;
            for (int i = 0; i < m; i++) {
                if (filters[j][i]) {
                    double[] classCounts = new double[countsT.length];
                    for (int wIdx = 0; wIdx < countsT.length; wIdx++) {
                        classCounts[wIdx] = countsT[wIdx][i];
                    }
                    double classMean = average(classCounts);
                    double classVariance = 0.0;
                    for (double v : classCounts) classVariance += (v - classMean) * (v - classMean);
                    classVariance /= classCounts.length;

                    double[] otherClassesCounts = new double[countsT.length];
                    for (int wIdx = 0; wIdx < countsT.length; wIdx++) {
                        double s = 0.0;
                        for (int classIdx = 0; classIdx < m; classIdx++) {
                            if (filters[j][classIdx] && classIdx != i) {
                                s += countsT[wIdx][classIdx];
                            }
                        }
                        otherClassesCounts[wIdx] = s;
                    }

                    double covariance = 0.0;
                    boolean any = false;
                    for (double v : otherClassesCounts) if (v != 0.0) { any = true; break; }
                    if (any) {
                        double otherMean = average(otherClassesCounts);
                        for (int idx = 0; idx < classCounts.length; idx++) {
                            covariance += (classCounts[idx] - classMean) * (otherClassesCounts[idx] - otherMean);
                        }
                        covariance /= classCounts.length;
                    }

                    result[resultIdx] = classVariance + covariance;
                    resultIdx++;
                }
            }
            gtc[j] = result;
        }

        return M3pp_interleave_fitc.m3pp_interleave_fitc(processRates, idcT, idcTinf,
                classRatesPerProcess, gtc, computedT, computedTinf, null, false);
    }

    private static double[][] computeMulticlassCountsAtScale(
            double[] interArrivalTimes, int[] classLabels, double scale) {
        TreeSet<Integer> uniqueLabelsSet = new TreeSet<Integer>();
        for (int label : classLabels) uniqueLabelsSet.add(label);
        Integer[] uniqueLabels = uniqueLabelsSet.toArray(new Integer[0]);
        int m = uniqueLabels.length;
        double totalTime = sum(interArrivalTimes);
        int numWindows = Math.max(1, (int) (totalTime / scale));

        double[][] counts = new double[numWindows][m];

        double currentTime = 0.0;
        int windowIndex = 0;

        for (int i = 0; i < interArrivalTimes.length; i++) {
            currentTime += interArrivalTimes[i];

            double targetWindowEnd = (windowIndex + 1) * scale;
            int classIdx = indexOf(uniqueLabels, classLabels[i]);
            if (currentTime <= targetWindowEnd) {
                if (classIdx >= 0) {
                    counts[windowIndex][classIdx] += 1.0;
                }
            } else {
                windowIndex++;
                if (windowIndex >= numWindows) break;
                if (classIdx >= 0) {
                    counts[windowIndex][classIdx] = 1.0;
                }
            }
        }
        return counts;
    }

    private static int indexOf(Integer[] arr, int v) {
        for (int i = 0; i < arr.length; i++) {
            if (arr[i].intValue() == v) return i;
        }
        return -1;
    }

    private static boolean[][] computeCorrelationBasedMapping(double[][] counts, int m, double threshold) {
        double[][] covariance = computeCovarianceMatrix(counts);
        double[] variance = new double[m];
        for (int i = 0; i < m; i++) {
            double[] classCounts = new double[counts.length];
            for (int j = 0; j < counts.length; j++) classCounts[j] = counts[j][i];
            double mean = average(classCounts);
            double v = 0.0;
            for (double cc : classCounts) v += (cc - mean) * (cc - mean);
            variance[i] = v / classCounts.length;
        }

        Set<Integer> pool = new HashSet<Integer>();
        for (int i = 0; i < m; i++) pool.add(i);
        List<List<Integer>> groups = new ArrayList<List<Integer>>();

        while (!pool.isEmpty()) {
            int pivot = -1;
            double maxVar = Double.NEGATIVE_INFINITY;
            for (int i : pool) {
                if (variance[i] > maxVar) {
                    maxVar = variance[i];
                    pivot = i;
                }
            }
            if (pivot < 0) break;

            List<Integer> currentGroup = new ArrayList<Integer>();
            currentGroup.add(pivot);
            pool.remove(pivot);

            Set<Integer> toRemove = new HashSet<Integer>();
            for (int h : pool) {
                double correlation = (variance[pivot] > 0 && variance[h] > 0)
                        ? covariance[pivot][h] / Math.sqrt(variance[pivot] * variance[h])
                        : 0.0;
                if (correlation >= threshold) {
                    currentGroup.add(h);
                    toRemove.add(h);
                }
            }
            pool.removeAll(toRemove);
            groups.add(currentGroup);
        }

        int k = groups.size();
        boolean[][] mapping = new boolean[m][k];
        for (int j = 0; j < k; j++) {
            for (int classIdx : groups.get(j)) {
                mapping[classIdx][j] = true;
            }
        }
        return mapping;
    }

    private static double[][] computeCovarianceMatrix(double[][] counts) {
        int m = counts[0].length;
        double[] means = new double[m];
        for (int i = 0; i < m; i++) {
            double s = 0.0;
            for (int j = 0; j < counts.length; j++) s += counts[j][i];
            means[i] = s / counts.length;
        }

        double[][] covariance = new double[m][m];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                double s = 0.0;
                for (int kk = 0; kk < counts.length; kk++) {
                    s += (counts[kk][i] - means[i]) * (counts[kk][j] - means[j]);
                }
                covariance[i][j] = s / counts.length;
            }
        }
        return covariance;
    }

    private static void validateMapping(boolean[][] mapping, int m) {
        if (mapping.length != m) {
            throw new IllegalArgumentException("Number of classes does not match mapping");
        }
        for (int i = 0; i < m; i++) {
            int count = 0;
            for (boolean b : mapping[i]) if (b) count++;
            if (count != 1) {
                throw new IllegalArgumentException(
                        "Invalid mapping: class " + i + " mapped to " + count + " processes");
            }
        }
    }

    private static double average(double[] arr) {
        double s = 0.0;
        for (double v : arr) s += v;
        return s / arr.length;
    }

    private static double sum(double[] arr) {
        double s = 0.0;
        for (double v : arr) s += v;
        return s;
    }
}
