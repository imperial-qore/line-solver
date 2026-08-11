/*
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
package jline.inference.api;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Random;

import jline.api.pfqn.mva.Pfqn_bs;
import jline.io.Ret;
import jline.util.matrix.Matrix;

/**
 * Gibbs Sampling demand estimation from trace data.
 */
public final class Infer_gibbs {
    private Infer_gibbs() {}

    public static double[] infer_gibbs(double[][][] data, int nbCores, double tol) {
        int datNeeded = 200000;
        int likelihoodSample = 5000;
        int nbSamples = 2000;

        int nbClasses = data[0].length - 1;
        int nbNodes = 2;
        int[] nbJobs = new int[nbClasses];

        AnalysisResult analysisResult = analyseData(data, nbJobs, nbClasses, nbNodes, datNeeded);
        Matrix prob = analysisResult.prob;
        int[] N = analysisResult.N;
        double[] N0 = analysisResult.N0;

        double usedCores = 0.0;
        for (int k = 0; k < prob.getNumRows(); k++) {
            double queueJobs = 0.0;
            for (int c = nbClasses; c < 2 * nbClasses; c++) {
                queueJobs += prob.get(k, c);
            }
            int lastCol = prob.getNumCols() - 1;
            double probVal = prob.get(k, lastCol);
            usedCores += (queueJobs > nbCores) ? nbCores * probVal : queueJobs * probVal;
        }
        int lastRow = prob.getNumRows() - 1;
        usedCores /= (1.0 - prob.get(lastRow, prob.getNumCols() - 1));

        double[] thinkTime = new double[nbClasses];
        for (int k = 0; k < nbClasses; k++) {
            double[] tput = data[6][k];
            if (tput != null && tput.length > 0) {
                double sum = 0.0;
                for (double v : tput) sum += v;
                double avg = sum / tput.length;
                thinkTime[k] = (N[k] - N0[k]) / avg;
            } else {
                thinkTime[k] = 1.0;
            }
        }

        double[] rangeSize = new double[nbClasses * (nbNodes - 1)];
        Arrays.fill(rangeSize, 1.0);

        double[] cumProb = new double[prob.getNumRows()];
        cumProb[0] = prob.get(0, prob.getNumCols() - 1);
        for (int k = 1; k < prob.getNumRows(); k++) {
            cumProb[k] = cumProb[k - 1] + prob.get(k, prob.getNumCols() - 1);
        }

        Random rng = new Random();
        double[][] testset = new double[likelihoodSample][nbClasses * nbNodes];
        for (int k = 0; k < likelihoodSample; k++) {
            double uniValue = rng.nextDouble();
            int index = -1;
            for (int i = 0; i < cumProb.length; i++) {
                if (cumProb[i] > uniValue) { index = i; break; }
            }
            if (index < 0) index = cumProb.length - 1;
            for (int c = 0; c < nbClasses * nbNodes; c++) {
                testset[k][c] = prob.get(index, c);
            }
        }

        int maxN = 0;
        for (int v : N) maxN += v;
        maxN += 1;
        double[] LV = new double[maxN + 1];
        LV[0] = 0.0;
        if (maxN >= 1) LV[1] = 0.0;
        for (int k = 2; k <= maxN; k++) {
            LV[k] = LV[k - 1] + Math.log(k - 1);
        }

        double sumA = 0.0;
        for (double[] row : testset) {
            for (int c = 0; c < nbClasses * nbNodes; c++) {
                int idx = (int) row[c];
                if (idx >= 0 && idx + 1 < LV.length) {
                    sumA += LV[idx + 1];
                }
            }
        }

        double logGInitial = 0.0;
        for (int k = 0; k < nbClasses; k++) {
            logGInitial += N[k] * Math.log(thinkTime[k] + 1e-15);
            for (int j = 1; j <= N[k]; j++) {
                logGInitial -= Math.log(j);
            }
        }

        double[][] smpl = new double[nbSamples][nbClasses * (nbNodes - 1)];
        int sampleIndex = 0;
        double[] demandOld = null;

        for (int k = 0; k < nbSamples / 50; k++) {
            for (int s = 0; s < 50; s++) {
                if (sampleIndex >= nbSamples) break;
                for (int h = 0; h < nbClasses * (nbNodes - 1); h++) {
                    double[] theta;
                    if (sampleIndex == 0) {
                        theta = new double[nbClasses * (nbNodes - 1)];
                        for (int i = 0; i < h; i++) theta[i] = smpl[0][i];
                    } else {
                        theta = new double[nbClasses * (nbNodes - 1)];
                        for (int i = 0; i < h; i++) theta[i] = smpl[sampleIndex][i];
                        for (int i = h; i < theta.length; i++) theta[i] = smpl[sampleIndex - 1][i];
                    }

                    GibbsResult result = gibbsSamplerSimple(
                            thinkTime, theta, testset, h, nbNodes, nbClasses,
                            N, logGInitial, tol, rangeSize[h], LV, sumA, rng);
                    smpl[sampleIndex][h] = result.value;
                    logGInitial = result.logGCurrent;
                    rangeSize[h] = result.rangeSizeDim * 2.0;
                }
                sampleIndex++;
            }

            if (k == 1) {
                demandOld = new double[nbClasses * (nbNodes - 1)];
                for (int dim = 0; dim < demandOld.length; dim++) {
                    double sum = 0.0;
                    for (int si = 50; si < sampleIndex; si++) sum += smpl[si][dim];
                    demandOld[dim] = sum / (sampleIndex - 50);
                }
            } else if (k > 1 && demandOld != null) {
                double[] demandNow = new double[nbClasses * (nbNodes - 1)];
                for (int dim = 0; dim < demandNow.length; dim++) {
                    double sum = 0.0;
                    for (int si = k * 50; si < sampleIndex; si++) sum += smpl[si][dim];
                    demandNow[dim] = sum / (sampleIndex - k * 50);
                }
                for (int dim = 0; dim < demandNow.length; dim++) {
                    demandNow[dim] = demandNow[dim] / (k + 1) + demandOld[dim] / (k + 1) * k;
                }

                double relChange = 0.0;
                for (int dim = 0; dim < demandNow.length; dim++) {
                    if (Math.abs(demandOld[dim]) > 1e-15) {
                        relChange += Math.abs((demandNow[dim] - demandOld[dim]) / demandOld[dim]);
                    }
                }
                relChange /= demandNow.length;

                if (relChange < tol) {
                    int cutoff = sampleIndex / 2;
                    double[] result = new double[nbClasses * (nbNodes - 1)];
                    for (int dim = 0; dim < result.length; dim++) {
                        double sum = 0.0;
                        for (int si = cutoff; si < sampleIndex; si++) sum += smpl[si][dim];
                        result[dim] = sum / (sampleIndex - cutoff) * usedCores;
                    }
                    return result;
                }
                demandOld = demandNow;
            }
        }

        int cutoff = sampleIndex / 2;
        double[] result = new double[nbClasses * (nbNodes - 1)];
        for (int dim = 0; dim < result.length; dim++) {
            double sum = 0.0;
            for (int si = cutoff; si < sampleIndex; si++) sum += smpl[si][dim];
            result[dim] = sum / (sampleIndex - cutoff) * usedCores;
        }
        return result;
    }

    public static double[] infer_gibbs(double[][][] data, int nbCores) {
        return infer_gibbs(data, nbCores, 1e-3);
    }

    private static final class AnalysisResult {
        final Matrix prob;
        final int[] N;
        final double[] N0;

        AnalysisResult(Matrix prob, int[] N, double[] N0) {
            this.prob = prob;
            this.N = N;
            this.N0 = N0;
        }
    }

    private static AnalysisResult analyseData(double[][][] data, int[] nbJobs, int nbClasses, int nbNodes, int datNeeded) {
        int K = nbClasses;
        int[] N = nbJobs.clone();
        double[] N0 = new double[K];

        List<Double> tempTS = new ArrayList<Double>();
        List<Integer> tempClass = new ArrayList<Integer>();
        List<Integer> tempLogger = new ArrayList<Integer>();

        for (int i = 0; i < K; i++) {
            double[] arvTimes = data[3][i];
            double[] respTimes = data[4][i];
            if (arvTimes == null || respTimes == null) continue;
            int tempLen = arvTimes.length;
            for (int j = 0; j < tempLen; j++) {
                tempTS.add(arvTimes[j]);
                tempClass.add(i);
                tempLogger.add(1);
                tempTS.add(arvTimes[j] + respTimes[j] * 1000);
                tempClass.add(i);
                tempLogger.add(2);
            }
        }

        Integer[] indices = new Integer[tempTS.size()];
        for (int i = 0; i < indices.length; i++) indices[i] = i;
        final List<Double> tsRef = tempTS;
        Arrays.sort(indices, new Comparator<Integer>() {
            @Override
            public int compare(Integer a, Integer b) {
                return Double.compare(tsRef.get(a), tsRef.get(b));
            }
        });

        double[] ts = new double[indices.length];
        int[] classId = new int[indices.length];
        int[] loggerId = new int[indices.length];
        for (int i = 0; i < indices.length; i++) {
            ts[i] = tempTS.get(indices[i]);
            classId[i] = tempClass.get(indices[i]);
            loggerId[i] = tempLogger.get(indices[i]);
        }

        int burnin = Math.max(0, ts.length - datNeeded);
        int totalLength = ts.length;

        int[][][] count = new int[totalLength][K][nbNodes];
        for (int k = 0; k < K; k++) count[0][k][0] = N[k];

        for (int i = 0; i < totalLength - 1; i++) {
            for (int k = 0; k < K; k++) {
                for (int nd = 0; nd < nbNodes; nd++) {
                    count[i + 1][k][nd] = count[i][k][nd];
                }
            }
            int cls = classId[i];
            int log = loggerId[i] - 1;
            count[i + 1][cls][log]--;
            if (log == nbNodes - 1) {
                count[i + 1][cls][0]++;
            } else {
                count[i + 1][cls][log + 1]++;
            }
        }

        int Nsum = 0;
        for (int v : N) Nsum += v;
        if (Nsum == 0) {
            for (int k = 0; k < K; k++) {
                int maxCount = 0;
                for (int i = 0; i < totalLength; i++) {
                    for (int nd = 0; nd < nbNodes; nd++) {
                        if (count[i][k][nd] > maxCount) maxCount = count[i][k][nd];
                    }
                }
                N[k] = maxCount;
            }
        }

        for (int i = 0; i < totalLength; i++) {
            for (int k = 0; k < K; k++) {
                count[i][k][0] += N[k];
            }
        }

        int flatSize = K * nbNodes;
        double[][] flatCount = new double[totalLength][flatSize];
        for (int i = 0; i < totalLength; i++) {
            for (int k = 0; k < K; k++) {
                for (int nd = 0; nd < nbNodes; nd++) {
                    flatCount[i][k + nd * K] = count[i][k][nd];
                }
            }
        }

        double[] timeInterval = new double[totalLength];
        timeInterval[0] = 0.0;
        for (int i = 1; i < totalLength; i++) {
            timeInterval[i] = ts[i] - ts[i - 1];
        }

        LinkedHashMap<String, double[]> stateMap = new LinkedHashMap<String, double[]>();
        for (int i = burnin; i < totalLength; i++) {
            StringBuilder sb = new StringBuilder();
            for (int j = 0; j < flatSize; j++) {
                if (j > 0) sb.append(",");
                sb.append(flatCount[i][j]);
            }
            String key = sb.toString();
            double[] entry = stateMap.get(key);
            if (entry == null) {
                entry = new double[flatSize + 1];
                stateMap.put(key, entry);
            }
            for (int j = 0; j < flatSize; j++) entry[j] = flatCount[i][j];
            entry[flatSize] += timeInterval[i];
        }

        List<double[]> sortedEntries = new ArrayList<double[]>(stateMap.values());
        final int fSize = flatSize;
        sortedEntries.sort(new Comparator<double[]>() {
            @Override
            public int compare(double[] a, double[] b) {
                for (int j = 0; j < fSize; j++) {
                    int cmp = Double.compare(a[j], b[j]);
                    if (cmp != 0) return cmp;
                }
                return 0;
            }
        });

        double obsLength = ts[ts.length - 1] - ts[burnin];
        Matrix probMatrix = new Matrix(sortedEntries.size(), flatSize + 1);
        int row = 0;
        for (double[] entry : sortedEntries) {
            for (int j = 0; j < flatSize; j++) {
                probMatrix.set(row, j, entry[j]);
            }
            probMatrix.set(row, flatSize, entry[flatSize] / obsLength);
            row++;
        }

        for (int k = 0; k < K; k++) {
            for (int r = 0; r < probMatrix.getNumRows(); r++) {
                N0[k] += probMatrix.get(r, probMatrix.getNumCols() - 1) * probMatrix.get(r, K + k);
            }
        }

        return new AnalysisResult(probMatrix, N, N0);
    }

    private static final class GibbsResult {
        final double value;
        final double logGCurrent;
        final double rangeSizeDim;

        GibbsResult(double value, double logGCurrent, double rangeSizeDim) {
            this.value = value;
            this.logGCurrent = logGCurrent;
            this.rangeSizeDim = rangeSizeDim;
        }
    }

    private static GibbsResult gibbsSamplerSimple(
            double[] thinkTime, double[] theta, double[][] testset, int index,
            int nbNodes, int nbClasses, int[] nbJobs, double logGInitial,
            double interval, double rangeSize, double[] LV, double sumA, Random rng) {
        List<Double> rangeList = new ArrayList<Double>();
        double cur = 0.0;
        while (cur <= rangeSize) {
            rangeList.add(cur);
            cur += interval;
        }
        double[] range = new double[rangeList.size()];
        for (int i = 0; i < range.length; i++) range[i] = rangeList.get(i);
        int rangeLen = range.length;

        double[][] x = new double[nbNodes][nbClasses];
        x[0] = thinkTime.clone();
        for (int i = 0; i < nbNodes - 1; i++) {
            for (int c = 0; c < nbClasses; c++) {
                x[i + 1][c] = theta[i * nbClasses + c];
            }
        }

        int indexI = index / nbClasses + 1;
        int indexJ = index % nbClasses;

        double[] logG = new double[rangeLen];

        Matrix demandMatrix = new Matrix(nbNodes - 1, nbClasses);
        for (int i = 0; i < nbNodes - 1; i++) {
            for (int c = 0; c < nbClasses; c++) {
                demandMatrix.set(i, c, x[i + 1][c]);
            }
        }
        Matrix popMatrix = new Matrix(1, nbClasses);
        for (int c = 0; c < nbClasses; c++) popMatrix.set(0, c, (double) nbJobs[c]);
        Matrix thinkMatrix = new Matrix(1, nbClasses);
        for (int c = 0; c < nbClasses; c++) thinkMatrix.set(0, c, thinkTime[c]);

        int indexPrevious = 0;
        double minDist = Double.MAX_VALUE;
        for (int i = 0; i < range.length; i++) {
            double dist = Math.abs(range[i] - theta[index]);
            if (dist < minDist) {
                minDist = dist;
                indexPrevious = i;
            }
        }
        logG[indexPrevious] = logGInitial;

        Ret.pfqnAMVA bsResult = Pfqn_bs.pfqn_bs(demandMatrix, popMatrix, thinkMatrix);
        Matrix QN = bsResult.Q;

        for (int i = indexPrevious - 1; i >= 0; i--) {
            x[indexI][indexJ] = range[i + 1];
            demandMatrix.set(indexI - 1, indexJ, x[indexI][indexJ]);
            Ret.pfqnAMVA result = Pfqn_bs.pfqn_bs(demandMatrix, popMatrix, thinkMatrix, interval, 1000, QN);
            QN = result.Q;
            double qnVal = QN.get(indexI - 1, indexJ);
            double ratio = 1.0 + qnVal / (range[i + 1] + 1e-15) * (-interval);
            logG[i] = ratio > 0 ? logG[i + 1] + Math.log(ratio) : logG[i + 1];
        }

        x[indexI][indexJ] = theta[index];
        demandMatrix.set(indexI - 1, indexJ, x[indexI][indexJ]);
        Ret.pfqnAMVA result2 = Pfqn_bs.pfqn_bs(demandMatrix, popMatrix, thinkMatrix);
        QN = result2.Q;

        for (int i = indexPrevious + 1; i < rangeLen; i++) {
            x[indexI][indexJ] = range[i - 1];
            demandMatrix.set(indexI - 1, indexJ, x[indexI][indexJ]);
            Ret.pfqnAMVA result = Pfqn_bs.pfqn_bs(demandMatrix, popMatrix, thinkMatrix, interval, 1000, QN);
            QN = result.Q;
            double qnVal = QN.get(indexI - 1, indexJ);
            double ratio = 1.0 + qnVal / (range[i - 1] + 1e-15) * interval;
            logG[i] = ratio > 0 ? logG[i - 1] + Math.log(ratio) : logG[i - 1];
        }

        double[] logProb = new double[rangeLen];
        for (int i = 0; i < rangeLen; i++) {
            double testSum = 0.0;
            for (double[] r : testset) {
                testSum += r[index + nbClasses];
            }
            logProb[i] = testSum * Math.log(range[i] + 1e-15) - logG[i] * testset.length;
        }

        double maxLogProb = Double.NEGATIVE_INFINITY;
        for (double v : logProb) if (v > maxLogProb) maxLogProb = v;
        double[] prob = new double[rangeLen];
        for (int i = 0; i < rangeLen; i++) prob[i] = Math.exp(logProb[i] - maxLogProb);
        double probSum = 0.0;
        for (double v : prob) probSum += v;
        for (int i = 0; i < prob.length; i++) prob[i] /= probSum;

        double[] cumProb = new double[rangeLen];
        cumProb[0] = prob[0];
        for (int i = 1; i < rangeLen; i++) cumProb[i] = cumProb[i - 1] + prob[i];

        int rangeSizeDimIdx = -1;
        for (int i = 0; i < rangeLen; i++) {
            if (cumProb[i] > 1.0 - 1e-10) { rangeSizeDimIdx = i; break; }
        }
        if (rangeSizeDimIdx < 0) rangeSizeDimIdx = rangeLen - 1;
        double rangeSizeDim = range[rangeSizeDimIdx] * 2.0;

        double randVar = rng.nextDouble();
        int indexProb = -1;
        for (int i = 0; i < cumProb.length; i++) {
            if (cumProb[i] > randVar) { indexProb = i; break; }
        }

        if (indexProb < 0) {
            return new GibbsResult(theta[index], logGInitial, rangeSizeDim);
        }
        return new GibbsResult(range[indexProb], logG[indexProb], rangeSizeDim);
    }
}
