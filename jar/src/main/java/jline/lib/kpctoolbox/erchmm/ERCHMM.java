package jline.lib.kpctoolbox.erchmm;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.util.FastMath;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * Extended Renewal Continuous-time Hidden Markov Model (ER-CHMM) functions.
 */
public final class ERCHMM {
    private ERCHMM() {}

    public static ERCHMMFitResult erchmm_emfit(double[] trace, int[] orders) {
        return erchmm_emfit(trace, orders, 300, 1e-7, false);
    }

    public static ERCHMMFitResult erchmm_emfit(double[] trace, int[] orders, int iterMax, double iterTol, boolean verbose) {
        if (orders.length == 1) {
            int totalOrder = orders[0];
            MatrixCell bestMAP = null;
            double bestLogLi = Double.NEGATIVE_INFINITY;
            int[] bestOrders = new int[0];
            for (int numBranches = 2; numBranches <= totalOrder; numBranches++) {
                List<int[]> allCombinations = allErlangCombinations(numBranches, totalOrder);
                for (int[] combination : allCombinations) {
                    if (verbose) {
                        StringBuilder sb = new StringBuilder();
                        for (int i = 0; i < combination.length; i++) {
                            if (i > 0) sb.append(",");
                            sb.append(combination[i]);
                        }
                        System.out.println("Calculating with orders " + sb);
                    }
                    ERCHMMFitResult result = erchmm_emfit(trace, combination, iterMax, iterTol, verbose);
                    if (result.logLikelihood > bestLogLi) {
                        bestMAP = result.MAP;
                        bestLogLi = result.logLikelihood;
                        bestOrders = combination;
                    }
                }
            }
            if (verbose) {
                StringBuilder sb = new StringBuilder();
                for (int i = 0; i < bestOrders.length; i++) {
                    if (i > 0) sb.append(",");
                    sb.append(bestOrders[i]);
                }
                System.out.println("Best solution: log-likelihood=" + bestLogLi + ", orders=" + sb);
            }
            return new ERCHMMFitResult(bestMAP != null ? bestMAP : new MatrixCell(2), bestLogLi, bestOrders);
        }

        int M = orders.length;
        int K = trace.length;
        double[] piV = new double[M];
        for (int i = 0; i < M; i++) piV[i] = 1.0 / M;
        double[] lambda = new double[M];
        for (int i = 0; i < M; i++) lambda[i] = orders[i] * (i + 1);

        double traceMean = 0.0;
        for (double v : trace) traceMean += v;
        traceMean /= K;
        double piMean = 0.0;
        for (int i = 0; i < M; i++) piMean += piV[i] / (i + 1);
        for (int i = 0; i < M; i++) lambda[i] = lambda[i] * piMean / traceMean;

        double[][] T = new double[M][];
        for (int i = 0; i < M; i++) T[i] = piV.clone();

        double[][] F = new double[M][K];
        double[][] aLikelihoods = new double[K][M];
        double[][] bLikelihoods = new double[M][K];
        double[] aLikelihoodsScale = new double[K];
        double[] bLikelihoodsScale = new double[K];

        double ologli = 1.0;
        double logL = 0.0;
        int steps = 1;

        while (Math.abs((ologli - logL) / logL) > iterTol && steps < iterMax) {
            ologli = logL;
            for (int i = 0; i < M; i++) {
                for (int k = 0; k < K; k++) {
                    double lt = lambda[i] * trace[k];
                    int order = orders[i];
                    F[i][k] = (FastMath.pow(lt, order - 1) / factorial(order - 1) * lambda[i])
                            * FastMath.exp(-lambda[i] * trace[k]);
                }
            }

            double[] prevPi = piV.clone();
            double scaledPrev = 0.0;
            for (int k = 0; k < K; k++) {
                double[] temp = new double[M];
                for (int j = 0; j < M; j++) {
                    double sum = 0.0;
                    for (int i = 0; i < M; i++) sum += prevPi[i] * F[i][k] * T[i][j];
                    temp[j] = sum;
                }
                double sumTemp = 0.0;
                for (double v : temp) sumTemp += v;
                double scale = log2(sumTemp);
                double scaleFactor = FastMath.pow(2.0, -scale);
                for (int j = 0; j < M; j++) {
                    temp[j] *= scaleFactor;
                    prevPi[j] = temp[j];
                }
                aLikelihoodsScale[k] = scaledPrev + scale;
                aLikelihoods[k] = temp.clone();
                scaledPrev = aLikelihoodsScale[k];
            }

            double[][] aForwardLikelihoods = new double[K][];
            for (int k = 0; k < K; k++) {
                aForwardLikelihoods[k] = (k == 0) ? piV.clone() : aLikelihoods[k - 1].clone();
            }
            double[] aScaleV = new double[K];
            for (int k = 0; k < K; k++) aScaleV[k] = (k == 0) ? 0.0 : aLikelihoodsScale[k - 1];

            double[] nextB = new double[M];
            Arrays.fill(nextB, 1.0);
            scaledPrev = 0.0;
            for (int k = K - 1; k >= 0; k--) {
                double[] temp = new double[M];
                for (int i = 0; i < M; i++) {
                    double sum = 0.0;
                    for (int j = 0; j < M; j++) sum += F[i][k] * T[i][j] * nextB[j];
                    temp[i] = sum;
                }
                double sumTemp = 0.0;
                for (double v : temp) sumTemp += v;
                double scale = log2(sumTemp);
                double scaleFactor = FastMath.pow(2.0, -scale);
                for (int i = 0; i < M; i++) {
                    temp[i] *= scaleFactor;
                    nextB[i] = temp[i];
                }
                bLikelihoodsScale[k] = scaledPrev + scale;
                for (int i = 0; i < M; i++) bLikelihoods[i][k] = nextB[i];
                scaledPrev = bLikelihoodsScale[k];
            }

            double[][] bBackwardLikelihoods = new double[M][K];
            for (int i = 0; i < M; i++) {
                for (int k = 0; k < K; k++) {
                    bBackwardLikelihoods[i][k] = (k == K - 1) ? 1.0 : bLikelihoods[i][k + 1];
                }
            }
            double[] bScaleV = new double[K];
            for (int k = 0; k < K; k++) bScaleV[k] = (k == K - 1) ? 0.0 : bLikelihoodsScale[k + 1];

            double likelihoodValue = 0.0;
            for (int i = 0; i < M; i++) likelihoodValue += piV[i] * bLikelihoods[i][0];

            logL = (Math.log(likelihoodValue) + bLikelihoodsScale[0] * Math.log(2.0)) / K;
            double iLikelihood = 1.0 / likelihoodValue;

            double[][] likelihoodsMultiplied = new double[K][M];
            for (int k = 0; k < K; k++) {
                for (int m = 0; m < M; m++) {
                    likelihoodsMultiplied[k][m] = aForwardLikelihoods[k][m] * bLikelihoods[m][k];
                }
            }
            for (int k = 0; k < K; k++) {
                double sumRow = 0.0;
                for (double v : likelihoodsMultiplied[k]) sumRow += v;
                if (sumRow > 0) {
                    for (int m = 0; m < M; m++) likelihoodsMultiplied[k][m] /= sumRow;
                }
            }

            double[] numeratorEstimation = new double[M];
            double[] denominatorEstimation = new double[M];
            for (int m = 0; m < M; m++) {
                for (int k = 0; k < K; k++) {
                    numeratorEstimation[m] += likelihoodsMultiplied[k][m];
                    denominatorEstimation[m] += trace[k] * likelihoodsMultiplied[k][m];
                }
            }
            for (int m = 0; m < M; m++) {
                piV[m] = numeratorEstimation[m] / K;
                if (denominatorEstimation[m] > 0) {
                    lambda[m] = orders[m] * numeratorEstimation[m] / denominatorEstimation[m];
                }
            }

            double[][] densMutALikelihood = new double[K][M];
            for (int k = 0; k < K; k++) {
                for (int m = 0; m < M; m++) {
                    densMutALikelihood[k][m] = aForwardLikelihoods[k][m] * F[m][k];
                }
            }
            double[] summedLm = new double[K];
            for (int k = 0; k < K; k++) {
                summedLm[k] = iLikelihood * FastMath.pow(2.0, aScaleV[k] + bScaleV[k] - bLikelihoodsScale[0]);
            }
            for (int k = 0; k < K; k++) {
                for (int m = 0; m < M; m++) densMutALikelihood[k][m] *= summedLm[k];
            }

            for (int i = 0; i < M; i++) {
                for (int j = 0; j < M; j++) {
                    double sum = 0.0;
                    for (int k = 0; k < K; k++) sum += densMutALikelihood[k][i] * bBackwardLikelihoods[j][k];
                    T[i][j] = sum * T[i][j];
                }
            }
            for (int i = 0; i < M; i++) {
                double rowSum = 0.0;
                for (double v : T[i]) rowSum += v;
                if (rowSum > 0) {
                    for (int j = 0; j < M; j++) T[i][j] /= rowSum;
                }
            }
            steps++;
            if (verbose && steps % 50 == 0) {
                System.out.println("Num of iterations: " + steps + ", log-likelihood: " + logL);
            }
        }

        if (verbose) {
            System.out.println("Num of iterations: " + steps + ", log-likelihood: " + logL);
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < orders.length; i++) {
                if (i > 0) sb.append(",");
                sb.append(orders[i]);
            }
            System.out.println("EM algorithm terminated. (orders=" + sb + ")");
        }

        Matrix D0 = generateD0FromErlangs(lambda, orders);
        Matrix D1 = generateD1FromErlangs(lambda, orders, T);
        MatrixCell MAP = new MatrixCell(2);
        MAP.set(0, D0);
        MAP.set(1, D1);
        return new ERCHMMFitResult(MAP, logL, orders);
    }

    public static MatrixCell erchmm_emfit_simple(double[] trace, int[] orders) {
        return erchmm_emfit_simple(trace, orders, 300, 1e-7);
    }

    public static MatrixCell erchmm_emfit_simple(double[] trace, int[] orders, int iterMax, double iterTol) {
        return erchmm_emfit(trace, orders, iterMax, iterTol, false).MAP;
    }

    private static List<int[]> allErlangCombinations(int branches, int sumErlangs) {
        if (branches == 1) {
            return Arrays.asList(new int[]{sumErlangs});
        }
        List<int[]> result = new ArrayList<int[]>();
        for (int k1 = 1; k1 <= sumErlangs - branches + 1; k1++) {
            List<int[]> subCombinations = allErlangCombinations(branches - 1, sumErlangs - k1);
            for (int[] sub : subCombinations) {
                int[] combined = new int[sub.length + 1];
                System.arraycopy(sub, 0, combined, 0, sub.length);
                combined[sub.length] = k1;
                Arrays.sort(combined);
                boolean exists = false;
                for (int[] existing : result) {
                    if (Arrays.equals(existing, combined)) { exists = true; break; }
                }
                if (!exists) result.add(combined);
            }
        }
        return result;
    }

    private static Matrix generateD0FromErlangs(double[] lambda, int[] orders) {
        int n = 0;
        for (int o : orders) n += o;
        Matrix D0 = new Matrix(n, n);
        int startIdx = 0;
        for (int i = 0; i < lambda.length; i++) {
            int order = orders[i];
            double lam = lambda[i];
            for (int j = 0; j < order; j++) {
                D0.set(startIdx + j, startIdx + j, -lam);
                if (j < order - 1) D0.set(startIdx + j, startIdx + j + 1, lam);
            }
            startIdx += order;
        }
        return D0;
    }

    private static Matrix generateD1FromErlangs(double[] lambda, int[] orders, double[][] T) {
        int n = 0;
        for (int o : orders) n += o;
        int M = orders.length;
        Matrix D1 = new Matrix(n, n);
        int[] indicesTo = new int[M];
        int[] indicesFrom = new int[M];
        int cumSum = 0;
        for (int i = 0; i < M; i++) {
            indicesTo[i] = cumSum;
            cumSum += orders[i];
            indicesFrom[i] = cumSum - 1;
        }
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < M; j++) {
                D1.set(indicesFrom[i], indicesTo[j], lambda[i] * T[i][j]);
            }
        }
        return D1;
    }

    private static double factorial(int n) {
        if (n <= 1) return 1.0;
        double result = 1.0;
        for (int i = 2; i <= n; i++) result *= (double) i;
        return result;
    }

    private static double log2(double x) {
        return x > 0 ? Math.log(x) / Math.log(2.0) : 0.0;
    }
}
