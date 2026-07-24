/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.pfqn.mva;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import jline.io.Ret;
import jline.lang.constant.SchedStrategy;
import jline.util.PopulationLattice;
import jline.util.matrix.Matrix;
import org.apache.commons.math3.util.CombinatoricsUtils;

/**
 * Akyildiz-Bolch AMVA method for multi-server BCMP networks.
 */
public final class Pfqn_ab_amva {
    private Pfqn_ab_amva() {}

    public static Ret.pfqnAMVAMS ab_amva(
            Matrix serviceTimes,
            Matrix N,
            Matrix V,
            Matrix nservers,
            List<SchedStrategy> schedStrategies,
            boolean fcfsSchmidt,
            String marginalProbMethod) {
        int M = serviceTimes.getNumRows();
        int K = serviceTimes.getNumCols();

        Ret.pfqnAMVA ret = ab_linearizer(K, M, N, nservers, schedStrategies, V, serviceTimes, fcfsSchmidt, marginalProbMethod);
        return new Ret.pfqnAMVAMS(ret.Q, ret.U, ret.R, ret.C, ret.X, ret.totiter);
    }

    public static Ret.pfqnAMVA ab_linearizer(
            int K, int M, Matrix population, Matrix nservers,
            List<SchedStrategy> type, Matrix v, Matrix s,
            boolean fcfsSchmidt, String marginalProbMethod) {
        Matrix L = new Matrix(M, K);
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < K; r++) {
                L.set(i, r, population.get(r) / M);
            }
        }

        Matrix[] lWithoutR = new Matrix[K];
        for (int i = 0; i < K; i++) lWithoutR[i] = new Matrix(M, K);

        for (int i = 0; i < M; i++) {
            for (int r = 0; r < K; r++) {
                for (int t = 0; t < K; t++) {
                    if (r == t) {
                        lWithoutR[r].set(i, t, (population.get(r) - 1) / M);
                    } else {
                        lWithoutR[r].set(i, t, L.get(i, r));
                    }
                }
            }
        }

        double[][][] D = new double[M][K][K];

        Ret.pfqnAMVA ret = pfqn_ab_core(K, M, population, nservers, type, v, s, 100, D, L, fcfsSchmidt, marginalProbMethod);
        Matrix LUpdated = ret.Q;

        for (int r = 0; r < K; r++) {
            Matrix populationWithoutC = new Matrix(population);
            populationWithoutC.set(r, population.get(r) - 1);

            Matrix lWithoutC = new Matrix(M, K);
            for (int j = 0; j < M; j++) {
                for (int c = 0; c < K; c++) {
                    lWithoutC.set(j, c, lWithoutR[c].get(j, c));
                }
            }

            Ret.pfqnAMVA retWithoutC = pfqn_ab_core(K, M, populationWithoutC, nservers, type, v, s, 100, D, lWithoutC, fcfsSchmidt, marginalProbMethod);

            for (int j = 0; j < M; j++) {
                for (int c = 0; c < K; c++) {
                    lWithoutR[c].set(j, r, retWithoutC.Q.get(j, c));
                }
            }
        }

        for (int i = 0; i < M; i++) {
            for (int r = 0; r < K; r++) {
                double F_ir = LUpdated.get(i, r) / population.get(r);
                for (int t = 0; t < K; t++) {
                    double divisor = (r == t) ? population.get(r) - 1 : population.get(r);
                    double F_irt = lWithoutR[r].get(i, t) / divisor;
                    if (Double.isNaN(F_irt)) F_irt = 0.0;
                    D[i][r][t] = F_irt - F_ir;
                }
            }
        }

        return pfqn_ab_core(K, M, population, nservers, type, v, s, 100, D, LUpdated, fcfsSchmidt, marginalProbMethod);
    }

    public static Ret.pfqnAMVA pfqn_ab_core(
            int K, int M, Matrix population, Matrix nservers,
            List<SchedStrategy> type, Matrix v, Matrix s, int maxiter,
            double[][][] D, Matrix lIn, boolean fcfsSchmidt, String marginalProbMethod) {
        Matrix L = new Matrix(lIn);
        double tol = 1.0 / (4000.0 + (16 * population.sumRows(0)));

        int totalIterations = 0;
        Matrix W = new Matrix(M, K);

        while (totalIterations < maxiter) {
            Matrix F = new Matrix(M, K);
            double[][][] lWithoutJ = new double[M][K][K];

            for (int i = 0; i < M; i++) {
                for (int r = 0; r < K; r++) {
                    if (population.get(r) > 0) {
                        F.set(i, r, L.get(i, r) / population.get(r));
                    } else {
                        F.set(i, r, 0.0);
                    }
                }
            }

            for (int i = 0; i < M; i++) {
                for (int r = 0; r < K; r++) {
                    for (int t = 0; t < K; t++) {
                        double scalar = (r == t) ? population.get(r) - 1 : population.get(r);
                        lWithoutJ[i][r][t] = scalar * (F.get(i, r) + D[i][r][t]);
                    }
                }
            }

            for (int i = 0; i < M; i++) {
                for (int r = 0; r < K; r++) {
                    if (type.get(i) == SchedStrategy.INF) {
                        W.set(i, r, s.get(i, r));
                    } else if (nservers.get(i) == 1.0) {
                        double totalQueueLengthKvecC = 0.0;
                        for (int c = 0; c < K; c++) {
                            totalQueueLengthKvecC += lWithoutJ[i][c][r];
                        }
                        W.set(i, r, s.get(i, r) * (1 + totalQueueLengthKvecC));
                    } else if (fcfsSchmidt && type.get(i) == SchedStrategy.FCFS) {
                        double waitTime = 0.0;
                        Matrix nvec = PopulationLattice.pprod(population);
                        while (nvec != null && nvec.get(0) >= 0) {
                            if (nvec.get(r) > 0) {
                                double Bcn = getBcnForAB(s, i, r, nvec, K, nservers.get(i));
                                double prob = getMarginalProb(oner(nvec, r), oner(population, r), population, lIn.get(i, r), K, i);
                                waitTime += (Bcn * prob);
                            }
                            nvec = PopulationLattice.pprod(nvec, population);
                        }
                        if (waitTime <= 1e-3) waitTime = 0.0;
                        W.set(i, r, waitTime);
                    } else {
                        double queueLength = 0.0;
                        for (int j = 0; j < K; j++) {
                            queueLength += lWithoutJ[i][j][r];
                        }
                        int numServers = (int) nservers.get(i);
                        double multiServerStationWeightedQueueLength = 0.0;
                        if (numServers > 1) {
                            Matrix populationWithoutR = new Matrix(population);
                            populationWithoutR.set(r, populationWithoutR.get(r) - 1);
                            Map<Integer, Double> marginalProbs = findMarginalProbs(queueLength, numServers, populationWithoutR, r, marginalProbMethod);
                            for (int j = 1; j < numServers; j++) {
                                Double mp = marginalProbs.get(j - 1);
                                multiServerStationWeightedQueueLength += (mp != null ? mp : 0.0) * (numServers - j);
                            }
                        }
                        double waitTime = (s.get(i, r) / numServers) * (1 + queueLength + multiServerStationWeightedQueueLength);
                        W.set(i, r, waitTime);
                    }
                }
            }

            Matrix C = Matrix.zeros(1, K);
            for (int r = 0; r < K; r++) {
                double cycleTime = 0.0;
                for (int i = 0; i < M; i++) cycleTime += v.get(i, r) * W.get(i, r);
                C.set(0, r, cycleTime);
            }

            Matrix iterationQueueLength = new Matrix(M, K);
            for (int i = 0; i < M; i++) {
                for (int r = 0; r < K; r++) {
                    double queueLength;
                    if (C.get(r) > 0) {
                        queueLength = population.get(r) * (v.get(i, r) * W.get(i, r) / C.get(r));
                    } else {
                        queueLength = 0.0;
                    }
                    iterationQueueLength.set(i, r, queueLength);
                }
            }

            double maxDifference = 0.0;
            for (int i = 0; i < M; i++) {
                for (int r = 0; r < K; r++) {
                    double difference = Math.abs(L.get(i, r) - iterationQueueLength.get(i, r)) / population.get(r);
                    maxDifference = Math.max(maxDifference, difference);
                }
            }

            totalIterations++;
            L = iterationQueueLength;
            if (maxDifference < tol) break;
        }

        Matrix X = Matrix.zeros(1, K);
        for (int r = 0; r < K; r++) {
            if (W.get(0, r) > 0) {
                X.set(0, r, L.get(0, r) / W.get(0, r));
            } else {
                X.set(0, r, 0.0);
            }
        }

        Matrix U = new Matrix(M, K);
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < K; r++) {
                if (s.get(i, r) > 0) {
                    if (type.get(i) == SchedStrategy.INF) {
                        U.set(i, r, X.get(r) * s.get(i, r));
                    } else {
                        U.set(i, r, (X.get(r) * s.get(i, r)) / nservers.get(i));
                    }
                } else {
                    U.set(i, r, 0.0);
                }
            }
        }

        Matrix CResult = Matrix.zeros(1, K);

        return new Ret.pfqnAMVA(L, U, W, CResult, CResult, X, totalIterations);
    }

    private static Matrix oner(Matrix kvec, int c) {
        Matrix result = new Matrix(kvec);
        result.set(c, result.get(c) - 1);
        return result;
    }

    private static double getBcnForAB(Matrix D, int i, int c, Matrix nvec, int K, double ns) {
        double Bcn = D.get(i, c);
        if (nvec.elementSum() > 1.0) {
            double eps = 1e-12;
            double sum = 0.0;
            for (int t = 0; t < K; t++) {
                sum += nvec.get(t) * D.get(i, t);
            }
            Bcn += (Math.max(0.0, nvec.elementSum() - ns)
                    / Math.max(ns * (nvec.elementSum() - 1), eps) * (sum - D.get(i, c)));
        }
        return Bcn;
    }

    private static double getMarginalProb(Matrix n, Matrix k, Matrix K, double L_jr, int R, int j) {
        double prob = 1.0;
        for (int r = 0; r < R; r++) {
            double frac = L_jr / K.get(r);
            if (frac != 0.0) {
                double term1 = (double) CombinatoricsUtils.binomialCoefficient((int) K.get(r), (int) n.get(r));
                double term2 = Math.pow(frac, n.get(r));
                double term3 = Math.pow(1 - frac, K.get(r) - n.get(r));
                prob *= (term1 * term2 * term3);
            }
        }
        return prob;
    }

    private static Matrix weightFun(Matrix population, double alpha, double beta) {
        int maxClassPopulation = 0;
        for (int i = 0; i < population.length(); i++) {
            if (population.get(i) > maxClassPopulation) {
                maxClassPopulation = (int) population.get(i);
            }
        }

        double[] scalingFun = new double[maxClassPopulation + 1];
        if (maxClassPopulation >= 1) {
            scalingFun[1] = alpha;
            for (int n = 2; n <= maxClassPopulation; n++) {
                scalingFun[n] = beta * scalingFun[n - 1];
            }
        }

        double[][] weightFun = new double[maxClassPopulation + 1][maxClassPopulation + 1];
        weightFun[0][0] = 1.0;

        for (int l = 1; l <= maxClassPopulation; l++) {
            for (int j = 0; j < l; j++) {
                weightFun[l][j] = weightFun[l - 1][j] - (weightFun[l - 1][j] * scalingFun[l]) / 100.0;
            }
            double sum = 0.0;
            for (int j = 0; j < l; j++) sum += weightFun[l][j];
            weightFun[l][l] = 1 - sum;
        }

        return new Matrix(weightFun);
    }

    public static Map<Integer, Double> findMarginalProbs(
            double avgJobs, int numServers, Matrix population, int classIdx, String marginalProbMethod) {
        Map<Integer, Double> marginalProbs = new HashMap<Integer, Double>();

        if ("scat".equals(marginalProbMethod)) {
            int floorVal = (int) Math.floor(avgJobs);
            int ceilVal = floorVal + 1;
            marginalProbs.put(floorVal, (double) ceilVal - avgJobs);
            marginalProbs.put(ceilVal, avgJobs - floorVal);
            return marginalProbs;
        }

        double ALPHA = 45.0;
        double BETA = 0.7;

        Matrix w = weightFun(population, ALPHA, BETA);

        int floorVal = (int) Math.floor(avgJobs);
        int ceiling = floorVal + 1;
        int maxVal = Math.min((2 * floorVal) + 1, numServers - 2);

        for (int j = 0; j <= maxVal; j++) {
            double prob;
            if (j <= floorVal) {
                int lDist = floorVal - j;
                int lowerVal = floorVal - lDist;
                int upperVal = ceiling + lDist;
                if (lDist > 25) {
                    prob = 0.0;
                } else {
                    if (floorVal < population.get(classIdx)) {
                        prob = w.get(floorVal, lDist) * (((double) upperVal - avgJobs) / (upperVal - lowerVal));
                    } else {
                        prob = 0.0;
                    }
                }
                marginalProbs.put(j, prob);
            } else {
                int uDist = j - ceiling;
                if (uDist > 25) {
                    marginalProbs.put(j, 0.0);
                } else if (j > population.get(classIdx) - 1 && uDist < 25) {
                    Double existingProb = marginalProbs.get((int) (population.get(classIdx) - 1));
                    if (existingProb == null) existingProb = 0.0;
                    Double mp = marginalProbs.get(floorVal - uDist);
                    double newProb = existingProb + (w.get(floorVal, uDist) - (mp == null ? 0.0 : mp));
                    marginalProbs.put((int) (population.get(classIdx) - 1), newProb);
                } else {
                    Double mp = marginalProbs.get(floorVal - uDist);
                    double newProb = w.get(floorVal, uDist) - (mp == null ? 0.0 : mp);
                    marginalProbs.put(j, newProb);
                }
            }
        }

        return marginalProbs;
    }
}
