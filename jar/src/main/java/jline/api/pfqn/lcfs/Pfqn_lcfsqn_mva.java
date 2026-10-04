/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.pfqn.lcfs;

import java.util.ArrayList;
import java.util.List;

import jline.util.PopulationLattice;
import jline.util.matrix.Matrix;

/**
 * Mean Value Analysis for LCFS Queueing Networks.
 *
 * Implements exact MVA for 2-station LCFS queueing networks with
 * LCFS and LCFS-PR scheduling disciplines using log-space arithmetic
 * for numerical stability.
 */
public final class Pfqn_lcfsqn_mva {
    private Pfqn_lcfsqn_mva() {}

    public static LcfsqnMvaResult pfqn_lcfsqn_mva(Matrix alpha, Matrix beta) {
        return pfqn_lcfsqn_mva(alpha, beta, null);
    }

    public static LcfsqnMvaResult pfqn_lcfsqn_mva(Matrix alpha, Matrix beta, Matrix N) {
        int R = alpha.length();
        Matrix populationVector = (N != null) ? N : Matrix.ones(1, R);
        int K = (int) populationVector.elementSum();

        if (K == 0) {
            return new LcfsqnMvaResult(new Matrix(1, R), new Matrix(2, R), new Matrix(2, R), new Matrix(2, R));
        }

        double[] logAlpha = new double[R];
        double[] logBeta = new double[R];
        for (int i = 0; i < R; i++) {
            logAlpha[i] = Math.log(alpha.get(i));
            logBeta[i] = Math.log(beta.get(i));
        }

        Matrix prods = new Matrix(1, R);
        for (int r = 0; r < R; r++) {
            double prod = 1.0;
            for (int j = 0; j < r; j++) {
                prod *= (populationVector.get(j) + 1);
            }
            prods.set(r, prod);
        }

        int totalSize = 1;
        for (int r = 0; r < R; r++) {
            totalSize *= ((int) populationVector.get(r) + 1);
        }

        Matrix[] QN = new Matrix[totalSize];
        Matrix[] TN = new Matrix[totalSize];
        Matrix[] BN_scaled = new Matrix[totalSize];
        double[][][] log_scale_BN = new double[totalSize][2][R];
        for (int i = 0; i < totalSize; i++) {
            QN[i] = new Matrix(2, R);
            TN[i] = new Matrix(1, R);
            BN_scaled[i] = new Matrix(2, R);
            for (int s = 0; s < 2; s++) {
                for (int r = 0; r < R; r++) log_scale_BN[i][s][r] = Double.NEGATIVE_INFINITY;
            }
        }

        Matrix n = PopulationLattice.pprod(populationVector);
        while (n != null && n.get(0) >= 0) {
            int idx = hashpopLcfs(n, populationVector, R, prods);
            int sumN = (int) n.elementSum();

            if (sumN == 1) {
                for (int r = 0; r < R; r++) {
                    if (n.get(r) > 0) {
                        double denom = alpha.get(r) + beta.get(r);
                        QN[idx].set(0, r, alpha.get(r) / denom);
                        QN[idx].set(1, r, beta.get(r) / denom);
                        BN_scaled[idx].set(0, r, alpha.get(r) / denom);
                        BN_scaled[idx].set(1, r, beta.get(r) / denom);
                        log_scale_BN[idx][0][r] = 0.0;
                        log_scale_BN[idx][1][r] = 0.0;
                        TN[idx].set(0, r, 1.0 / denom);
                    }
                }
            } else if (sumN > 1) {
                double log_scale_np = 0.0;
                for (int r = 0; r < R; r++) {
                    log_scale_np += n.get(r) * logAlpha[r];
                }

                for (int k = 0; k < R; k++) {
                    if (n.get(k) > 0) {
                        Matrix n_minus_k = oner(n, k);
                        int idx_k = hashpopLcfs(n_minus_k, populationVector, R, prods);

                        double Wnp_unscaled = 1.0 + QN[idx_k].get(0, k);
                        double Wpr_unscaled = 1.0 + QN[idx_k].get(1, k);

                        for (int r = 0; r < R; r++) {
                            if (r != k && n.get(r) > 0) {
                                Matrix n_minus_r = oner(n, r);
                                int idx_r = hashpopLcfs(n_minus_r, populationVector, R, prods);

                                if (BN_scaled[idx_k].get(0, r) > 0 && BN_scaled[idx_r].get(0, k) > 0) {
                                    double log_ratio_np = Math.log(BN_scaled[idx_k].get(0, r)) + log_scale_BN[idx_k][0][r]
                                            - Math.log(BN_scaled[idx_r].get(0, k)) - log_scale_BN[idx_r][0][k];
                                    double ratio_np = Math.exp(log_ratio_np);
                                    Wnp_unscaled += (alpha.get(k) / alpha.get(r)) * ratio_np * QN[idx_r].get(0, k);
                                }
                                if (BN_scaled[idx_k].get(1, r) > 0 && BN_scaled[idx_r].get(1, k) > 0) {
                                    double log_ratio_pr = Math.log(BN_scaled[idx_k].get(1, r)) + log_scale_BN[idx_k][1][r]
                                            - Math.log(BN_scaled[idx_r].get(1, k)) - log_scale_BN[idx_r][1][k];
                                    double ratio_pr = Math.exp(log_ratio_pr);
                                    Wpr_unscaled += (alpha.get(r) / alpha.get(k)) * ratio_pr * QN[idx_r].get(1, k);
                                }
                            }
                        }

                        double log_scale_pr_k = (sumN - 1) * logAlpha[k] + logBeta[k];
                        double log_Wnp = log_scale_np + Math.log(Wnp_unscaled);
                        double log_Wpr = log_scale_pr_k + Math.log(Wpr_unscaled);
                        double max_log = Math.max(log_Wnp, log_Wpr);
                        double log_sum_W = max_log + Math.log(Math.exp(log_Wnp - max_log) + Math.exp(log_Wpr - max_log));

                        BN_scaled[idx].set(0, k, n.get(k));
                        log_scale_BN[idx][0][k] = log_scale_np - log_sum_W;
                        BN_scaled[idx].set(1, k, n.get(k));
                        log_scale_BN[idx][1][k] = log_scale_pr_k - log_sum_W;
                    }
                }

                for (int k = 0; k < R; k++) {
                    if (n.get(k) > 0) {
                        for (int station = 0; station < 2; station++) {
                            List<Double> terms = new ArrayList<Double>();
                            List<Double> logScales = new ArrayList<Double>();
                            if (BN_scaled[idx].get(station, k) > 0) {
                                terms.add(BN_scaled[idx].get(station, k));
                                logScales.add(log_scale_BN[idx][station][k]);
                            }
                            for (int r = 0; r < R; r++) {
                                if (n.get(r) > 0) {
                                    Matrix n_minus_r = oner(n, r);
                                    int idx_r = hashpopLcfs(n_minus_r, populationVector, R, prods);
                                    if (BN_scaled[idx].get(station, r) > 0 && QN[idx_r].get(station, k) > 0) {
                                        terms.add(BN_scaled[idx].get(station, r) * QN[idx_r].get(station, k));
                                        logScales.add(log_scale_BN[idx][station][r]);
                                    }
                                }
                            }
                            QN[idx].set(station, k, logsumexp(terms, logScales));
                        }
                    }
                }

                for (int k = 0; k < R; k++) {
                    if (n.get(k) > 0) {
                        Matrix n_minus_k = oner(n, k);
                        int idx_k = hashpopLcfs(n_minus_k, populationVector, R, prods);
                        double Unp = 0.0;
                        for (int r = 0; r < R; r++) {
                            Unp += alpha.get(r) * TN[idx_k].get(0, r);
                        }
                        List<Double> terms = new ArrayList<Double>();
                        List<Double> logScales = new ArrayList<Double>();
                        for (int r = 0; r < R; r++) {
                            if (n.get(r) > 0) {
                                Matrix n_minus_r = oner(n, r);
                                int idx_r = hashpopLcfs(n_minus_r, populationVector, R, prods);
                                if (BN_scaled[idx].get(0, r) > 0 && TN[idx_r].get(0, k) > 0) {
                                    terms.add(BN_scaled[idx].get(0, r) * TN[idx_r].get(0, k));
                                    logScales.add(log_scale_BN[idx][0][r]);
                                }
                            }
                        }
                        if (BN_scaled[idx].get(0, k) > 0 && (1 - Unp) > 0) {
                            terms.add((1.0 / alpha.get(k)) * BN_scaled[idx].get(0, k) * (1 - Unp));
                            logScales.add(log_scale_BN[idx][0][k]);
                        }
                        TN[idx].set(0, k, logsumexp(terms, logScales));
                    }
                }
            }
            n = PopulationLattice.pprod(n, populationVector);
        }

        int finalIdx = totalSize - 1;
        Matrix Q = QN[finalIdx];
        Matrix T = TN[finalIdx];
        Matrix B = new Matrix(2, R);
        for (int r = 0; r < R; r++) {
            if (BN_scaled[finalIdx].get(0, r) > 0) {
                B.set(0, r, BN_scaled[finalIdx].get(0, r) * Math.exp(log_scale_BN[finalIdx][0][r]));
            }
            if (BN_scaled[finalIdx].get(1, r) > 0) {
                B.set(1, r, BN_scaled[finalIdx].get(1, r) * Math.exp(log_scale_BN[finalIdx][1][r]));
            }
        }

        Matrix U = new Matrix(2, R);
        for (int r = 0; r < R; r++) {
            U.set(0, r, T.get(0, r) * alpha.get(r));
            U.set(1, r, T.get(0, r) * beta.get(r));
        }

        return new LcfsqnMvaResult(T, Q, U, B);
    }

    private static int hashpopLcfs(Matrix n, Matrix N, int R, Matrix prods) {
        int idx = 0;
        for (int r = 0; r < R; r++) {
            idx += (int) (prods.get(r) * n.get(r));
        }
        return idx;
    }

    private static Matrix oner(Matrix n, int r) {
        Matrix result = n.copy();
        result.set(r, result.get(r) - 1);
        return result;
    }

    private static double logsumexp(List<Double> terms, List<Double> logScales) {
        if (terms.isEmpty()) return 0.0;
        List<Double> validLog = new ArrayList<Double>();
        for (int i = 0; i < terms.size(); i++) {
            double t = terms.get(i);
            double s = logScales.get(i);
            if (t > 0 && Double.isFinite(s)) {
                validLog.add(Math.log(t) + s);
            }
        }
        if (validLog.isEmpty()) return 0.0;
        double maxLog = Double.NEGATIVE_INFINITY;
        for (Double v : validLog) if (v > maxLog) maxLog = v;
        if (!Double.isFinite(maxLog)) return 0.0;
        double sum = 0.0;
        for (Double v : validLog) sum += Math.exp(v - maxLog);
        return Math.exp(maxLog) * sum;
    }
}
