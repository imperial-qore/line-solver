/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.pfqn.mva;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.util.CombinatoricsUtils;

import jline.io.Ret;
import jline.lang.constant.SchedStrategy;
import jline.util.PopulationLattice;
import jline.util.matrix.Matrix;

/**
 * Schmidt MVA algorithm for multi-class FCFS queueing networks.
 */
public final class Pfqn_schmidt_amva {

    private Pfqn_schmidt_amva() {}

    public static Ret.pfqnAMVASchmidt pfqn_schmidt(
            Matrix rates,
            Matrix N,
            Matrix S,
            Matrix v,
            List<SchedStrategy> sched) {
        int M = rates.getNumRows();
        int R = rates.getNumCols();

        // Convert rates to service demands (D = 1/rate)
        Matrix D = new Matrix(M, R);
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                if (rates.get(i, r) == 0.0) {
                    D.set(i, r, 0.0);
                } else {
                    D.set(i, r, 1.0 / rates.get(i, r));
                }
            }
        }

        int[] closedClasses = new int[R];
        for (int i = 0; i < R; i++) closedClasses[i] = i;

        Matrix XN = Matrix.zeros(1, R);
        Matrix UN = Matrix.zeros(M, R);
        Matrix CN = Matrix.zeros(M, R);
        Matrix QN = Matrix.zeros(M, R);
        List<Map<Object, Double>> PN = new ArrayList<Map<Object, Double>>();

        for (int i = 0; i < M; i++) {
            PN.add(new HashMap<Object, Double>());
        }

        int C = closedClasses.length;
        Matrix Nc = N;

        Matrix prods = Matrix.zeros(1, C);
        for (int r = 0; r < C; r++) {
            double prod = 1.0;
            for (int i = 0; i < r; i++) {
                prod *= (Nc.get(i) + 1);
            }
            prods.set(0, r, prod);
        }

        Matrix kvec = PopulationLattice.pprod(Nc);
        int totalStates = (int) prod(Nc);

        List<Matrix> L = new ArrayList<Matrix>();
        List<Matrix> Pc = new ArrayList<Matrix>();

        for (int i = 0; i < M; i++) {
            L.add(Matrix.zeros(R, totalStates));
            SchedStrategy ss = sched.get(i);
            if (ss == SchedStrategy.INF) {
                Pc.add(null);
            } else if (ss == SchedStrategy.PS) {
                boolean singleServerPS = S.get(i) == 1.0;
                if (singleServerPS) {
                    Pc.add(null);
                } else {
                    Pc.add(Matrix.zeros((int) S.get(i), totalStates));
                }
            } else if (ss == SchedStrategy.FCFS) {
                boolean classIndependent = true;
                for (int r = 1; r < R; r++) {
                    if (D.get(i, r) != D.get(i, 0)) {
                        classIndependent = false;
                        break;
                    }
                }
                boolean isSingleServer = S.get(i) == 1.0;

                if (classIndependent) {
                    if (isSingleServer) {
                        Pc.add(null);
                    } else {
                        Pc.add(Matrix.zeros((int) S.get(i), totalStates));
                    }
                } else {
                    Pc.add(Matrix.zeros(totalStates, totalStates));
                }
            } else {
                Pc.add(null);
            }
        }

        Matrix[] x = new Matrix[M];
        Matrix[] w = new Matrix[M];
        for (int i = 0; i < M; i++) {
            x[i] = Matrix.zeros(C, totalStates);
            w[i] = Matrix.zeros(C, totalStates);
        }

        for (int ist = 0; ist < M; ist++) {
            if (Pc.get(ist) != null) {
                Pc.get(ist).set(0, PopulationLattice.hashpop(kvec, Nc), 1.0);
            }
        }

        kvec = PopulationLattice.pprod(kvec, Nc);
        int hkvec = PopulationLattice.hashpop(kvec, Nc);

        while (allGE(kvec, 0) && allLE(kvec, Nc)) {
            Matrix kprods = Matrix.zeros(1, C);
            for (int r = 0; r < C; r++) {
                double prod = 1.0;
                for (int i = 0; i < r; i++) {
                    prod *= (kvec.get(i) + 1);
                }
                kprods.set(r, prod);
            }

            for (int i = 0; i < M; i++) {
                for (int c = 0; c < C; c++) {
                    double ns = S.get(i);
                    hkvec = PopulationLattice.hashpop(kvec, Nc);

                    Matrix kvec_c = oner(kvec, c);
                    int hkvec_c = PopulationLattice.hashpop(kvec_c, Nc);

                    if (kvec.get(c) == 0.0) continue;

                    SchedStrategy ss = sched.get(i);
                    if (ss == SchedStrategy.INF) {
                        w[i].set(c, hkvec, D.get(i, c));
                    } else if (ss == SchedStrategy.FCFS) {
                        boolean classIndependent = true;
                        for (int r = 1; r < R; r++) {
                            if (D.get(i, r) != D.get(i, 0)) {
                                classIndependent = false;
                                break;
                            }
                        }

                        if (!classIndependent) {
                            if (ns == 1.0) {
                                double totalQueueLengthKvecC = 0.0;
                                for (int r = 0; r < R; r++) {
                                    totalQueueLengthKvecC += L.get(i).get(r, hkvec_c);
                                }
                                w[i].set(c, hkvec, D.get(i, c) * (1 + totalQueueLengthKvecC));
                            } else {
                                double waitTime = 0.0;
                                Matrix nvec = PopulationLattice.pprod(kvec);
                                while (nvec != null && nvec.value() >= 0) {
                                    if (nvec.get(c) > 0) {
                                        Matrix nvec_c = oner(nvec, c);
                                        int hnvec_c = PopulationLattice.hashpop(nvec_c, Nc, C, prods);

                                        double Bcn = getBcn(D, i, c, nvec, C, ns);
                                        double prob = Pc.get(i).get(hnvec_c, hkvec_c);
                                        waitTime += (Bcn * prob);
                                    }
                                    nvec = PopulationLattice.pprod(nvec, kvec);
                                }
                                w[i].set(c, hkvec, waitTime);
                            }
                        }
                        if (classIndependent) {
                            if (ns == 1.0) {
                                double totalQueueLengthKvecC = 0.0;
                                for (int r = 0; r < R; r++) {
                                    totalQueueLengthKvecC += L.get(i).get(r, hkvec_c);
                                }
                                w[i].set(c, hkvec, D.get(i, c) * (1 + totalQueueLengthKvecC));
                            } else {
                                double totalQueueLengthKvecC = 0.0;
                                for (int r = 0; r < R; r++) {
                                    totalQueueLengthKvecC += L.get(i).get(r, hkvec_c);
                                }
                                double weightedQueueLength = 0.0;
                                int nsm1 = (int) (ns - 1);
                                for (int j = 0; j < nsm1; j++) {
                                    weightedQueueLength += ((ns - j - 1) * Pc.get(i).get(j, hkvec_c));
                                }
                                double totalWaitTime = (D.get(i, c) / ns) * (1 + totalQueueLengthKvecC + weightedQueueLength);
                                w[i].set(c, hkvec, totalWaitTime);
                            }
                        }
                    } else if (ss == SchedStrategy.PS) {
                        if (ns == 1.0) {
                            double totalQueueLengthKvecC = 0.0;
                            for (int r = 0; r < R; r++) {
                                totalQueueLengthKvecC += L.get(i).get(r, hkvec_c);
                            }
                            w[i].set(c, hkvec, D.get(i, c) * (1 + totalQueueLengthKvecC));
                        } else {
                            double totalQueueLengthKvecC = 0.0;
                            for (int r = 0; r < R; r++) {
                                totalQueueLengthKvecC += L.get(i).get(r, hkvec_c);
                            }
                            double weightedQueueLength = 0.0;
                            int nsm1 = (int) (ns - 1);
                            for (int j = 0; j < nsm1; j++) {
                                weightedQueueLength += ((ns - j - 1) * Pc.get(i).get(j, hkvec_c));
                            }
                            double totalWaitTime = (D.get(i, c) / ns) * (1 + totalQueueLengthKvecC + weightedQueueLength);
                            w[i].set(c, hkvec, totalWaitTime);
                        }
                    }
                }
            }

            // Compute throughputs
            for (int c = 0; c < C; c++) {
                double denom = 0.0;
                for (int i = 0; i < M; i++) {
                    denom += (v.get(i, c) * w[i].get(c, hkvec));
                }

                for (int i = 0; i < M; i++) {
                    if (denom > 0) {
                        x[i].set(c, hkvec, v.get(i, c) * kvec.get(c) / denom);
                    } else {
                        x[i].set(c, hkvec, 0.0);
                    }
                }
            }

            // Update queue lengths and state probabilities
            for (int i = 0; i < M; i++) {
                for (int c = 0; c < C; c++) {
                    L.get(i).set(c, hkvec, x[i].get(c, hkvec) * w[i].get(c, hkvec));
                }
                int ns = (int) S.get(i);
                int K_j = 0;
                for (int r = 0; r < R; r++) {
                    if (v.get(i, r) > 0) {
                        K_j += (int) Nc.get(r);
                    }
                }

                SchedStrategy ss = sched.get(i);
                if (ss == SchedStrategy.FCFS) {
                    boolean classIndependent = true;
                    for (int r = 1; r < R; r++) {
                        if (D.get(i, r) != D.get(i, 0)) {
                            classIndependent = false;
                            break;
                        }
                    }

                    if (!classIndependent) {
                        Matrix nvec = PopulationLattice.pprod(kvec);
                        nvec = PopulationLattice.pprod(nvec, kvec);
                        double sumOfAllProbs = 0.0;
                        while (nvec != null && nvec.value() >= 0) {
                            int hnvec = PopulationLattice.hashpop(nvec, Nc, C, prods);

                            double prob = 0.0;
                            for (int r = 0; r < C; r++) {
                                if (nvec.get(r) > 0) {
                                    int hnvec_c = PopulationLattice.hashpop(oner(nvec, r), Nc, C, prods);
                                    int hkvec_c = PopulationLattice.hashpop(oner(kvec, r), Nc);
                                    double Bcn = getBcn(D, i, r, nvec, C, (double) ns);
                                    double capacity_inv = 1 / nvec.elementSum();
                                    double x_ir = x[i].get(r, hkvec);
                                    double prob_c = Pc.get(i).get(hnvec_c, hkvec_c);
                                    double classProb = (Bcn * capacity_inv * x_ir * prob_c);
                                    prob += classProb;
                                }
                            }
                            Pc.get(i).set(hnvec, hkvec, prob);
                            sumOfAllProbs += prob;
                            nvec = PopulationLattice.pprod(nvec, kvec);
                        }
                        Pc.get(i).set(0, hkvec, Math.max(1e-12, 1.0 - sumOfAllProbs));
                    } else if (ns > 1) {
                        double meanQueueLength = 0.0;
                        for (int r = 0; r < R; r++) {
                            meanQueueLength += L.get(i).get(r, hkvec);
                        }
                        for (int n = 1; n < ns; n++) {
                            double x1 = (double) CombinatoricsUtils.binomialCoefficient(K_j, n);
                            double frac = meanQueueLength / K_j;
                            double x2 = Math.pow(frac, (double) n);
                            double x3 = Math.pow(1 - frac, (double) (K_j - n));
                            double prob = x1 * x2 * x3;
                            Pc.get(i).set(n, hkvec, prob);
                        }
                        double sum1 = 0.0;
                        double sum2 = 0.0;
                        for (int r = 0; r < R; r++) {
                            sum1 += D.get(i, r) * x[i].get(r, hkvec);
                        }
                        for (int n = 0; n < ns; n++) {
                            sum2 += (ns - n) * Pc.get(i).get(n, hkvec);
                        }
                        double term = (sum1 + sum2) / ns;
                        Pc.get(i).set(0, hkvec, Math.max(1e-12, 1 - term));
                    }
                } else if (ss == SchedStrategy.PS) {
                    if (ns > 1) {
                        double meanQueueLength = 0.0;
                        for (int r = 0; r < R; r++) {
                            meanQueueLength += L.get(i).get(r, hkvec);
                        }
                        for (int n = 1; n < ns; n++) {
                            double x1 = (double) CombinatoricsUtils.binomialCoefficient(K_j, n);
                            double frac = meanQueueLength / K_j;
                            double x2 = Math.pow(frac, (double) n);
                            double x3 = Math.pow(1 - frac, (double) (K_j - n));
                            double prob = x1 * x2 * x3;
                            Pc.get(i).set(n, hkvec, prob);
                        }
                        double sum1 = 0.0;
                        double sum2 = 0.0;
                        for (int r = 0; r < R; r++) {
                            sum1 += D.get(i, r) * x[i].get(r, hkvec);
                        }
                        for (int n = 0; n < ns; n++) {
                            sum2 += (ns - n) * Pc.get(i).get(n, hkvec);
                        }
                        double term = (sum1 + sum2) / ns;
                        Pc.get(i).set(0, hkvec, Math.max(1e-12, 1 - term));
                    }
                }
            }

            kvec = PopulationLattice.pprod(kvec, Nc);
        }

        int hkvecFinal = hkvec;

        // Extract final results
        for (int c = 0; c < C; c++) {
            double totalResponseTime = 0.0;
            for (int i = 0; i < M; i++) {
                totalResponseTime += w[i].get(c, hkvecFinal);
            }
            XN.set(c, Nc.get(c) / totalResponseTime);
        }

        for (int m = 0; m < M; m++) {
            for (int c = 0; c < C; c++) {
                UN.set(m, c, (D.get(m, c) * XN.get(c)) / S.get(m));
            }
        }

        for (int m = 0; m < M; m++) {
            for (int c = 0; c < C; c++) {
                CN.set(m, c, w[m].get(c, hkvecFinal));
            }
        }

        for (int m = 0; m < M; m++) {
            Matrix Lmat = L.get(m);
            for (int c = 0; c < C; c++) {
                QN.set(m, c, Lmat.get(c, hkvecFinal));
            }
        }

        return new Ret.pfqnAMVASchmidt(XN, QN, UN, CN, PN);
    }

    public static Ret.pfqnAMVASchmidt pfqn_schmidt_ext(
            Matrix rates,
            Matrix N,
            Matrix S,
            Matrix v,
            List<SchedStrategy> sched) {
        int M = rates.getNumRows();
        int R = rates.getNumCols();

        Matrix D = new Matrix(M, R);
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                if (rates.get(i, r) == 0.0 || Double.isNaN(rates.get(i, r))) {
                    D.set(i, r, 0.0);
                } else {
                    D.set(i, r, 1.0 / rates.get(i, r));
                }
            }
        }

        // see _kb/03-api-layer.md for rationale
        Map<String, Matrix> alphas = new HashMap<String, Matrix>();

        for (int i = 0; i < M; i++) {
            if (sched.get(i) == SchedStrategy.FCFS) {
                for (int r = 0; r < R; r++) {
                    Matrix rates_mod = new Matrix(M, R + 1);
                    Matrix N_mod = new Matrix(1, R + 1);
                    Matrix v_mod = new Matrix(M, R + 1);

                    for (int k = 0; k < R; k++) {
                        N_mod.set(k, (k == r) ? N.get(k) - 1 : N.get(k));
                        for (int j = 0; j < M; j++) {
                            rates_mod.set(j, k, rates.get(j, k));
                            v_mod.set(j, k, v.get(j, k));
                            rates_mod.set(j, R, (i == j) ? rates.get(j, r) : 0.0);
                            v_mod.set(j, R, (i == j) ? 1.0 : 0.0);
                        }
                    }
                    N_mod.set(0, R, 1.0);

                    Ret.pfqnAMVASchmidt result = pfqn_schmidt(rates_mod, N_mod, S, v_mod, sched);
                    alphas.put(i + ":" + r, result.U);
                }
            }
        }

        int[] closedClasses = new int[R];
        for (int i = 0; i < R; i++) closedClasses[i] = i;

        Matrix XN = Matrix.zeros(1, R);
        Matrix UN = Matrix.zeros(M, R);
        Matrix CN = Matrix.zeros(M, R);
        Matrix QN = Matrix.zeros(M, R);
        List<Map<Object, Double>> PN = new ArrayList<Map<Object, Double>>();

        for (int i = 0; i < M; i++) {
            PN.add(new HashMap<Object, Double>());
        }

        int C = closedClasses.length;
        Matrix Nc = N;

        Matrix prods = Matrix.zeros(1, C);
        for (int r = 0; r < C; r++) {
            double prod = 1.0;
            for (int i = 0; i < r; i++) {
                prod *= (Nc.get(i) + 1);
            }
            prods.set(0, r, prod);
        }

        Matrix kvec = PopulationLattice.pprod(Nc);
        int totalStates = (int) prod(Nc);

        List<Matrix> L = new ArrayList<Matrix>();
        List<Matrix> Pc = new ArrayList<Matrix>();

        for (int i = 0; i < M; i++) {
            L.add(Matrix.zeros(R, totalStates));
            SchedStrategy ss = sched.get(i);
            if (ss == SchedStrategy.INF) {
                Pc.add(null);
            } else if (ss == SchedStrategy.PS) {
                boolean singleServerPS = S.get(i) == 1.0;
                if (singleServerPS) {
                    Pc.add(null);
                } else {
                    Pc.add(Matrix.zeros((int) S.get(i), totalStates));
                }
            } else if (ss == SchedStrategy.FCFS) {
                boolean classIndependent = true;
                for (int r = 1; r < R; r++) {
                    if (D.get(i, r) != D.get(i, 0)) {
                        classIndependent = false;
                        break;
                    }
                }
                boolean isSingleServer = S.get(i) == 1.0;
                if (classIndependent) {
                    if (isSingleServer) {
                        Pc.add(null);
                    } else {
                        Pc.add(Matrix.zeros((int) S.get(i), totalStates));
                    }
                } else {
                    Pc.add(Matrix.zeros(totalStates, totalStates));
                }
            } else {
                Pc.add(null);
            }
        }

        Matrix[] x = new Matrix[M];
        Matrix[] w = new Matrix[M];
        for (int i = 0; i < M; i++) {
            x[i] = Matrix.zeros(C, totalStates);
            w[i] = Matrix.zeros(C, totalStates);
        }

        for (int ist = 0; ist < M; ist++) {
            if (Pc.get(ist) != null) {
                Pc.get(ist).set(0, PopulationLattice.hashpop(kvec, Nc), 1.0);
            }
        }

        kvec = PopulationLattice.pprod(kvec, Nc);
        int hkvec = PopulationLattice.hashpop(kvec, Nc);

        while (allGE(kvec, 0) && allLE(kvec, Nc)) {
            Matrix kprods = Matrix.zeros(1, C);
            for (int r = 0; r < C; r++) {
                double prod = 1.0;
                for (int i = 0; i < r; i++) {
                    prod *= (kvec.get(i) + 1);
                }
                kprods.set(r, prod);
            }

            for (int i = 0; i < M; i++) {
                for (int c = 0; c < C; c++) {
                    double ns = S.get(i);
                    hkvec = PopulationLattice.hashpop(kvec, Nc);

                    Matrix kvec_c = oner(kvec, c);
                    int hkvec_c = PopulationLattice.hashpop(kvec_c, Nc);

                    if (kvec.get(c) == 0.0) continue;

                    SchedStrategy ss = sched.get(i);
                    if (ss == SchedStrategy.INF) {
                        w[i].set(c, hkvec, D.get(i, c));
                    } else if (ss == SchedStrategy.FCFS) {
                        boolean classIndependent = true;
                        for (int r = 1; r < R; r++) {
                            if (D.get(i, r) != D.get(i, 0)) {
                                classIndependent = false;
                                break;
                            }
                        }

                        if (!classIndependent) {
                            if (ns == 1.0) {
                                double totalQueueLengthKvecC = 0.0;
                                for (int r = 0; r < R; r++) {
                                    totalQueueLengthKvecC += L.get(i).get(r, hkvec_c);
                                }
                                w[i].set(c, hkvec, D.get(i, c) * (1 + totalQueueLengthKvecC));
                            } else {
                                double waitTime = 0.0;
                                Matrix nvec = PopulationLattice.pprod(kvec);
                                while (nvec != null && nvec.value() >= 0) {
                                    if (nvec.get(c) > 0) {
                                        Matrix nvec_c = oner(nvec, c);
                                        int hnvec_c = PopulationLattice.hashpop(nvec_c, Nc, C, prods);

                                        double Bcn;
                                        if (Nc.get(c) > 1.0) {
                                            Bcn = getBcnExt(alphas.get(i + ":" + c), D, i, c, nvec, C, ns);
                                        } else {
                                            Bcn = getBcn(D, i, c, nvec, C, ns);
                                        }

                                        double prob = Pc.get(i).get(hnvec_c, hkvec_c);
                                        waitTime += (Bcn * prob);
                                    }
                                    nvec = PopulationLattice.pprod(nvec, kvec);
                                }
                                w[i].set(c, hkvec, waitTime);
                            }
                        }
                        if (classIndependent) {
                            if (ns == 1.0) {
                                double totalQueueLengthKvecC = 0.0;
                                for (int r = 0; r < R; r++) {
                                    totalQueueLengthKvecC += L.get(i).get(r, hkvec_c);
                                }
                                w[i].set(c, hkvec, D.get(i, c) * (1 + totalQueueLengthKvecC));
                            } else {
                                double totalQueueLengthKvecC = 0.0;
                                for (int r = 0; r < R; r++) {
                                    totalQueueLengthKvecC += L.get(i).get(r, hkvec_c);
                                }
                                double weightedQueueLength = 0.0;
                                int nsm1 = (int) (ns - 1);
                                for (int j = 0; j < nsm1; j++) {
                                    weightedQueueLength += ((ns - j - 1) * Pc.get(i).get(j, hkvec_c));
                                }
                                double totalWaitTime = (D.get(i, c) / ns) * (1 + totalQueueLengthKvecC + weightedQueueLength);
                                w[i].set(c, hkvec, totalWaitTime);
                            }
                        }
                    } else if (ss == SchedStrategy.PS) {
                        if (ns == 1.0) {
                            double totalQueueLengthKvecC = 0.0;
                            for (int r = 0; r < R; r++) {
                                totalQueueLengthKvecC += L.get(i).get(r, hkvec_c);
                            }
                            w[i].set(c, hkvec, D.get(i, c) * (1 + totalQueueLengthKvecC));
                        } else {
                            double totalQueueLengthKvecC = 0.0;
                            for (int r = 0; r < R; r++) {
                                totalQueueLengthKvecC += L.get(i).get(r, hkvec_c);
                            }
                            double weightedQueueLength = 0.0;
                            int nsm1 = (int) (ns - 1);
                            for (int j = 0; j < nsm1; j++) {
                                weightedQueueLength += ((ns - j - 1) * Pc.get(i).get(j, hkvec_c));
                            }
                            double totalWaitTime = (D.get(i, c) / ns) * (1 + totalQueueLengthKvecC + weightedQueueLength);
                            w[i].set(c, hkvec, totalWaitTime);
                        }
                    }
                }
            }

            for (int c = 0; c < C; c++) {
                double denom = 0.0;
                for (int i = 0; i < M; i++) {
                    denom += (v.get(i, c) * w[i].get(c, hkvec));
                }

                for (int i = 0; i < M; i++) {
                    if (denom > 0) {
                        x[i].set(c, hkvec, v.get(i, c) * kvec.get(c) / denom);
                    } else {
                        x[i].set(c, hkvec, 0.0);
                    }
                }
            }

            for (int i = 0; i < M; i++) {
                for (int c = 0; c < C; c++) {
                    L.get(i).set(c, hkvec, x[i].get(c, hkvec) * w[i].get(c, hkvec));
                }
                int ns = (int) S.get(i);
                int K_j = 0;
                for (int r = 0; r < R; r++) {
                    if (v.get(i, r) > 0) {
                        K_j += (int) Nc.get(r);
                    }
                }

                SchedStrategy ss = sched.get(i);
                if (ss == SchedStrategy.FCFS) {
                    boolean classIndependent = true;
                    for (int r = 1; r < R; r++) {
                        if (D.get(i, r) != D.get(i, 0)) {
                            classIndependent = false;
                            break;
                        }
                    }

                    if (!classIndependent) {
                        Matrix nvec = PopulationLattice.pprod(kvec);
                        nvec = PopulationLattice.pprod(nvec, kvec);
                        double sumOfAllProbs = 0.0;
                        while (nvec != null && nvec.value() >= 0) {
                            int hnvec = PopulationLattice.hashpop(nvec, Nc, C, prods);

                            double prob = 0.0;
                            for (int r = 0; r < C; r++) {
                                if (nvec.get(r) > 0) {
                                    int hnvec_c = PopulationLattice.hashpop(oner(nvec, r), Nc, C, prods);
                                    int hkvec_c = PopulationLattice.hashpop(oner(kvec, r), Nc);
                                    double Bcn = getBcnExt(alphas.get(i + ":" + r), D, i, r, nvec, C, (double) ns);
                                    double capacity_inv = 1 / nvec.elementSum();
                                    double x_ir = x[i].get(r, hkvec);
                                    double prob_c = Pc.get(i).get(hnvec_c, hkvec_c);
                                    double classProb = (Bcn * capacity_inv * x_ir * prob_c);
                                    prob += classProb;
                                }
                            }
                            Pc.get(i).set(hnvec, hkvec, prob);
                            sumOfAllProbs += prob;
                            nvec = PopulationLattice.pprod(nvec, kvec);
                        }
                        Pc.get(i).set(0, hkvec, Math.max(1e-12, 1.0 - sumOfAllProbs));
                    } else if (ns > 1) {
                        double meanQueueLength = 0.0;
                        for (int r = 0; r < R; r++) {
                            meanQueueLength += L.get(i).get(r, hkvec);
                        }
                        for (int n = 1; n < ns; n++) {
                            double x1 = (double) CombinatoricsUtils.binomialCoefficient(K_j, n);
                            double frac = meanQueueLength / K_j;
                            double x2 = Math.pow(frac, (double) n);
                            double x3 = Math.pow(1 - frac, (double) (K_j - n));
                            double prob = x1 * x2 * x3;
                            Pc.get(i).set(n, hkvec, prob);
                        }
                        double sum1 = 0.0;
                        double sum2 = 0.0;
                        for (int r = 0; r < R; r++) {
                            sum1 += D.get(i, r) * x[i].get(r, hkvec);
                        }
                        for (int n = 0; n < ns; n++) {
                            sum2 += (ns - n) * Pc.get(i).get(n, hkvec);
                        }
                        double term = (sum1 + sum2) / ns;
                        Pc.get(i).set(0, hkvec, Math.max(1e-12, 1 - term));
                    }
                } else if (ss == SchedStrategy.PS) {
                    if (ns > 1) {
                        double meanQueueLength = 0.0;
                        for (int r = 0; r < R; r++) {
                            meanQueueLength += L.get(i).get(r, hkvec);
                        }
                        for (int n = 1; n < ns; n++) {
                            double x1 = (double) CombinatoricsUtils.binomialCoefficient(K_j, n);
                            double frac = meanQueueLength / K_j;
                            double x2 = Math.pow(frac, (double) n);
                            double x3 = Math.pow(1 - frac, (double) (K_j - n));
                            double prob = x1 * x2 * x3;
                            Pc.get(i).set(n, hkvec, prob);
                        }
                        double sum1 = 0.0;
                        double sum2 = 0.0;
                        for (int r = 0; r < R; r++) {
                            sum1 += D.get(i, r) * x[i].get(r, hkvec);
                        }
                        for (int n = 0; n < ns; n++) {
                            sum2 += (ns - n) * Pc.get(i).get(n, hkvec);
                        }
                        double term = (sum1 + sum2) / ns;
                        Pc.get(i).set(0, hkvec, Math.max(1e-12, 1 - term));
                    }
                }
            }

            kvec = PopulationLattice.pprod(kvec, Nc);
        }

        int hkvecFinal = hkvec;

        for (int c = 0; c < C; c++) {
            double totalResponseTime = 0.0;
            for (int i = 0; i < M; i++) {
                totalResponseTime += w[i].get(c, hkvecFinal);
            }
            XN.set(c, Nc.get(c) / totalResponseTime);
        }

        for (int m = 0; m < M; m++) {
            for (int c = 0; c < C; c++) {
                UN.set(m, c, (D.get(m, c) * XN.get(c)) / S.get(m));
            }
        }

        for (int m = 0; m < M; m++) {
            for (int c = 0; c < C; c++) {
                CN.set(m, c, w[m].get(c, hkvecFinal));
            }
        }

        for (int m = 0; m < M; m++) {
            Matrix Lmat = L.get(m);
            for (int c = 0; c < C; c++) {
                QN.set(m, c, Lmat.get(c, hkvecFinal));
            }
        }

        return new Ret.pfqnAMVASchmidt(XN, QN, UN, CN, PN);
    }

    // Helper functions

    private static boolean allLE(Matrix kvec, Matrix nc) {
        for (int i = 0; i < kvec.length(); i++) {
            if (kvec.get(i) > nc.get(i)) {
                return false;
            }
        }
        return true;
    }

    private static boolean allGE(Matrix kvec, int value) {
        for (int j = 0; j < kvec.length(); j++) {
            if (kvec.get(j) < value) {
                return false;
            }
        }
        return true;
    }

    private static double prod(Matrix n) {
        if (n.isEmpty()) {
            return 1.0;
        }
        double product = 1.0;
        for (int i = 0; i < n.length(); i++) {
            if (Double.isInfinite(n.get(i))) {
                return (double) Integer.MAX_VALUE;
            }
            product *= (n.get(i) + 1);
        }
        return product;
    }

    private static Matrix oner(Matrix kvec, int c) {
        Matrix result = new Matrix(kvec);
        result.set(c, result.get(c) - 1);
        return result;
    }

    private static double getBcn(Matrix D, int i, int c, Matrix nvec, int C, double ns) {
        double Bcn = D.get(i, c);
        if (nvec.elementSum() > 1.0) {
            double eps = 1e-12;
            double sum = 0.0;
            for (int t = 0; t < C; t++) {
                sum += nvec.get(t) * D.get(i, t);
            }
            Bcn += (Math.max(0.0, nvec.elementSum() - ns) / Math.max(ns * (nvec.elementSum() - 1), eps) * (sum - D.get(i, c)));
        }
        return Bcn;
    }

    private static double getBcnExt(Matrix u, Matrix D, int i, int c, Matrix nvec, int C, double ns) {
        double weightedProb = 0.0;
        double totalNonPinnedTime = u.sumRows(i) - u.get(i, C);

        for (int s = 0; s < C; s++) {
            double prob = u.get(i, s) / totalNonPinnedTime;
            weightedProb += (prob / D.get(i, s));
        }

        double meanInterdepartureTime = 1.0 / (ns * weightedProb);
        double Bcn = D.get(i, c);

        if (nvec.elementSum() > 1.0) {
            Bcn += (Math.max(0.0, nvec.elementSum() - ns) * meanInterdepartureTime);
        }

        if (Double.isNaN(Bcn) || Double.isInfinite(Bcn)) {
            Bcn = 0.0;
        }

        return Bcn;
    }
}
