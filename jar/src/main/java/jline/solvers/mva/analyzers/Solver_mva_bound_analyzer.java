/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.mva.analyzers;

import jline.api.pfqn.mva.Pfqn_qzgblow;
import jline.api.pfqn.mva.Pfqn_qzgbup;
import jline.api.pfqn.mva.Pfqn_xzgsblow;
import jline.api.pfqn.mva.Pfqn_xzgsbup;
import jline.api.pfqn.Pfqn_harel_bounds;
import jline.api.pfqn.Pfqn_mwrbb;
import jline.api.pfqn.mva.Pfqn_pbh;
import jline.api.pfqn.mva.Pfqn_cbh;
import jline.api.pfqn.mva.Pfqn_pbk;
import jline.api.pfqn.mva.Pfqn_bjbk;
import jline.api.pfqn.mva.Pfqn_ssd;
import jline.api.pfqn.mva.Pfqn_sib;
import jline.api.pfqn.mva.Pfqn_ldbcmp;
import jline.api.pfqn.mva.Pfqn_scb;
import jline.api.pfqn.mva.Pfqn_mcub;
import jline.api.pfqn.mva.Pfqn_looping;
import jline.api.sn.SnGetDemandsChain;
import jline.api.sn.SnDeaggregateChainResults;
import jline.io.Ret;
import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.solvers.SolverOptions;
import jline.solvers.mva.MVAResult;
import jline.util.Maths;
import jline.util.matrix.Matrix;
import org.apache.commons.math3.util.FastMath;

/**
 * MVA Analyzer class for bounding methods.
 */
public final class Solver_mva_bound_analyzer {

    private Solver_mva_bound_analyzer() {
    }

    public static MVAResult solver_mva_bound_analyzer(NetworkStruct sn, SolverOptions options) {
        MVAResult res = new MVAResult();
        long startTime = System.nanoTime();
        long endTime = startTime;
        int iter = 1;
        String method = options.method;
        if ("auto.upper".equals(method) || "auto.lower".equals(method)) {
            return baAuto(sn, options, method);
        }
        Matrix QN = new Matrix(0, 0);
        Matrix UN = new Matrix(0, 0);
        Matrix RN = new Matrix(0, 0);
        Matrix TN = new Matrix(0, 0);
        Matrix CN = new Matrix(0, 0);
        Matrix WN = new Matrix(0, 0);
        Matrix XN = new Matrix(0, 0);
        double lG = Double.NaN;

        if ("aba.upper".equals(method)) {
            if (sn.nclasses == 1 && sn.nclosedjobs > 0) {
                {
                    int i = 0;
                    while (i < sn.nservers.getNumRows()) {
                        if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF && sn.nservers.get(i) > 1) {
                            throw new RuntimeException("Unsupported method for a model with multi-server stations.");
                        }
                        i++;
                    }
                }
                Matrix V = sn.visits.get(0);
                double Z = 0.0;
                int nondelays = 0;
                {
                    int i = 0;
                    while (i < V.getNumRows()) {
                        if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                            Z += V.get(i) / sn.rates.get(i);
                        } else {
                            nondelays++;
                        }
                        i++;
                    }
                }
                Matrix D = new Matrix(nondelays, 1);
                int idx = 0;
                {
                    int i = 0;
                    while (i < V.getNumRows()) {
                        if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF) {
                            D.set(idx, 0, V.get(i) / sn.rates.get(i));
                            idx++;
                        }
                        i++;
                    }
                }
                double N = (double) sn.nclosedjobs;
                double Dmax = D.elementMax();
                double Dsum = D.elementSum();
                CN = new Matrix(1, 1);
                CN.set(0, 0, Z + N * Dsum);
                XN = new Matrix(1, 1);
                XN.set(0, 0, Maths.min(1 / Dmax, N / (Z + Dsum)));
                TN = new Matrix(V.getNumRows(), 1);
                {
                    int i = 0;
                    while (i < V.getNumRows()) {
                        TN.set(i, 0, V.get(i, 0) * XN.value());
                        i++;
                    }
                }
                RN = new Matrix(sn.rates.getNumRows(), 1);
                {
                    int i = 0;
                    while (i < RN.getNumRows()) {
                        if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                            RN.set(i, 0, 1.0 / sn.rates.get(i));
                        } else {
                            RN.set(i, 0, 1 / sn.rates.get(i) * N);
                        }
                        i++;
                    }
                }
                QN = new Matrix(TN.getNumRows(), 1);
                {
                    int i = 0;
                    while (i < TN.getNumRows()) {
                        QN.set(i, 0, TN.get(i, 0) * RN.get(i, 0));
                        i++;
                    }
                }
                UN = new Matrix(TN.getNumRows(), 1);
                int i = 0;
                while (i < TN.getNumRows()) {
                    if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                        UN.set(i, 0, QN.get(i, 0));
                    } else {
                        UN.set(i, 0, TN.get(i, 0) / sn.rates.get(i));
                    }
                    i++;
                }
                lG = -N * FastMath.log(XN.value());
            }
            endTime = System.nanoTime();
        } else if ("aba.lower".equals(method)) {
            if (sn.nclasses == 1 && sn.nclosedjobs > 0) {
                {
                    int i = 0;
                    while (i < sn.nservers.getNumRows()) {
                        if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF && sn.nservers.get(i) > 1) {
                            throw new RuntimeException("Unsupported method for a model with multi-server stations.");
                        }
                        i++;
                    }
                }
                Matrix V = sn.visits.get(0);
                double Z = 0.0;
                int nondelays = 0;
                {
                    int i = 0;
                    while (i < V.getNumRows()) {
                        if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                            Z += V.get(i) / sn.rates.get(i);
                        } else {
                            nondelays++;
                        }
                        i++;
                    }
                }
                Matrix D = new Matrix(nondelays, 1);
                int idx = 0;
                {
                    int i = 0;
                    while (i < V.getNumRows()) {
                        if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF) {
                            D.set(idx, 0, V.get(i) / sn.rates.get(i));
                            idx++;
                        }
                        i++;
                    }
                }
                double N = (double) sn.nclosedjobs;
                double Dsum = D.elementSum();
                XN = new Matrix(1, 1);
                XN.set(0, 0, N / (Z + N * Dsum));
                CN = new Matrix(1, 1);
                CN.set(0, 0, Z + Dsum);
                TN = new Matrix(V.getNumRows(), 1);
                {
                    int i = 0;
                    while (i < V.getNumRows()) {
                        TN.set(i, 0, V.get(i, 0) * XN.value());
                        i++;
                    }
                }
                RN = new Matrix(sn.rates.getNumRows(), 1);
                {
                    int i = 0;
                    while (i < RN.getNumRows()) {
                        RN.set(i, 0, 1.0 / sn.rates.get(i));
                        i++;
                    }
                }
                QN = new Matrix(TN.getNumRows(), 1);
                {
                    int i = 0;
                    while (i < TN.getNumRows()) {
                        QN.set(i, 0, TN.get(i, 0) * RN.get(i, 0));
                        i++;
                    }
                }
                UN = new Matrix(TN.getNumRows(), 1);
                int i = 0;
                while (i < TN.getNumRows()) {
                    if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                        UN.set(i, 0, QN.get(i, 0));
                    } else {
                        UN.set(i, 0, TN.get(i, 0) / sn.rates.get(i));
                    }
                    i++;
                }
                lG = -N * FastMath.log(XN.value());
            }
            endTime = System.nanoTime();
        } else if ("bjb.upper".equals(method)) {
            if (sn.nclasses == 1 && sn.nclosedjobs > 0) {
                {
                    int i = 0;
                    while (i < sn.nservers.getNumRows()) {
                        if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF && sn.nservers.get(i) > 1) {
                            throw new RuntimeException("Unsupported method for a model with multi-server stations.");
                        }
                        i++;
                    }
                }
                Matrix V = sn.visits.get(0);
                double Z = 0.0;
                int nondelays = 0;
                {
                    int i = 0;
                    while (i < V.getNumRows()) {
                        if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                            Z += V.get(i) / sn.rates.get(i);
                        } else {
                            nondelays++;
                        }
                        i++;
                    }
                }
                Matrix D = new Matrix(nondelays, 1);
                int idx = 0;
                {
                    int i = 0;
                    while (i < V.getNumRows()) {
                        if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF) {
                            D.set(idx, 0, V.get(i) / sn.rates.get(i));
                            idx++;
                        }
                        i++;
                    }
                }
                double N = (double) sn.nclosedjobs;
                double Dmax = D.elementMax();
                double Dsum = D.elementSum();

                double Xaba_upper_1 = Maths.min(1 / Dmax, (N - 1) / (Z + Dsum));
                double Xaba_lower_1 = (N - 1) / (Z + (N - 1) * Dsum);

                CN = new Matrix(1, 1);
                CN.set(0, 0, Z + Dsum + Dmax * (N - 1 - Z * Xaba_lower_1));
                XN = new Matrix(1, 1);
                XN.set(0, 0, Maths.min(1 / Dmax, N / (Z + Dsum + D.meanCol().value() * (N - 1 - Z * Xaba_upper_1))));
                TN = new Matrix(V.getNumRows(), 1);
                {
                    int i = 0;
                    while (i < V.getNumRows()) {
                        TN.set(i, 0, V.get(i, 0) * XN.value());
                        i++;
                    }
                }
                RN = new Matrix(sn.rates.getNumRows(), 1);
                {
                    int i = 0;
                    while (i < RN.getNumRows()) {
                        if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                            RN.set(i, 0, 1.0 / sn.rates.get(i));
                        } else {
                            RN.set(i, 0, 1 / sn.rates.get(i) * N);
                        }
                        i++;
                    }
                }
                QN = new Matrix(TN.getNumRows(), 1);
                {
                    int i = 0;
                    while (i < TN.getNumRows()) {
                        QN.set(i, 0, TN.get(i, 0) * RN.get(i, 0));
                        i++;
                    }
                }
                UN = new Matrix(TN.getNumRows(), 1);
                int i = 0;
                while (i < TN.getNumRows()) {
                    if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                        UN.set(i, 0, QN.get(i, 0));
                    } else {
                        UN.set(i, 0, TN.get(i, 0) / sn.rates.get(i));
                    }
                    i++;
                }
                lG = -N * FastMath.log(XN.value());
            }
            endTime = System.nanoTime();
        } else if ("bjb.lower".equals(method)) {
            if (sn.nclasses == 1 && sn.nclosedjobs > 0) {
                {
                    int i = 0;
                    while (i < sn.nservers.getNumRows()) {
                        if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF && sn.nservers.get(i) > 1) {
                            throw new RuntimeException("Unsupported method for a model with multi-server stations.");
                        }
                        i++;
                    }
                }
                Matrix V = sn.visits.get(0);
                double Z = 0.0;
                int nondelays = 0;
                {
                    int i = 0;
                    while (i < V.getNumRows()) {
                        if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                            Z += V.get(i) / sn.rates.get(i);
                        } else {
                            nondelays++;
                        }
                        i++;
                    }
                }
                Matrix D = new Matrix(nondelays, 1);
                int idx = 0;
                {
                    int i = 0;
                    while (i < V.getNumRows()) {
                        if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF) {
                            D.set(idx, 0, V.get(i) / sn.rates.get(i));
                            idx++;
                        }
                        i++;
                    }
                }
                double N = (double) sn.nclosedjobs;
                double Dmax = D.elementMax();
                double Dsum = D.elementSum();

                double Xaba_upper_1 = Maths.min(1 / Dmax, (N - 1) / (Z + Dsum));
                double Xaba_lower_1 = (N - 1) / (Z + (N - 1) * Dsum);

                CN = new Matrix(1, 1);
                CN.set(0, 0, Z + Dsum + D.meanCol().value() * (N - 1 - Z * Xaba_upper_1));
                XN = new Matrix(1, 1);
                XN.set(0, 0, N / (Z + Dsum + Dmax * (N - 1 - Z * Xaba_lower_1)));
                TN = new Matrix(V.getNumRows(), 1);
                {
                    int i = 0;
                    while (i < V.getNumRows()) {
                        TN.set(i, 0, V.get(i, 0) * XN.value());
                        i++;
                    }
                }
                RN = new Matrix(sn.rates.getNumRows(), 1);
                {
                    int i = 0;
                    while (i < RN.getNumRows()) {
                        RN.set(i, 0, 1.0 / sn.rates.get(i));
                        i++;
                    }
                }
                QN = new Matrix(TN.getNumRows(), 1);
                {
                    int i = 0;
                    while (i < TN.getNumRows()) {
                        QN.set(i, 0, TN.get(i, 0) * RN.get(i, 0));
                        i++;
                    }
                }
                UN = new Matrix(TN.getNumRows(), 1);
                int i = 0;
                while (i < TN.getNumRows()) {
                    if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                        UN.set(i, 0, QN.get(i, 0));
                    } else {
                        UN.set(i, 0, TN.get(i, 0) / sn.rates.get(i));
                    }
                    i++;
                }
                lG = -N * FastMath.log(XN.value());
            }
            endTime = System.nanoTime();
        } else if ("pb.upper".equals(method)) {
            if (sn.nclasses == 1 && sn.nclosedjobs > 0) {
                {
                    int i = 0;
                    while (i < sn.nservers.getNumRows()) {
                        if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF && sn.nservers.get(i) > 1) {
                            throw new RuntimeException("Unsupported method for a model with multi-server stations.");
                        }
                        i++;
                    }
                }
                Matrix V = sn.visits.get(0);
                double Z = 0.0;
                int nondelays = 0;
                {
                    int i = 0;
                    while (i < V.getNumRows()) {
                        if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                            Z += V.get(i) / sn.rates.get(i);
                        } else {
                            nondelays++;
                        }
                        i++;
                    }
                }
                Matrix D = new Matrix(nondelays, 1);
                int idx = 0;
                {
                    int i = 0;
                    while (i < V.getNumRows()) {
                        if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF) {
                            D.set(idx, 0, V.get(i) / sn.rates.get(i));
                            idx++;
                        }
                        i++;
                    }
                }
                double N = (double) sn.nclosedjobs;
                double Dmax = D.elementMax();
                double Dsum = D.elementSum();

                double Xaba_upper_1 = Maths.min(1 / Dmax, (N - 1) / (Z + Dsum));
                double Xaba_lower_1 = (N - 1) / (Z + (N - 1) * Dsum);

                double Dpb2 = D.elementPower(2.0).elementSum() / Dsum;
                double DpbN = D.elementPower(N).elementSum() / D.elementPower(N - 1).elementSum();

                CN = new Matrix(1, 1);
                CN.set(0, 0, Z + Dsum + DpbN * (N - 1 - Z * Xaba_lower_1));
                XN = new Matrix(1, 1);
                XN.set(0, 0, Maths.min(1 / Dmax, N / (Z + Dsum + Dpb2 * (N - 1 - Z * Xaba_upper_1))));
                TN = new Matrix(V.getNumRows(), 1);
                {
                    int i = 0;
                    while (i < V.getNumRows()) {
                        TN.set(i, 0, V.get(i, 0) * XN.value());
                        i++;
                    }
                }
                RN = new Matrix(sn.rates.getNumRows(), 1);
                {
                    int i = 0;
                    while (i < RN.getNumRows()) {
                        if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                            RN.set(i, 0, 1.0 / sn.rates.get(i));
                        } else {
                            RN.set(i, 0, 1 / sn.rates.get(i) * N);
                        }
                        i++;
                    }
                }
                QN = new Matrix(TN.getNumRows(), 1);
                {
                    int i = 0;
                    while (i < TN.getNumRows()) {
                        QN.set(i, 0, TN.get(i, 0) * RN.get(i, 0));
                        i++;
                    }
                }
                UN = new Matrix(TN.getNumRows(), 1);
                int i = 0;
                while (i < TN.getNumRows()) {
                    if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                        UN.set(i, 0, QN.get(i, 0));
                    } else {
                        UN.set(i, 0, TN.get(i, 0) / sn.rates.get(i));
                    }
                    i++;
                }
                lG = -N * FastMath.log(XN.value());
            }
            endTime = System.nanoTime();
        } else if ("pb.lower".equals(method)) {
            if (sn.nclasses == 1 && sn.nclosedjobs > 0) {
                {
                    int i = 0;
                    while (i < sn.nservers.getNumRows()) {
                        if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF && sn.nservers.get(i) > 1) {
                            throw new RuntimeException("Unsupported method for a model with multi-server stations.");
                        }
                        i++;
                    }
                }
                Matrix V = sn.visits.get(0);
                double Z = 0.0;
                int nondelays = 0;
                {
                    int i = 0;
                    while (i < V.getNumRows()) {
                        if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                            Z += V.get(i) / sn.rates.get(i);
                        } else {
                            nondelays++;
                        }
                        i++;
                    }
                }
                Matrix D = new Matrix(nondelays, 1);
                int idx = 0;
                {
                    int i = 0;
                    while (i < V.getNumRows()) {
                        if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF) {
                            D.set(idx, 0, V.get(i) / sn.rates.get(i));
                            idx++;
                        }
                        i++;
                    }
                }
                double N = (double) sn.nclosedjobs;
                double Dmax = D.elementMax();
                double Dsum = D.elementSum();

                double Xaba_upper_1 = Maths.min(1 / Dmax, (N - 1) / (Z + Dsum));
                double Xaba_lower_1 = (N - 1) / (Z + (N - 1) * Dsum);

                double Dpb2 = D.elementPower(2.0).elementSum() / Dsum;
                double DpbN = D.elementPower(N).elementSum() / D.elementPower(N - 1).elementSum();

                CN = new Matrix(1, 1);
                CN.set(0, 0, Z + Dsum + Dpb2 * (N - 1 - Z * Xaba_upper_1));
                XN = new Matrix(1, 1);
                XN.set(0, 0, N / (Z + Dsum + DpbN * (N - 1 - Z * Xaba_lower_1)));
                TN = new Matrix(V.getNumRows(), 1);
                {
                    int i = 0;
                    while (i < V.getNumRows()) {
                        TN.set(i, 0, V.get(i, 0) * XN.value());
                        i++;
                    }
                }
                RN = new Matrix(sn.rates.getNumRows(), 1);
                {
                    int i = 0;
                    while (i < RN.getNumRows()) {
                        RN.set(i, 0, 1.0 / sn.rates.get(i));
                        i++;
                    }
                }
                QN = new Matrix(TN.getNumRows(), 1);
                {
                    int i = 0;
                    while (i < TN.getNumRows()) {
                        QN.set(i, 0, TN.get(i, 0) * RN.get(i, 0));
                        i++;
                    }
                }
                UN = new Matrix(TN.getNumRows(), 1);
                int i = 0;
                while (i < TN.getNumRows()) {
                    if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                        UN.set(i, 0, QN.get(i, 0));
                    } else {
                        UN.set(i, 0, TN.get(i, 0) / sn.rates.get(i));
                    }
                    i++;
                }
                lG = -N * FastMath.log(XN.value());
            }
            endTime = System.nanoTime();
        } else if ("sb.upper".equals(method)) {
            if (sn.nclasses == 1 && sn.nclosedjobs > 0) {
                {
                    int i = 0;
                    while (i < sn.nservers.getNumRows()) {
                        if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF && sn.nservers.get(i) > 1) {
                            throw new RuntimeException("Unsupported method for a model with multi-server stations.");
                        }
                        if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                            throw new RuntimeException("Unsupported method for a model with infinite-server stations.");
                        }
                        i++;
                    }
                }
                Matrix V = sn.visits.get(0);
                double Z = 0.0;
                int nondelays = 0;
                {
                    int i = 0;
                    while (i < V.getNumRows()) {
                        if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                            Z += V.get(i) / sn.rates.get(i);
                        } else {
                            nondelays++;
                        }
                        i++;
                    }
                }
                Matrix D = new Matrix(nondelays, 1);
                int idx = 0;
                {
                    int i = 0;
                    while (i < V.getNumRows()) {
                        if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF) {
                            D.set(idx, 0, V.get(i) / sn.rates.get(i));
                            idx++;
                        }
                        i++;
                    }
                }
                double N = (double) sn.nclosedjobs;
                double Dmax = D.elementMax();

                double A3 = D.elementPower(3.0).elementSum();
                double A2 = D.elementPower(2.0).elementSum();
                double A1 = D.elementPower(1.0).elementSum();

                // Harel UB(n) is defined for n <= N only; the level-3 coefficient is
                // not a bound at N < 3, so fall back to UB(2) (Dallery), exact there
                // since UB(N) = TH(N).
                double cub3 = (N >= 3) ? (A1 * A2 + A3) / (A1 * A1 + A2) : A2 / A1;
                CN = new Matrix(1, 1);
                CN.set(0, 0, Z + A1 + (N - 1) * cub3);
                XN = new Matrix(1, 1);
                XN.set(0, 0, Maths.min(1 / Dmax, N / CN.value()));
                TN = new Matrix(V.getNumRows(), 1);
                {
                    int i = 0;
                    while (i < V.getNumRows()) {
                        TN.set(i, 0, V.get(i, 0) * XN.value());
                        i++;
                    }
                }
                // RN is undefined in the literature for this bound, so it carries
                // the ABA PESSIMISTIC residence: an upper method must publish an
                // upper-consistent R, else QN = TN*RN lands below exact
                // (E[n_i] <= N*U_i makes N/mu_i a valid residence bound).
                RN = new Matrix(sn.rates.getNumRows(), 1);
                {
                    int i = 0;
                    while (i < RN.getNumRows()) {
                        RN.set(i, 0, N / sn.rates.get(i));
                        i++;
                    }
                }
                QN = new Matrix(TN.getNumRows(), 1);
                {
                    int i = 0;
                    while (i < TN.getNumRows()) {
                        QN.set(i, 0, TN.get(i, 0) * RN.get(i, 0));
                        i++;
                    }
                }
                UN = new Matrix(TN.getNumRows(), 1);
                int i = 0;
                while (i < TN.getNumRows()) {
                    if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                        UN.set(i, 0, QN.get(i, 0));
                    } else {
                        UN.set(i, 0, TN.get(i, 0) / sn.rates.get(i));
                    }
                    i++;
                }
                lG = -N * FastMath.log(XN.value());
            }
            endTime = System.nanoTime();
        } else if ("sb.lower".equals(method)) {
            if (sn.nclasses == 1 && sn.nclosedjobs > 0) {
                {
                    int i = 0;
                    while (i < sn.nservers.getNumRows()) {
                        if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF && sn.nservers.get(i) > 1) {
                            throw new RuntimeException("Unsupported method for a model with multi-server stations.");
                        }
                        if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                            throw new RuntimeException("Unsupported method for a model with infinite-server stations.");
                        }
                        i++;
                    }
                }
                Matrix V = sn.visits.get(0);
                double Z = 0.0;
                int nondelays = 0;
                {
                    int i = 0;
                    while (i < V.getNumRows()) {
                        if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                            Z += V.get(i) / sn.rates.get(i);
                        } else {
                            nondelays++;
                        }
                        i++;
                    }
                }
                Matrix D = new Matrix(nondelays, 1);
                int idx = 0;
                {
                    int i = 0;
                    while (i < V.getNumRows()) {
                        if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF) {
                            D.set(idx, 0, V.get(i) / sn.rates.get(i));
                            idx++;
                        }
                        i++;
                    }
                }
                double N = (double) sn.nclosedjobs;
                D.elementMax();

                double AN = D.elementPower(N).elementSum();
                double A1 = D.elementPower(1.0).elementSum();

                // (N-1)*(AN/A1)^(1/(N-1)) -> 0 as N -> 1, leaving the exact
                // single-job cycle time; evaluated directly it divides by zero
                double cterm = (N == 1) ? 0.0 : (N - 1) * FastMath.pow(AN / A1, 1 / (N - 1));
                CN = new Matrix(1, 1);
                CN.set(0, 0, Z + A1 + cterm);
                XN = new Matrix(1, 1);
                XN.set(0, 0, N / CN.value());
                TN = new Matrix(V.getNumRows(), 1);
                {
                    int i = 0;
                    while (i < V.getNumRows()) {
                        TN.set(i, 0, V.get(i, 0) * XN.value());
                        i++;
                    }
                }
                RN = new Matrix(sn.rates.getNumRows(), 1);
                {
                    int i = 0;
                    while (i < RN.getNumRows()) {
                        RN.set(i, 0, 1.0 / sn.rates.get(i));
                        i++;
                    }
                }
                QN = new Matrix(TN.getNumRows(), 1);
                {
                    int i = 0;
                    while (i < TN.getNumRows()) {
                        QN.set(i, 0, TN.get(i, 0) * RN.get(i, 0));
                        i++;
                    }
                }
                UN = new Matrix(TN.getNumRows(), 1);
                int i = 0;
                while (i < TN.getNumRows()) {
                    if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                        UN.set(i, 0, QN.get(i, 0));
                    } else {
                        UN.set(i, 0, TN.get(i, 0) / sn.rates.get(i));
                    }
                    i++;
                }
                lG = -N * FastMath.log(XN.value());
            }
            endTime = System.nanoTime();
        } else if ("gb.upper".equals(method)) {
            if (sn.nclasses == 1 && sn.nclosedjobs > 0) {
                {
                    int i = 0;
                    while (i < sn.nservers.getNumRows()) {
                        if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF && sn.nservers.get(i) > 1) {
                            throw new RuntimeException("Unsupported method for a model with multi-server stations.");
                        }
                        i++;
                    }
                }
                Matrix V = sn.visits.get(0);
                double Z = 0.0;
                int nondelays = 0;
                {
                    int i = 0;
                    while (i < V.getNumRows()) {
                        if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                            Z += V.get(i) / sn.rates.get(i);
                        } else {
                            nondelays++;
                        }
                        i++;
                    }
                }
                Matrix D = new Matrix(nondelays, 1);
                int idx = 0;
                {
                    int i = 0;
                    while (i < V.getNumRows()) {
                        if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF) {
                            D.set(idx, 0, V.get(i) / sn.rates.get(i));
                            idx++;
                        }
                        i++;
                    }
                }
                double N = (double) sn.nclosedjobs;
                double Dmax = D.elementMax();
                XN = new Matrix(1, 1);
                XN.set(0, 0, Maths.min(1 / Dmax, Pfqn_xzgsbup.pfqn_xzgsbup(D, N, Z)));
                double ret = Pfqn_xzgsblow.pfqn_xzgsblow(D, N, Z);
                CN = new Matrix(1, 1);
                CN.set(0, 0, N / ret);
                TN = new Matrix(V.getNumRows(), 1);
                {
                    int i = 0;
                    while (i < V.getNumRows()) {
                        TN.set(i, 0, V.get(i, 0) * XN.value());
                        i++;
                    }
                }
                double XNlow = ret;
                int k = 0;
                RN = new Matrix(sn.sched.size(), 1);
                QN = new Matrix(sn.sched.size(), 1);
                {
                    int i = 0;
                    while (i < sn.sched.size()) {
                        if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                            RN.set(i, 0, 1.0 / sn.rates.get(i));
                            // TN, not XN: the delay is visited V(i) times per
                            // cycle, and dropping that factor lets QN exceed
                            // the population
                            QN.set(i, 0, TN.get(i, 0) * RN.get(i, 0));
                        } else {
                            QN.set(i, 0, Pfqn_qzgbup.pfqn_qzgbup(D, N, Z, k));
                            RN.set(i, 0, QN.get(i, 0) / XNlow / V.get(i));
                            k++;
                        }
                        i++;
                    }
                }
                UN = new Matrix(TN.getNumRows(), 1);
                {
                    int i = 0;
                    while (i < UN.getNumRows()) {
                        UN.set(i, 0, TN.get(i, 0) / sn.rates.get(i));
                        i++;
                    }
                }
                int i = 0;
                while (i < sn.sched.size()) {
                    if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                        UN.set(i, 0, QN.get(i, 0));
                    }
                    i++;
                }
                lG = -N * FastMath.log(XN.value());
            }
            endTime = System.nanoTime();
        } else if ("gb.lower".equals(method)) {
            if (sn.nclasses == 1 && sn.nclosedjobs > 0) {
                {
                    int i = 0;
                    while (i < sn.nservers.getNumRows()) {
                        if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF && sn.nservers.get(i) > 1) {
                            throw new RuntimeException("Unsupported method for a model with multi-server stations.");
                        }
                        i++;
                    }
                }
                Matrix V = sn.visits.get(0);
                double Z = 0.0;
                int nondelays = 0;
                {
                    int i = 0;
                    while (i < V.getNumRows()) {
                        if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                            Z += V.get(i) / sn.rates.get(i);
                        } else {
                            nondelays++;
                        }
                        i++;
                    }
                }
                Matrix D = new Matrix(nondelays, 1);
                int idx = 0;
                {
                    int i = 0;
                    while (i < V.getNumRows()) {
                        if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF) {
                            D.set(idx, 0, V.get(i) / sn.rates.get(i));
                            idx++;
                        }
                        i++;
                    }
                }
                double N = (double) sn.nclosedjobs;
                XN = new Matrix(1, 1);
                XN.set(0, 0, Pfqn_xzgsblow.pfqn_xzgsblow(D, N, Z));
                double ret = Pfqn_xzgsbup.pfqn_xzgsbup(D, N, Z);
                CN = new Matrix(1, 1);
                CN.set(0, 0, N / ret);
                TN = new Matrix(V.getNumRows(), 1);
                {
                    int i = 0;
                    while (i < V.getNumRows()) {
                        TN.set(i, 0, V.get(i, 0) * XN.value());
                        i++;
                    }
                }
                double XNup = ret;
                int k = 0;
                RN = new Matrix(sn.sched.size(), 1);
                QN = new Matrix(sn.sched.size(), 1);
                {
                    int i = 0;
                    while (i < sn.sched.size()) {
                        if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                            RN.set(i, 0, 1.0 / sn.rates.get(i));
                            // TN, not XN: the delay is visited V(i) times per
                            // cycle, and dropping that factor lets QN exceed
                            // the population
                            QN.set(i, 0, TN.get(i, 0) * RN.get(i, 0));
                        } else {
                            QN.set(i, 0, Pfqn_qzgblow.pfqn_qzgblow(D, N, Z, k));
                            RN.set(i, 0, QN.get(i, 0) / XNup / V.get(i));
                            k++;
                        }
                        i++;
                    }
                }
                UN = new Matrix(TN.getNumRows(), 1);
                {
                    int i = 0;
                    while (i < UN.getNumRows()) {
                        UN.set(i, 0, TN.get(i, 0) / sn.rates.get(i));
                        i++;
                    }
                }
                int i = 0;
                while (i < sn.sched.size()) {
                    if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                        UN.set(i, 0, QN.get(i, 0));
                    }
                    i++;
                }
                lG = -N * FastMath.log(XN.value());
            }
            endTime = System.nanoTime();
        } else if ("harel.lower".equals(method) || "harel.upper".equals(method)) {
            if (sn.nclasses == 1 && sn.nclosedjobs > 0) {
                {
                    int i = 0;
                    while (i < sn.nservers.getNumRows()) {
                        if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF && sn.nservers.get(i) > 1) {
                            throw new RuntimeException("Unsupported method for a model with multi-server stations.");
                        }
                        if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                            throw new RuntimeException("Harel bounds do not support models with think times (infinite-server stations).");
                        }
                        i++;
                    }
                }
                Matrix V = sn.visits.get(0);
                int nondelays = 0;
                {
                    int i = 0;
                    while (i < V.getNumRows()) {
                        if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF) {
                            nondelays++;
                        }
                        i++;
                    }
                }
                Matrix rho = new Matrix(nondelays, 1);
                int idx = 0;
                {
                    int i = 0;
                    while (i < V.getNumRows()) {
                        if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF) {
                            rho.set(idx, 0, V.get(i) / sn.rates.get(i));
                            idx++;
                        }
                        i++;
                    }
                }
                int N = sn.nclosedjobs;
                double Dmax = rho.elementMax();

                int maxN = Math.min(N, 7);
                Ret.pfqnHarelBounds bounds = Pfqn_harel_bounds.pfqn_harel_bounds(rho, N, 0.0, maxN);

                double throughputBound;
                if ("harel.lower".equals(method)) {
                    throughputBound = bounds.LB;
                } else {
                    if (maxN >= 2) {
                        throughputBound = bounds.UB[maxN];
                    } else {
                        throughputBound = bounds.TH[1];
                    }
                }

                XN = new Matrix(1, 1);
                if ("harel.upper".equals(method)) {
                    XN.set(0, 0, Maths.min(1 / Dmax, throughputBound));
                } else {
                    XN.set(0, 0, throughputBound);
                }

                // Same fill as every other bound family here: the optimistic
                // side charges the full-contention residence time N/mu, the
                // pessimistic side the no-contention 1/mu, and the cycle time
                // is the ABA one at that side. Harel forbids think times, so Z
                // is zero in both expressions.
                boolean upperSide = "harel.upper".equals(method);
                CN = new Matrix(1, 1);
                CN.set(0, 0, upperSide ? N * rho.elementSum() : rho.elementSum());

                TN = new Matrix(V.getNumRows(), 1);
                {
                    int i = 0;
                    while (i < V.getNumRows()) {
                        TN.set(i, 0, V.get(i, 0) * XN.value());
                        i++;
                    }
                }
                RN = new Matrix(sn.rates.getNumRows(), 1);
                {
                    int i = 0;
                    while (i < RN.getNumRows()) {
                        RN.set(i, 0, upperSide ? N / sn.rates.get(i) : 1.0 / sn.rates.get(i));
                        i++;
                    }
                }
                QN = new Matrix(TN.getNumRows(), 1);
                {
                    int i = 0;
                    while (i < TN.getNumRows()) {
                        QN.set(i, 0, TN.get(i, 0) * RN.get(i, 0));
                        i++;
                    }
                }
                UN = new Matrix(TN.getNumRows(), 1);
                int i = 0;
                while (i < TN.getNumRows()) {
                    if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                        UN.set(i, 0, QN.get(i, 0));
                    } else {
                        UN.set(i, 0, TN.get(i, 0) / sn.rates.get(i));
                    }
                    i++;
                }
                lG = -N * FastMath.log(XN.value());
            }
            endTime = System.nanoTime();
        }

        if ("mwba.upper".equals(method) || "mwba.lower".equals(method)) {
            boolean isOpen = false;
            for (int r = 0; r < sn.njobs.getNumCols(); r++) {
                if (Double.isInfinite(sn.njobs.get(0, r))) {
                    isOpen = true;
                }
            }
            if (sn.nclosedjobs > 0 && !isOpen) {
                for (int i = 0; i < sn.nservers.getNumRows(); i++) {
                    if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF && sn.nservers.get(i) > 1) {
                        throw new RuntimeException("Unsupported method for a model with multi-server stations.");
                    }
                }
                Ret.snGetDemands chainReturn = SnGetDemandsChain.snGetDemandsChain(sn);
                Matrix Lchain = chainReturn.Dchain;
                Matrix STchain = chainReturn.STchain;
                Matrix Vchain = chainReturn.Vchain;
                Matrix alpha = chainReturn.alpha;
                Matrix Nchain = chainReturn.Nchain;
                int M = sn.nstations;
                int Cc = sn.nchains;

                // FIFO queueing stations vs pure-delay (INF) stations
                boolean[] isdelay = new boolean[M];
                int Kq = 0;
                for (int i = 0; i < M; i++) {
                    isdelay[i] = (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF);
                    if (!isdelay[i]) {
                        Kq++;
                    }
                }
                int[] qrow = new int[Kq];
                Matrix Vq = new Matrix(Kq, Cc);
                Matrix Sq = new Matrix(Kq, Cc);
                Matrix Zc = new Matrix(1, Cc);
                Matrix schedq = new Matrix(Kq, 1);
                int kk = 0;
                for (int i = 0; i < M; i++) {
                    if (isdelay[i]) {
                        for (int c = 0; c < Cc; c++) {
                            Zc.set(0, c, Zc.get(0, c) + Lchain.get(i, c));
                        }
                    } else {
                        qrow[kk] = i;
                        for (int c = 0; c < Cc; c++) {
                            Vq.set(kk, c, Vchain.get(i, c));
                            Sq.set(kk, c, STchain.get(i, c));
                        }
                        schedq.set(kk, 0, mwrbbDiscCode(sn.sched.get(sn.stations.get(i))));
                        kk++;
                    }
                }
                // per-chain priority (lower value = higher priority): use the
                // highest-priority (minimum) classprio among the chain's classes
                Matrix prioc = new Matrix(1, Cc);
                for (int c = 0; c < Cc; c++) {
                    Matrix inchain_c = sn.inchain.get(c);
                    double p = Double.POSITIVE_INFINITY;
                    for (int idx = 0; idx < inchain_c.getNumCols(); idx++) {
                        int cls = (int) inchain_c.get(idx);
                        p = Math.min(p, sn.classprio.get(0, cls));
                    }
                    prioc.set(0, c, Double.isInfinite(p) ? 0.0 : p);
                }

                Pfqn_mwrbb.Result rbb = Pfqn_mwrbb.pfqn_mwrbb(Vq, Sq, Nchain, Zc, schedq, prioc);
                boolean upper = "mwba.upper".equals(method);
                Matrix Xchain = new Matrix(1, Cc);
                for (int c = 0; c < Cc; c++) {
                    Xchain.set(0, c, upper ? rbb.Xup.get(0, c) : rbb.Xlo.get(0, c));
                }

                Matrix Tchain = new Matrix(M, Cc);
                Matrix Uchain = new Matrix(M, Cc);
                Matrix Qchain = new Matrix(M, Cc);
                int[] rowToQ = new int[M];
                for (int i = 0; i < M; i++) {
                    rowToQ[i] = -1;
                }
                for (int k = 0; k < Kq; k++) {
                    rowToQ[qrow[k]] = k;
                }
                for (int c = 0; c < Cc; c++) {
                    for (int i = 0; i < M; i++) {
                        double xc = Xchain.get(0, c);
                        Tchain.set(i, c, xc * Vchain.get(i, c));
                        Uchain.set(i, c, xc * Lchain.get(i, c));   // utilization law
                    }
                }
                // Neither the no-contention residence Sq nor the Theorem-1
                // residence Wlo yields a queue length on the declared side.
                baChainQfill(Qchain, Uchain, Xchain, Lchain, Nchain, isdelay, upper);

                Ret.snDeaggregateChainResults dre = SnDeaggregateChainResults.snDeaggregateChainResults(
                        sn, Lchain, null, STchain, Vchain, alpha, Qchain, Uchain, null, Tchain, null, Xchain);
                QN = dre.Q;
                UN = dre.U;
                RN = dre.R;
                TN = dre.T;
                CN = dre.C;
                XN = dre.X;
                lG = Double.NaN;
            }
            endTime = System.nanoTime();
        }

        // Single-class hierarchical bound families (SolverBA). Each returns a
        // scalar chain-throughput bound; ba_fill builds per-station Q/U/R/T/C.
        if (isScHierMethod(method)) {
            Matrix V = sn.visits.get(0);
            double N = (double) sn.nclosedjobs;
            if (sn.nclasses != 1 || sn.nclosedjobs <= 0) {
                throw new RuntimeException("Method '" + method + "' supports single-class closed networks only.");
            }
            double Zt = baThink(sn, V);
            Matrix D = baQueueDemands(sn, V);
            int level = options.level;
            boolean up = method.endsWith(".upper");
            double X;
            if (method.startsWith("ssd")) {
                Matrix cvec = baQueueServers(sn, D.getNumRows());
                double[] xb = Pfqn_ssd.pfqn_ssd(D, N, Zt, cvec);
                X = up ? xb[1] : xb[0];
            } else if (method.startsWith("ldbcmp")) {
                double[] xb = Pfqn_ldbcmp.pfqn_ldbcmp(D, N, Zt, new Matrix(D.getNumRows(), 1));
                if (Double.isNaN(xb[0])) {
                    throw new RuntimeException(String.format("Method '%s' requires the asymptotic regime N >= Qhat (Qhat=%.4f > N=%d).", method, xb[2], (int) N));
                }
                X = xb[0]; up = false;
            } else if (method.startsWith("scb")) {
                // Single-class bounds of Dowdy et al. (1992). THE BRACKETED OBJECT IS
                // NOT THIS MODEL: scb brackets the multiclass system that this
                // single-class model aggregates, so scb.lower is the EXACT single-class
                // throughput and scb.upper adds the demand-free Expression-(3) gap.
                // That is why scb is absent from BA_AUTO_* -- mixing it with families
                // that bracket this model's own solution would compare two different
                // quantities.
                baRejectMultiserver(sn, method);
                if (Zt > 0) {
                    throw new RuntimeException("Method '" + method + "' supports Z=0 (no delay station) only; Theorem 3 rests on the delay-free balanced-network throughput.");
                }
                double[] xb = Pfqn_scb.pfqn_scb(D, (int) N);
                X = up ? xb[1] : xb[0];
            } else {
                baRejectMultiserver(sn, method);
                double[] xb;
                if (method.startsWith("pbh")) {
                    xb = Pfqn_pbh.pfqn_pbh(D, N, Zt, level);
                } else if (method.startsWith("cbh")) {
                    xb = Pfqn_cbh.pfqn_cbh(D, N, Zt, level);
                } else if (method.startsWith("pbk")) {
                    xb = Pfqn_pbk.pfqn_pbk(D, N, Zt, level);
                } else if (method.startsWith("bjbk")) {
                    xb = Pfqn_bjbk.pfqn_bjbk(D, N, Zt, level);
                } else { // sib
                    if (Zt > 0) {
                        throw new RuntimeException("Method '" + method + "' supports Z=0 (no delay station) only; delay needs the SIB Section-3.2 extension.");
                    }
                    xb = Pfqn_sib.pfqn_sib(D, N, 0.0, level);
                }
                X = up ? xb[1] : xb[0];
            }
            Matrix[] fill = baFill(sn, V, N, X, Zt, D, up);
            QN = fill[0]; UN = fill[1]; RN = fill[2]; TN = fill[3]; CN = fill[4];
            XN = new Matrix(1, 1);
            XN.set(0, 0, X);
            lG = -N * FastMath.log(X);
            endTime = System.nanoTime();
        }

        // Multiclass composite/BJB bounds (cub upper-only; mbjb lower).
        if ("cub.upper".equals(method) || "mbjb.lower".equals(method)) {
            boolean isOpen = false;
            for (int r = 0; r < sn.njobs.getNumCols(); r++) {
                if (Double.isInfinite(sn.njobs.get(0, r))) isOpen = true;
            }
            if (sn.nclosedjobs <= 0 || isOpen) {
                throw new RuntimeException("Method '" + method + "' supports fully closed networks only.");
            }
            for (int i = 0; i < sn.nservers.getNumRows(); i++) {
                if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF && sn.nservers.get(i) > 1) {
                    throw new RuntimeException("Method '" + method + "' does not support multi-server stations (use 'ssd').");
                }
            }
            Ret.snGetDemands cr = SnGetDemandsChain.snGetDemandsChain(sn);
            Matrix Lchain = cr.Dchain;
            Matrix STchain = cr.STchain;
            Matrix Vchain = cr.Vchain;
            Matrix alpha = cr.alpha;
            Matrix Nchain = cr.Nchain;
            int M = sn.nstations;
            int Cc = sn.nchains;
            boolean[] isdelay = new boolean[M];
            int Kq = 0;
            for (int i = 0; i < M; i++) {
                isdelay[i] = (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF);
                if (!isdelay[i]) Kq++;
            }
            Matrix Lq = new Matrix(Kq, Cc);
            Matrix Zc = new Matrix(1, Cc);
            int kk = 0;
            for (int i = 0; i < M; i++) {
                if (isdelay[i]) {
                    for (int c = 0; c < Cc; c++) Zc.set(0, c, Zc.get(0, c) + Lchain.get(i, c));
                } else {
                    for (int c = 0; c < Cc; c++) Lq.set(kk, c, Lchain.get(i, c));
                    kk++;
                }
            }
            Matrix[] mc = Pfqn_mcub.pfqn_mcub(Lq, Nchain, Zc);
            boolean up = "cub.upper".equals(method);
            Matrix Xchain = up ? mc[0] : mc[1];
            Matrix Tchain = new Matrix(M, Cc);
            Matrix Uchain = new Matrix(M, Cc);
            Matrix Qchain = new Matrix(M, Cc);
            for (int c = 0; c < Cc; c++) {
                for (int i = 0; i < M; i++) {
                    double xc = Xchain.get(0, c);
                    Tchain.set(i, c, xc * Vchain.get(i, c));
                    Uchain.set(i, c, xc * Lchain.get(i, c));
                }
            }
            baChainQfill(Qchain, Uchain, Xchain, Lchain, Nchain, isdelay, up);
            Ret.snDeaggregateChainResults dre = SnDeaggregateChainResults.snDeaggregateChainResults(
                    sn, Lchain, null, STchain, Vchain, alpha, Qchain, Uchain, null, Tchain, null, Xchain);
            QN = dre.Q; UN = dre.U; RN = dre.R; TN = dre.T; CN = dre.C; XN = dre.X;
            lG = Double.NaN;
            endTime = System.nanoTime();
        }

        // Eager Looping: the multiclass bracket that initializes the
        // multiple-class PBH. Pessimistic side from the heap-inflated response
        // time, optimistic side from the response-time lower bound.
        if ("looping.upper".equals(method) || "looping.lower".equals(method)) {
            boolean isOpen = false;
            for (int r = 0; r < sn.njobs.getNumCols(); r++) {
                if (Double.isInfinite(sn.njobs.get(0, r))) isOpen = true;
            }
            if (sn.nclosedjobs <= 0 || isOpen) {
                throw new RuntimeException("Method '" + method + "' supports fully closed networks only.");
            }
            for (int i = 0; i < sn.nservers.getNumRows(); i++) {
                if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF && sn.nservers.get(i) > 1) {
                    throw new RuntimeException("Method '" + method + "' does not support multi-server stations (use 'ssd').");
                }
            }
            Ret.snGetDemands cr = SnGetDemandsChain.snGetDemandsChain(sn);
            Matrix Lchain = cr.Dchain;
            Matrix STchain = cr.STchain;
            Matrix Vchain = cr.Vchain;
            Matrix alpha = cr.alpha;
            Matrix Nchain = cr.Nchain;
            int M = sn.nstations;
            int Cc = sn.nchains;
            boolean[] isdelay = new boolean[M];
            int Kq = 0;
            for (int i = 0; i < M; i++) {
                isdelay[i] = (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF);
                if (!isdelay[i]) Kq++;
            }
            Matrix Lq = new Matrix(Kq, Cc);
            Matrix Zc = new Matrix(1, Cc);
            int kk = 0;
            for (int i = 0; i < M; i++) {
                if (isdelay[i]) {
                    for (int c = 0; c < Cc; c++) Zc.set(0, c, Zc.get(0, c) + Lchain.get(i, c));
                } else {
                    for (int c = 0; c < Cc; c++) Lq.set(kk, c, Lchain.get(i, c));
                    kk++;
                }
            }
            Pfqn_looping.Result lp = Pfqn_looping.pfqn_looping(Lq, Nchain, Zc);
            boolean up = "looping.upper".equals(method);
            Matrix Xchain = up ? lp.Xup : lp.Xlo;
            Matrix Tchain = new Matrix(M, Cc);
            Matrix Uchain = new Matrix(M, Cc);
            Matrix Qchain = new Matrix(M, Cc);
            for (int c = 0; c < Cc; c++) {
                for (int i = 0; i < M; i++) {
                    double xc = Xchain.get(0, c);
                    Tchain.set(i, c, xc * Vchain.get(i, c));
                    Uchain.set(i, c, xc * Lchain.get(i, c));
                }
            }
            baChainQfill(Qchain, Uchain, Xchain, Lchain, Nchain, isdelay, up);
            Ret.snDeaggregateChainResults dre = SnDeaggregateChainResults.snDeaggregateChainResults(
                    sn, Lchain, null, STchain, Vchain, alpha, Qchain, Uchain, null, Tchain, null, Xchain);
            QN = dre.Q; UN = dre.U; RN = dre.R; TN = dre.T; CN = dre.C; XN = dre.X;
            lG = Double.NaN;
            endTime = System.nanoTime();
        }

        res.QN = QN;
        res.UN = UN;
        res.RN = RN;
        res.TN = TN;
        res.CN = CN;
        res.XN = XN;
        res.AN = new Matrix(0, 0);
        res.WN = WN;
        res.logNormConstAggr = lG;
        res.runtime = (endTime - startTime) / 1000000000.0;
        res.iter = iter;
        res.method = method;
        return res;
    }

    // Per-chain queue lengths from a chain-throughput bound, on the declared
    // side. Lower: Q_ic >= U_ic, since the station holds a class-c job whenever
    // it serves one. Upper: E[n_ic] <= N_c*P(station busy) = N_c*min(1,sum_c
    // U_ic), and a delay station queues nothing, so there Q_ic = X_c*L_ic.
    private static void baChainQfill(Matrix Qchain, Matrix Uchain, Matrix Xchain,
                                     Matrix Lchain, Matrix Nchain, boolean[] isdelay,
                                     boolean isUpper) {
        int M = Qchain.getNumRows();
        int C = Qchain.getNumCols();
        for (int i = 0; i < M; i++) {
            double utot = 0.0;
            for (int c = 0; c < C; c++) utot += Uchain.get(i, c);
            if (utot > 1.0) utot = 1.0;
            for (int c = 0; c < C; c++) {
                if (isdelay[i]) {
                    Qchain.set(i, c, Xchain.get(0, c) * Lchain.get(i, c));
                } else if (isUpper) {
                    Qchain.set(i, c, Nchain.get(0, c) * utot);
                } else {
                    Qchain.set(i, c, Uchain.get(i, c));
                }
            }
        }
    }

    // Map a SchedStrategy to a Majumdar-Woodside discipline code:
    // 0=FIFO, 1=PS, 2=non-preemptive priority, 3=preemptive priority,
    // 4=ABA full-contention (discipline-independent).
    private static int mwrbbDiscCode(SchedStrategy s) {
        switch (s) {
            case FCFS:
                return Pfqn_mwrbb.FIFO;   // Theorem 1 / Lemma 1
            case PS:
            case DPS:
            case GPS:
            case PSPRIO:
            case DPSPRIO:
            case GPSPRIO:
                return Pfqn_mwrbb.PS;
            case HOL:            // HOL == FCFSPRIO: non-preemptive priority
            case FCFSPRIO:
                return Pfqn_mwrbb.NPPRIO;
            case FCFSPRPRIO:
            case LCFSPRPRIO:
                return Pfqn_mwrbb.PPPRIO;
            default:             // non-FCFS work-conserving: ABA full-contention
                return Pfqn_mwrbb.ABA;
        }
    }

    // Candidate lists for the AUTO composite. Noniterative families only: the
    // level-parameterized hierarchies (pbh/cbh/pbk/bjbk/sib) and the LP
    // reductions are excluded because their cost is not O(K) and their accuracy
    // is a user choice, not a fixed property.
    private static final String[] BA_AUTO_UPPER = {"aba.upper", "bjb.upper", "pb.upper",
            "gb.upper", "sb.upper", "mwba.upper", "ssd.upper", "cub.upper"};
    private static final String[] BA_AUTO_LOWER = {"aba.lower", "bjb.lower", "pb.lower",
            "gb.lower", "sb.lower", "mwba.lower", "ssd.lower", "mbjb.lower", "ldbcmp.lower"};

    /**
     * AUTO composite: evaluate every noniterative bound and keep the tightest
     * side. Feasibility is probed by execution -- a candidate that rejects the
     * model (multiserver, delay station, regime gate) throws and is skipped --
     * so the list stays correct as families are added.
     *
     * @param sn      network structure (single-class closed)
     * @param options caller options; only the method is overridden per candidate
     * @param method  "auto.upper" or "auto.lower"
     * @return per-station metrics built from the tightest feasible bound
     */
    private static MVAResult baAuto(NetworkStruct sn, SolverOptions options, String method) {
        long tstart = System.nanoTime();
        if (sn.nclasses != 1 || sn.nclosedjobs <= 0) {
            throw new RuntimeException("Method '" + method + "' supports single-class closed networks only.");
        }
        Matrix V = sn.visits.get(0);
        double N = (double) sn.nclosedjobs;
        double Zt = baThink(sn, V);
        Matrix D = baQueueDemands(sn, V);
        boolean up = method.endsWith(".upper");
        String[] cand = up ? BA_AUTO_UPPER : BA_AUTO_LOWER;
        double Xbest = Double.NaN;
        for (int c = 0; c < cand.length; c++) {
            SolverOptions oc = options.copy();
            oc.method = cand[c];
            double Xc;
            try {
                MVAResult r = solver_mva_bound_analyzer(sn, oc);
                if (r.XN == null || r.XN.isEmpty()) continue;
                Xc = r.XN.get(0, 0);
            } catch (RuntimeException e) {
                continue;
            }
            if (!Double.isFinite(Xc) || Xc <= 0) continue;
            if (Double.isNaN(Xbest) || (up && Xc < Xbest) || (!up && Xc > Xbest)) {
                Xbest = Xc;
            }
        }
        if (Double.isNaN(Xbest)) {
            throw new RuntimeException("Method '" + method + "' found no feasible bound for this model.");
        }
        Matrix[] fill = baFill(sn, V, N, Xbest, Zt, D, up);
        MVAResult res = new MVAResult();
        res.QN = fill[0]; res.UN = fill[1]; res.RN = fill[2]; res.TN = fill[3]; res.CN = fill[4];
        res.XN = new Matrix(1, 1);
        res.XN.set(0, 0, Xbest);
        res.AN = new Matrix(0, 0);
        res.WN = new Matrix(0, 0);
        res.logNormConstAggr = -N * FastMath.log(Xbest);
        res.runtime = (System.nanoTime() - tstart) / 1000000000.0;
        res.iter = 1;
        res.method = method;
        return res;
    }

    // ---- SolverBA single-class hierarchical bound helpers ----

    private static boolean isScHierMethod(String m) {
        return m.startsWith("pbh") || m.startsWith("cbh") || m.startsWith("pbk")
                || m.startsWith("bjbk") || m.startsWith("ssd") || m.startsWith("sib")
                || m.startsWith("scb") || m.startsWith("ldbcmp");
    }

    // Aggregate think time: sum over infinite-server (delay) stations of V/rate.
    private static double baThink(NetworkStruct sn, Matrix V) {
        double Z = 0.0;
        for (int i = 0; i < V.getNumRows(); i++) {
            if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                Z += V.get(i) / sn.rates.get(i);
            }
        }
        return Z;
    }

    // Per-queue demand column (non-delay stations): V/rate.
    private static Matrix baQueueDemands(NetworkStruct sn, Matrix V) {
        int n = 0;
        for (int i = 0; i < V.getNumRows(); i++) {
            if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF) n++;
        }
        Matrix D = new Matrix(n, 1);
        int idx = 0;
        for (int i = 0; i < V.getNumRows(); i++) {
            if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF) {
                D.set(idx++, 0, V.get(i) / sn.rates.get(i));
            }
        }
        return D;
    }

    // Per-queue server-count column (non-delay stations), for ssd.
    private static Matrix baQueueServers(NetworkStruct sn, int Kq) {
        Matrix c = new Matrix(Kq, 1);
        int idx = 0;
        for (int i = 0; i < sn.nservers.getNumRows(); i++) {
            if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF) {
                c.set(idx++, 0, sn.nservers.get(i));
            }
        }
        return c;
    }

    private static void baRejectMultiserver(NetworkStruct sn, String method) {
        for (int i = 0; i < sn.nservers.getNumRows(); i++) {
            if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF && sn.nservers.get(i) > 1) {
                throw new RuntimeException("Method '" + method + "' does not support multi-server stations (use 'ssd').");
            }
        }
    }

    // Fill per-station {Q,U,R,T,C} from a scalar chain-throughput bound X,
    // using the ABA optimistic (isUpper) / pessimistic residence construction.
    private static Matrix[] baFill(NetworkStruct sn, Matrix V, double N, double X, double Zt, Matrix D, boolean isUpper) {
        int M = V.getNumRows();
        Matrix TN = new Matrix(M, 1);
        Matrix RN = new Matrix(M, 1);
        Matrix QN = new Matrix(M, 1);
        Matrix UN = new Matrix(M, 1);
        for (int i = 0; i < M; i++) {
            boolean isINF = sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF;
            double rate = sn.rates.get(i);
            double t = V.get(i) * X;
            TN.set(i, 0, t);
            double r;
            if (isUpper) {
                r = isINF ? (1.0 / rate) : (N / rate);
            } else {
                r = 1.0 / rate;
            }
            RN.set(i, 0, r);
            QN.set(i, 0, t * r);
            // Utilization law per SERVER: without the nservers divisor a
            // multiserver station reports U > 1 (ssd/ldbcmp/auto reach here)
            double c = Math.max(1.0, sn.nservers.get(i));
            UN.set(i, 0, isINF ? (t * r) : (t / (c * rate)));
        }
        Matrix CN = new Matrix(1, 1);
        CN.set(0, 0, isUpper ? (Zt + N * D.elementSum()) : (Zt + D.elementSum()));
        return new Matrix[]{QN, UN, RN, TN, CN};
    }
}
