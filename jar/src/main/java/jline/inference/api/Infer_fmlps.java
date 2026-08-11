/*
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
package jline.inference.api;

import java.util.ArrayList;
import java.util.List;
import java.util.TreeSet;

import jline.util.Pair;
import java.util.function.Function;

import jline.inference.util.OptimUtils;
import jline.lang.JobClass;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Station;
import jline.util.matrix.Matrix;

public final class Infer_fmlps {
    private Infer_fmlps() {}

    /**
     * FMLPS demand estimation using fluid-based likelihood.
     */
    public static double[] infer_fmlps(Network model,
                                       Queue node,
                                       final double[] rt,
                                       final int[] classVec,
                                       final Matrix ql,
                                       final int W) {
        final NetworkStruct sn = model.getStruct(false);
        final int R = sn.nclasses;
        node.getNumberOfServers(); // V (unused but mirrors original)
        final int stIdx = node.getStationIdx();

        double rtMin = Double.POSITIVE_INFINITY;
        double rtMax = Double.NEGATIVE_INFINITY;
        for (double v : rt) {
            if (v < rtMin) rtMin = v;
            if (v > rtMax) rtMax = v;
        }
        final double[] xLB = new double[R];
        final double[] xUB = new double[R];
        for (int i = 0; i < R; i++) {
            xLB[i] = rtMin / W;
            xUB[i] = rtMax;
        }

        double meanQLsum = 0.0;
        for (int i = 0; i < ql.getNumRows(); i++) {
            double rowSum = 0.0;
            for (int j = 0; j < ql.getNumCols(); j++) {
                rowSum += ql.get(i, j);
            }
            meanQLsum += rowSum;
        }
        final double meanQL = (ql.getNumRows() > 0) ? meanQLsum / ql.getNumRows() : 0.0;
        int Vsrv = node.getNumberOfServers();
        final double Vtilde = Math.min(meanQL, (double) Vsrv);
        double[] x0 = new double[R];
        for (int j = 0; j < R; j++) {
            double sum = 0.0;
            int count = 0;
            for (int i = 0; i < rt.length; i++) {
                if (classVec[i] == j) {
                    sum += rt[i];
                    count++;
                }
            }
            if (count > 0 && meanQL > 0) {
                x0[j] = Vtilde * (sum / count) / meanQL;
            } else {
                x0[j] = xLB[j];
            }
        }

        int delayIdxFound = -1;
        int refIdxFound = -1;
        for (int ii = 0; ii < sn.nstations; ii++) {
            Station station = sn.stations.get(ii);
            if (sn.sched.get(station) == SchedStrategy.INF) {
                delayIdxFound = ii;
            } else {
                refIdxFound = ii;
            }
        }
        final int delayIdx = delayIdxFound;
        final int refIdx = refIdxFound;

        final double[] delayRates = new double[R];
        for (int k = 0; k < R; k++) {
            Station station = sn.stations.get(delayIdx);
            JobClass jobclass = sn.jobclasses.get(k);
            delayRates[k] = sn.mu.get(station).get(jobclass).get(0, 0);
        }

        Function<double[], Double> objFun = new Function<double[], Double>() {
            @Override
            public Double apply(double[] x) {
                final double TOL = 1e-6;

                for (int r = 0; r < R; r++) {
                    Sn_set_service_coc.sn_set_service_coc(sn, stIdx, r, 1.0 / x[r]);
                }

                TreeSet<Integer> uniqueTC = new TreeSet<Integer>();
                for (int v : classVec) uniqueTC.add(v);

                double[] ftemp = new double[rt.length];

                for (Integer tc : uniqueTC) {
                    List<Integer> mask = new ArrayList<Integer>();
                    for (int i = 0; i < classVec.length; i++) {
                        if (classVec[i] == tc.intValue()) mask.add(i);
                    }
                    double[] rtTc = new double[mask.size()];
                    Matrix qlTc = new Matrix(mask.size(), R);
                    for (int i = 0; i < mask.size(); i++) {
                        rtTc[i] = rt[mask.get(i)];
                        for (int r = 0; r < R; r++) {
                            qlTc.set(i, r, ql.get(mask.get(i), r));
                        }
                    }

                    Infer_fluid_ps_rt_likelihood.FluidPsRtResult fluidResult =
                            Infer_fluid_ps_rt_likelihood.infer_fluid_ps_rt_likelihood(sn, tc.intValue());

                    for (int rr = 0; rr < rtTc.length; rr++) {
                        double[] aQueue = new double[R];
                        for (int r = 0; r < R; r++) {
                            aQueue[r] = qlTc.get(rr, r);
                        }

                        int M = fluidResult.augPhases.getNumRows();
                        Matrix y0Levels = new Matrix(M, R);

                        double totalDelayRate = 0.0;
                        for (double dr : delayRates) totalDelayRate += dr;
                        double aQueueSum = 0.0;
                        for (double q : aQueue) aQueueSum += q;
                        double delayJobsTotal = W - aQueueSum;
                        for (int k = 0; k < R; k++) {
                            double delayJob = (totalDelayRate > 0) ? delayJobsTotal * delayRates[k] / totalDelayRate : 0.0;
                            y0Levels.set(delayIdx, k, delayJob);
                        }

                        for (int k = 0; k < R; k++) {
                            y0Levels.set(refIdx, k, aQueue[k]);
                        }

                        double like = Infer_fluid_ps_rt_likelihood.infer_fluid_ps_rt_solve(
                                fluidResult, y0Levels, rtTc[rr], tc.intValue());
                        ftemp[mask.get(rr)] = Math.log(TOL + like);
                    }
                }
                double sum = 0.0;
                for (double v : ftemp) sum += v;
                return -sum;
            }
        };

        Pair<double[], Double> result = OptimUtils.fmincon(objFun, x0, xLB, xUB, 10000);
        return result.getFirst();
    }
}
