/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.nc.handlers;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import jline.io.Ret;
import jline.api.sn.SnDeaggregateChainResults;
import jline.api.sn.SnGetDemandsChain;
import jline.util.SerializableFunction;
import jline.api.pfqn.ld.CdPeakScaling;
import jline.api.pfqn.ld.Pfqn_conv;
import jline.io.InputOutput;
import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.nodes.Station;
import jline.solvers.SolverOptions;
import jline.solvers.nc.SolverNC;
import jline.util.Utils;
import jline.util.matrix.Matrix;

/**
 * Exact normalizing constant solver for closed networks with Limited
 * class-dependent (cdscaling) service rates, using the multichain
 * convolution algorithm of Sauer (1983), Section 5.2.
 *
 * Ports MATLAB solver_nc_conv.m
 */
public final class Solver_nc_conv {
    private Solver_nc_conv() {}

    public static SolverNC.SolverNCLDReturn solver_nc_conv(NetworkStruct sn, SolverOptions options) {
        long startTime = System.nanoTime();
        final String method = "conv";
        final int iter = 1;

        int M = sn.nstations;
        Matrix nservers = sn.nservers;

        // The convolution runs on CHAINS, not classes: a class-switching model
        // splits one circulating population across several classes, so sn.njobs
        // has zeros in the classes that hold no reference jobs and the
        // class-level recursion charges those stations nothing at all. Every
        // other closed NC and MVA path aggregates the same way and deaggregates
        // at the end.
        Ret.snGetDemands chainReturn = SnGetDemandsChain.snGetDemandsChain(sn);
        Matrix Lchain = chainReturn.Dchain;
        Matrix STchain = chainReturn.STchain;
        Matrix Vchain = chainReturn.Vchain;
        Matrix alpha = chainReturn.alpha;
        Matrix Nchain = chainReturn.Nchain;
        int K = sn.nchains;
        Matrix V = Vchain;
        Matrix ST = STchain;
        Matrix Ldemand = Lchain;

        // Separate delay and queue stations
        boolean[] isDelay = new boolean[M];
        for (int i = 0; i < M; i++) {
            isDelay[i] = Utils.isInf(nservers.get(i));
        }
        List<Integer> delayIdx = new ArrayList<Integer>();
        List<Integer> queueIdx = new ArrayList<Integer>();
        for (int i = 0; i < M; i++) {
            if (isDelay[i]) {
                delayIdx.add(Integer.valueOf(i));
            } else {
                queueIdx.add(Integer.valueOf(i));
            }
        }
        int nQueues = queueIdx.size();

        // Z_conv: total delay demand per class
        double[] Z_conv = new double[K];
        for (Integer ist : delayIdx) {
            for (int k = 0; k < K; k++) {
                Z_conv[k] += Ldemand.get(ist.intValue(), k);
            }
        }

        // L_conv: demands for queue stations only
        Matrix L_conv = new Matrix(nQueues, K);
        for (int qi = 0; qi < nQueues; qi++) {
            for (int k = 0; k < K; k++) {
                L_conv.set(qi, k, Ldemand.get(queueIdx.get(qi).intValue(), k));
            }
        }

        // NK: population vector, per chain
        int[] NK = new int[K];
        for (int k = 0; k < K; k++) {
            NK[k] = (int) Math.round(Nchain.get(k));
        }
        // Class-dependence functions beta_{i,r}(n) for the queue stations. Each
        // takes the per-class population vector at its station and returns either
        // a scalar (chain-independent) or a length-R vector of per-class rates.
        List<Station> stations = sn.stations;
        List<JobClass> jobClasses = sn.jobclasses;
        // Joint-dependence handles (sn.jdscaling, non-product-form eta_i) are
        // folded into the same per-station handle used by the convolution
        // recursion. cd and jd are evaluated identically; the product reproduces
        // the single-mechanism case when only one is present.
        boolean hasCd = sn.cdscaling != null && !sn.cdscaling.isEmpty();
        boolean hasJd = sn.jdscaling != null && !sn.jdscaling.isEmpty();
        List<SerializableFunction<Matrix, Matrix>> cdscalingConv = null;
        if (hasCd || hasJd) {
            cdscalingConv = new ArrayList<SerializableFunction<Matrix, Matrix>>();
            for (int qi = 0; qi < nQueues; qi++) {
                Station station = stations.get(queueIdx.get(qi).intValue());
                SerializableFunction<Matrix, Matrix> cdh = hasCd ? sn.cdscaling.get(station) : null;
                SerializableFunction<Matrix, Matrix> jdh = hasJd ? sn.jdscaling.get(station) : null;
                cdscalingConv.add(combineCdJd(cdh, jdh));
            }
        }

        // Compute G(N)
        double[] result = Pfqn_conv.pfqn_conv(L_conv, NK, Z_conv, cdscalingConv);
        double G_N = result[0];
        double lG = result[1];

        // Compute G(N - e_k) for each class -> throughput
        double[] XN = new double[K];
        for (int k = 0; k < K; k++) {
            if (NK[k] > 0) {
                int[] NKminus = NK.clone();
                NKminus[k]--;
                double[] resultK = Pfqn_conv.pfqn_conv(L_conv, NKminus, Z_conv, cdscalingConv);
                XN[k] = resultK[0] / G_N;
            }
        }

        // Per-station throughput TN = V * XN
        Matrix TN = new Matrix(M, K);
        for (int i = 0; i < M; i++) {
            for (int k = 0; k < K; k++) {
                TN.set(i, k, V.get(i, k) * XN[k]);
            }
        }

        // Compute queue lengths QN
        Matrix QN = new Matrix(M, K);

        // Delay stations: Q = L * X
        for (Integer ist : delayIdx) {
            for (int k = 0; k < K; k++) {
                QN.set(ist.intValue(), k, Ldemand.get(ist.intValue(), k) * XN[k]);
            }
        }

        // Queue stations: marginal distribution
        int stateSpaceSize = 1;
        for (int n : NK) stateSpaceSize *= (n + 1);

        for (int qi = 0; qi < nQueues; qi++) {
            int ist = queueIdx.get(qi).intValue();

            double[] Xm = new double[stateSpaceSize];
            Xm[0] = 1.0;
            boolean isCdStation = cdscalingConv != null && cdscalingConv.get(qi) != null;

            int[] n = new int[K];
            while (true) {
                int idx = hashpopConv(n, NK);
                int sumN = 0;
                for (int v : n) sumN += v;
                if (sumN > 0) {
                    if (isCdStation) {
                        for (int r = 0; r < K; r++) {
                            if (n[r] > 0) {
                                // beta_{qi,r}(n): DIMENSIONLESS scaling of the demand,
                                // so the effective demand is L/beta. The function
                                // returns a scalar (shared by all classes) or a
                                // length-R vector.
                                Matrix nvec = new Matrix(1, K);
                                int tot = 0;
                                for (int k = 0; k < K; k++) {
                                    nvec.set(0, k, (double) n[k]);
                                    tot += n[k];
                                }
                                Matrix bval = cdscalingConv.get(qi).apply(nvec);
                                double beta = (bval.length() > 1) ? bval.get(r) : bval.get(0);
                                // X_m(n) = (|n|/n_r) * (L/beta) * X_m(n-e_r); at
                                // beta=1 this is the LI multinomial recurrence.
                                double nr = (double) n[r];
                                n[r]--;
                                int idxPrev = hashpopConv(n, NK);
                                n[r]++;
                                if (beta > 0) {
                                    Xm[idx] = ((double) tot / nr) * (L_conv.get(qi, r) / beta) * Xm[idxPrev];
                                }
                                break;
                            }
                        }
                    } else {
                        for (int r = 0; r < K; r++) {
                            if (n[r] > 0) {
                                n[r]--;
                                int idxPrev = hashpopConv(n, NK);
                                n[r]++;
                                Xm[idx] += L_conv.get(qi, r) * Xm[idxPrev];
                            }
                        }
                    }
                }
                if (!pprodNextConv(n, NK)) break;
            }

            // Build complement: all stations except qi
            Matrix L_comp = new Matrix(nQueues - 1, K);
            List<SerializableFunction<Matrix, Matrix>> cdComp = new ArrayList<SerializableFunction<Matrix, Matrix>>();
            int ci = 0;
            for (int qj = 0; qj < nQueues; qj++) {
                if (qj != qi) {
                    for (int k = 0; k < K; k++) {
                        L_comp.set(ci, k, L_conv.get(qj, k));
                    }
                    if (cdscalingConv != null && qj < cdscalingConv.size()) {
                        cdComp.add(cdscalingConv.get(qj));
                    } else {
                        cdComp.add(null);
                    }
                    ci++;
                }
            }

            int[] nMarg = new int[K];
            while (true) {
                boolean anyPos = false;
                for (int v : nMarg) {
                    if (v > 0) { anyPos = true; break; }
                }
                if (anyPos) {
                    int idx = hashpopConv(nMarg, NK);
                    int[] nmi = new int[K];
                    boolean allNonNeg = true;
                    for (int k = 0; k < K; k++) {
                        nmi[k] = NK[k] - nMarg[k];
                        if (nmi[k] < 0) { allNonNeg = false; break; }
                    }
                    if (allNonNeg) {
                        double[] resultComp = Pfqn_conv.pfqn_conv(
                                L_comp, nmi, Z_conv,
                                cdComp.isEmpty() ? null : cdComp);
                        double prob = Xm[idx] * resultComp[0] / G_N;
                        for (int k = 0; k < K; k++) {
                            QN.set(ist, k, QN.get(ist, k) + nMarg[k] * prob);
                        }
                    }
                }
                if (!pprodNextConv(nMarg, NK)) break;
            }
        }

        // Remaining metrics. RN is the PER-VISIT response time Qchain/Tchain: the
        // deaggregation below multiplies the visit ratio back in, so dividing by
        // Xchain would count it twice
        Matrix RN = new Matrix(M, K);
        for (int i = 0; i < M; i++) {
            for (int k = 0; k < K; k++) {
                RN.set(i, k, (TN.get(i, k) != 0.0) ? QN.get(i, k) / TN.get(i, k) : 0.0);
            }
        }

        Matrix UN = new Matrix(M, K);
        for (int i = 0; i < M; i++) {
            for (int k = 0; k < K; k++) {
                UN.set(i, k, TN.get(i, k) * ST.get(i, k));
            }
        }

        // Utilization at a class-dependent station is the fraction of the
        // station's PEAK service capacity in use, not T*ST: the scaling
        // beta_{i,r}(n) multiplies the nominal rate, so T*ST measures capacity
        // used in units of the nominal rate and reaches max_n beta(n), not 1,
        // at saturation. Normalizing by that peak matches what Solver_ncld
        // does with max(lldscaling(ist,:)) for the load-dependent case. On a
        // beta emulating c servers this returns E[busy]/c, matching the true
        // c-server station; without it a 2-server beta reports U = 2*(true
        // utilization) and exceeds 1.
        //
        // Utilization is normalized by the declared per-class peak rate scaling
        // (sn.cdscalingpeak), giving the T*S/peak convention of multiserver
        // stations (mirrors MATLAB solver_nc_conv.m).
        if (cdscalingConv != null) {
            for (int qi = 0; qi < nQueues; qi++) {
                SerializableFunction<Matrix, Matrix> beta = cdscalingConv.get(qi);
                if (beta == null) continue;
                int ist = queueIdx.get(qi).intValue();
                Station st = sn.stations.get(ist);
                // Effective peak = product of the class- and joint-dependence
                // peaks declared at the station (a missing one contributes 1).
                Matrix cdPeakVec = (hasCd && sn.cdscalingpeak != null) ? sn.cdscalingpeak.get(st) : null;
                Matrix jdPeakVec = (hasJd && sn.jdscalingpeak != null) ? sn.jdscalingpeak.get(st) : null;
                for (int c = 0; c < K; c++) {
                    // The peaks are declared per class, so the chain takes the
                    // largest peak among its classes: utilization is a
                    // per-station quantity with one normalizer.
                    double bmax = 1.0;
                    if (cdPeakVec != null) bmax *= chainPeak(cdPeakVec, sn.chains, c);
                    if (jdPeakVec != null) bmax *= chainPeak(jdPeakVec, sn.chains, c);
                    if (bmax > 0) {
                        UN.set(ist, c, UN.get(ist, c) / bmax);
                    }
                }
            }
        }

        Matrix Xchain = new Matrix(1, K);
        for (int k = 0; k < K; k++) Xchain.set(0, k, XN[k]);

        // Deaggregate the chain solution onto the classes
        Ret.snDeaggregateChainResults deagg = SnDeaggregateChainResults.snDeaggregateChainResults(
                sn, Lchain, null, STchain, Vchain, alpha, null, UN, RN, TN, null, Xchain);

        double runtime = (System.nanoTime() - startTime) / 1e9;

        return new SolverNC.SolverNCLDReturn(deagg.Q, deagg.U, deagg.R, deagg.T, deagg.C, deagg.X,
                lG, runtime, iter, method);
    }

    /** Largest declared peak among the classes of chain c, the station's single normalizer. */
    private static double chainPeak(Matrix peakPerClass, Matrix chains, int c) {
        double bmax = 0.0;
        for (int r = 0; r < chains.getNumCols(); r++) {
            if (chains.get(c, r) > 0) {
                bmax = Math.max(bmax, peakPerClass.get(0, r));
            }
        }
        return bmax > 0 ? bmax : 1.0;
    }

    // The class-dependence lattice peak now lives in
    // jline.api.pfqn.ld.CdPeakScaling (shared by NC, CTMC and MVA so all
    // solvers report class-dependent utilization under one convention).

    private static int hashpopConv(int[] n, int[] N) {
        int idx = 0;
        int stride = 1;
        for (int r = 0; r < n.length; r++) {
            idx += stride * n[r];
            stride *= (N[r] + 1);
        }
        return idx;
    }

    private static boolean pprodNextConv(int[] n, int[] N) {
        int R = n.length;
        boolean atMax = true;
        for (int i = 0; i < R; i++) {
            if (n[i] != N[i]) { atMax = false; break; }
        }
        if (atMax) return false;
        int s = R - 1;
        while (s >= 0 && n[s] == N[s]) {
            n[s] = 0;
            s--;
        }
        if (s >= 0) n[s]++;
        return true;
    }

    /**
     * Combine a class-dependence handle beta and a joint-dependence handle eta
     * into a single per-station handle whose value is the element-wise product
     * beta(ni).*eta(ni), broadcasting a 1x1 result against a 1xR result (as
     * MATLAB's scalar .* vector does). Returns null if both are null, or the
     * single non-null handle unchanged when only one is present, so the common
     * single-mechanism case is numerically identical to the pre-jd behaviour.
     */
    private static SerializableFunction<Matrix, Matrix> combineCdJd(
            final SerializableFunction<Matrix, Matrix> beta,
            final SerializableFunction<Matrix, Matrix> eta) {
        if (beta == null) return eta;
        if (eta == null) return beta;
        return new SerializableFunction<Matrix, Matrix>() {
            public Matrix apply(Matrix ni) {
                Matrix a = beta.apply(ni);
                Matrix b = eta.apply(ni);
                int la = a.length();
                int lb = b.length();
                int n = Math.max(la, lb);
                Matrix out = new Matrix(1, n);
                for (int j = 0; j < n; j++) {
                    double av = (la == 1) ? a.get(0) : a.get(j);
                    double bv = (lb == 1) ? b.get(0) : b.get(j);
                    out.set(0, j, av * bv);
                }
                return out;
            }
        };
    }
}
