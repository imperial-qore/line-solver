package jline.solvers.ag.handlers;

import jline.solvers.mam.handlers.MetricsResult;

import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.util.matrix.Matrix;

import java.util.List;

public final class Solver_ag_metrics {
    private Solver_ag_metrics() {}

    public static MetricsResult solver_ag_metrics(NetworkStruct sn, INAPResult inapResult, RCATModel rcat) {
        int M = sn.nstations;
        int K = sn.nclasses;
        Matrix processMap = rcat.processMap;
        List<Matrix> pi = inapResult.pi;
        int[] N = rcat.N;

        Matrix QN = new Matrix(M, K);
        Matrix UN = new Matrix(M, K);
        Matrix RN = new Matrix(M, K);
        Matrix TN = new Matrix(M, K);

        // Per-process geometric-tail decay, set by the matrix-geometric
        // 'inapinf' method; null for finite-state inap/inapplus.
        double[] rhoProc = inapResult.rhoProc;
        boolean[] isGeomProc = inapResult.isGeomProc;
        QbdTail[] geomData = inapResult.geomData;

        for (int ist = 0; ist < M; ist++) {
            for (int r = 0; r < K; r++) {
                int p = (int) processMap.get(ist, r);
                if (p >= 0 && p < pi.size() && pi.get(p).getNumCols() > 0) {
                    int mph = rcat.mph[p];

                    if (isGeomProc != null && p < isGeomProc.length && isGeomProc[p]) {
                        if (mph == 1) {
                            // Infinite geometric marginal pi_n = (1-rho) rho^n:
                            //   E[N] = rho/(1-rho),  P(N>0) = rho.
                            double rho = rhoProc[p];
                            double muIr = sn.rates.get(ist, r);
                            QN.set(ist, r, rho / (1.0 - rho));
                            UN.set(ist, r, rho);
                            if (!Double.isNaN(muIr) && muIr > 0) {
                                TN.set(ist, r, muIr * rho);
                            }
                        } else {
                            // Matrix-geometric tail pi_(n+1) = pi_n R.
                            QbdTail g = geomData[p];
                            double busySum = 0.0;
                            double tput = 0.0;
                            for (int i = 0; i < mph; i++) {
                                busySum += g.busy.get(0, i);
                                tput += g.busy.get(0, i) * rcat.svcdown[p][i];
                            }
                            QN.set(ist, r, g.qlen);
                            UN.set(ist, r, busySum);
                            TN.set(ist, r, tput);
                        }
                    } else {
                        int Np = N[p];
                        Matrix piP = pi.get(p);

                        // Queue length: E[N] = sum over states of level*pi
                        double queueLength = 0.0;
                        // Throughput: the rate of service completions, i.e. the
                        // phase-dependent departure rate averaged over the
                        // marginal. With one phase this is the mean rate times
                        // P(N>0).
                        double tput = 0.0;
                        for (int n = 0; n < Np; n++) {
                            queueLength += rcat.level[p][n] * piP.get(0, n);
                            tput += rcat.svcrate[p][n] * piP.get(0, n);
                        }
                        QN.set(ist, r, queueLength);

                        // Utilization: P(N > 0) = 1 - P(level 0)
                        double level0 = 0.0;
                        for (int i = 0; i < mph; i++) level0 += piP.get(0, i);
                        UN.set(ist, r, 1.0 - level0);
                        TN.set(ist, r, tput);
                    }
                }
            }
        }

        for (int ist = 0; ist < M; ist++) {
            for (int r = 0; r < K; r++) {
                double tput = TN.get(ist, r);
                if (tput > 0) {
                    RN.set(ist, r, QN.get(ist, r) / tput);
                } else {
                    RN.set(ist, r, 0.0);
                }
            }
        }

        Matrix CN = new Matrix(1, K);
        Matrix XN = new Matrix(1, K);

        for (int r = 0; r < K; r++) {
            double njobs = sn.njobs.get(r);

            if (Double.isInfinite(njobs)) {
                for (int ist = 0; ist < M; ist++) {
                    int nodeIdx = (int) sn.stationToNode.get(ist);
                    if (sn.nodetype.get(nodeIdx) == NodeType.Source) {
                        double rate = sn.rates.get(ist, r);
                        if (!Double.isNaN(rate)) {
                            XN.set(0, r, rate);
                        }
                        break;
                    }
                }
                double totalRespTime = 0.0;
                for (int ist = 0; ist < M; ist++) {
                    totalRespTime += RN.get(ist, r);
                }
                CN.set(0, r, totalRespTime);
            } else {
                int refst = (int) sn.refstat.get(r);
                if (refst >= 0 && refst < M) {
                    double xr = TN.get(refst, r);
                    XN.set(0, r, xr);
                    if (xr > 0) {
                        CN.set(0, r, njobs / xr);
                    }
                }
            }
        }

        return new MetricsResult(QN, UN, RN, TN, CN, XN);
    }
}
