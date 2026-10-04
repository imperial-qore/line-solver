package jline.api.sn;

import java.util.Map;

import jline.GlobalConstants;
import jline.lang.JobClass;
import jline.lang.Metric;
import jline.lang.NetworkStruct;
import jline.lang.nodes.Station;
import jline.solvers.AvgHandle;
import jline.util.matrix.Matrix;

public final class SnGetResidTFromRespT {
    private SnGetResidTFromRespT() {}

    /**
     * Calculates the residence times at each station from the response times.
     *
     * @param sn      the NetworkStruct object for the queueing network model
     * @param RNclass the matrix of resp times
     * @param WH      residt handles for the solver
     * @return a matrix of residence times
     */
    public static Matrix snGetResidTFromRespT(NetworkStruct sn, Matrix RNclass, AvgHandle WH) {
        int M = sn.nstations;
        int K = sn.nclasses;
        int nstateful = sn.nstateful;
        int Vcells = sn.visits.size();

        Matrix V = new Matrix(M, K);
        V.zero();

        for (int chainIdx = 0; chainIdx < Vcells; chainIdx++) {
            Matrix visits = sn.visits.get(chainIdx);
            if (visits == null) continue;
            for (int sf = 0; sf < nstateful; sf++) {
                double stationIdxRaw;
                if (sn.statefulToStation != null && sf < sn.statefulToStation.getNumCols()) {
                    stationIdxRaw = sn.statefulToStation.get(0, sf);
                } else {
                    stationIdxRaw = Double.NaN;
                }
                if (!Double.isNaN(stationIdxRaw)) {
                    int stationIdx = (int) stationIdxRaw;
                    if (stationIdx >= 0 && stationIdx < M) {
                        for (int r = 0; r < K; r++) {
                            if (r < visits.getNumCols()) {
                                V.set(stationIdx, r, V.get(stationIdx, r) + visits.get(sf, r));
                            }
                        }
                    }
                }
            }
        }

        Matrix WNclass = RNclass.copy();
        WNclass.zero();

        for (int ist = 0; ist < M; ist++) {
            for (int k = 0; k < K; k++) {
                boolean isHandleDisabled = true;
                if (WH != null && !WH.isEmpty()) {
                    Station station = sn.stations.get(ist);
                    Map<JobClass, Metric> innerMap = WH.get(station);
                    if (innerMap != null) {
                        Metric metric = innerMap.get(sn.jobclasses.get(k));
                        if (metric != null) {
                            isHandleDisabled = metric.isDisabled;
                        }
                    }
                }

                if (isHandleDisabled) {
                    WNclass.set(ist, k, Double.NaN);
                } else if (!RNclass.isEmpty() && RNclass.get(ist, k) > 0) {
                    int c = -1;
                    for (int chain = 0; chain < sn.chains.getNumRows(); chain++) {
                        if (sn.chains.get(chain, k) > 0) {
                            c = chain;
                            break;
                        }
                    }
                    if (RNclass.get(ist, k) < GlobalConstants.FineTol) {
                        WNclass.set(ist, k, RNclass.get(ist, k));
                    } else {
                        int refClass = (int) sn.refclass.get(0, c);
                        if (refClass >= 0) {
                            double denom = V.get((int) sn.refstat.get(k, 0), refClass);
                            WNclass.set(ist, k, RNclass.get(ist, k) * V.get(ist, k) / denom);
                        } else {
                            int Vrow = (int) sn.refstat.get(k, 0);
                            double Vsum = 0.0;
                            Matrix inchainC = sn.inchain.get(c);
                            for (int col = 0; col < inchainC.getNumCols(); col++) {
                                int classIdx = (int) inchainC.get(0, col);
                                Vsum += V.get(Vrow, classIdx);
                            }
                            WNclass.set(ist, k, RNclass.get(ist, k) * V.get(ist, k) / Vsum);
                        }
                    }
                }
            }
        }
        for (int ist = 0; ist < M; ist++) {
            for (int k = 0; k < K; k++) {
                double v = WNclass.get(ist, k);
                if (Double.isNaN(v) || v < 10 * GlobalConstants.FineTol || v < GlobalConstants.Zero) {
                    WNclass.set(ist, k, 0);
                }
            }
        }
        return WNclass;
    }
}
