/**
 * Map a state distribution to mean performance metrics (CTMC reduction).
 *
 * Given an arbitrary probability vector over the enumerated CTMC state space of
 * {@code sn} (rows of StateSpace / StateSpaceAggr), returns the per-(station,
 * class) mean queue length QN, utilization UN, response time RN, throughput TN,
 * system response time CN and system throughput XN. The discipline-aware mapping
 * is identical to the steady-state reduction performed by Solver_ctmc_analyzer;
 * it is factored here so that callers holding their own distribution (e.g. the
 * SolverENV state-vector analyzer, which time-averages a transient distribution)
 * can reuse it without re-solving for the stationary vector.
 *
 * Mirrors matlab/src/solvers/CTMC/solver_ctmc_avg_from_pi.m.
 *
 * @since LINE 3.0
 */
package jline.solvers.ctmc.handlers;

import java.util.Map;

import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.nodes.Station;
import jline.lang.constant.NodeType;
import jline.lang.constant.SchedStrategy;
import jline.lang.state.State;
import jline.lang.state.ToMarginal;
import jline.api.mam.Map_mean;
import jline.api.pfqn.ld.CdPeakScaling;
import jline.io.InputOutput;
import jline.solvers.SolverOptions;
import jline.util.Maths;
import jline.util.SerializableFunction;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Ctmc_avg_from_pi {
    private Ctmc_avg_from_pi() {}

    public static final class Result {
        public Matrix QN;
        public Matrix UN;
        public Matrix RN;
        public Matrix TN;
        public Matrix CN;
        public Matrix XN;
    }

    public static Result ctmc_avg_from_pi(NetworkStruct sn, Matrix pivec, Matrix StateSpace,
                                          Matrix StateSpaceAggr, double[][][] arvRates, double[][][] depRates) {
        // see _kb/06-solver-catalog.md for rationale
        return ctmc_avg_from_pi(sn, pivec, StateSpace, StateSpaceAggr, arvRates, depRates, null);
    }

    @SuppressWarnings("unchecked")
    public static Result ctmc_avg_from_pi(NetworkStruct sn, Matrix pivec, Matrix StateSpace,
                                          Matrix StateSpaceAggr, double[][][] arvRates, double[][][] depRates,
                                          SolverOptions options) {
        int M = sn.nstations;
        int K = sn.nclasses;
        Matrix S = sn.nservers;
        Matrix NK = sn.njobs;
        Map<Station, SchedStrategy> schedid = sn.sched;
        Object PH = sn.proc;

        // Normalise the supplied distribution (drop sub-zero noise).
        Matrix probSysState = pivec.copy();
        if (probSysState.getNumCols() == 1 && probSysState.getNumRows() > 1) {
            probSysState = probSysState.transpose();
        }
        for (int col = 0; col < probSysState.getNumCols(); col++) {
            if (probSysState.get(col) < jline.GlobalConstants.Zero) {
                probSysState.set(col, 0);
            }
        }
        double psum = probSysState.elementSum();
        if (psum > 0) {
            probSysState.divideEq(psum);
        }

        int n = StateSpace.getNumRows();
        Matrix wset = new Matrix(1, n);
        for (int col = 0; col < n; col++) {
            wset.set(0, col, col);
        }
        Matrix StateSpaceWork = StateSpace;
        Matrix StateSpaceAggrWork = StateSpaceAggr;

        Matrix XN = new Matrix(1, K); XN.zero();
        Matrix UN = new Matrix(M, K); UN.zero();
        Matrix QN = new Matrix(M, K); QN.zero();
        Matrix RN = new Matrix(M, K); RN.zero();
        Matrix TN = new Matrix(M, K); TN.zero();
        Matrix CN = new Matrix(1, K); CN.zero();

        Matrix istSpaceShift = new Matrix(1, M);
        istSpaceShift.zero();
        for (int i = 0; i < M; i++) {
            if (i == 0) {
                istSpaceShift.set(0, i, 0);
            } else {
                double temp = istSpaceShift.get(0, i - 1) + sn.space.get(sn.stateful.get(i - 1)).getNumCols();
                istSpaceShift.set(0, i, temp);
            }
        }

        for (int k = 0; k < K; k++) {
            int refsf = (int) sn.stationToStateful.get((int) sn.refstat.get(k));
            double sumValue = 0.0;
            for (int i = 0; i < wset.getNumCols(); i++) {
                int index = (int) wset.get(i);
                sumValue += probSysState.get(index) * arvRates[index][refsf][k];
            }
            XN.set(0, k, sumValue);
        }

        for (int i = 0; i < M; i++) {
            int isf = (int) sn.stationToStateful.get(i);
            int ind = (int) sn.stationToNode.get(i);
            for (int k = 0; k < K; k++) {
                double sumTN = 0.0;
                double sumQN = 0.0;
                for (int index = 0; index < wset.getNumCols(); index++) {
                    int wstIdx = (int) wset.get(index);
                    double depRate = depRates[wstIdx][isf][k];
                    double probState = probSysState.get(wstIdx);
                    sumTN += probState * depRate;
                    sumQN += probState * StateSpaceAggrWork.get(wstIdx, i * K + k);
                }
                TN.set(i, k, sumTN);
                QN.set(i, k, sumQN);
            }
            if (sn.nodetype.get(ind) != NodeType.Source) {
                // see _kb/06-solver-catalog.md for rationale
                boolean stationCapFinite = !Double.isInfinite(sn.cap.get(i)) && sn.cap.get(i) < Integer.MAX_VALUE;
                boolean[] canDropClass = new boolean[K];
                for (int r = 0; r < K; r++) {
                    boolean classCapFinite = !Double.isInfinite(sn.classcap.get(i, r)) && sn.classcap.get(i, r) < Integer.MAX_VALUE;
                    canDropClass[r] = Double.isInfinite(sn.njobs.get(r)) && (stationCapFinite || classCapFinite);
                }
                // see _kb/06-solver-catalog.md for rationale
                double[] arvAtStation = new double[K];
                for (int r = 0; r < K; r++) {
                    for (int idx = 0; idx < wset.length(); idx++) {
                        arvAtStation[r] += probSysState.get(idx) * arvRates[(int) wset.get(idx)][isf][r];
                    }
                }
                boolean[] signalLossy = CtmcSignalLossy.signalLossyClasses(sn, arvAtStation);
                for (int r = 0; r < K; r++) {
                    canDropClass[r] = canDropClass[r] || signalLossy[r];
                }
                SchedStrategy schedStrategy = schedid.get(sn.stations.get(i));
                if (schedStrategy == SchedStrategy.INF) {
                    for (int k = 0; k < K; k++) {
                        UN.set(i, k, QN.get(i, k));
                    }
                } else if (schedStrategy == SchedStrategy.PS || schedStrategy == SchedStrategy.DPS || schedStrategy == SchedStrategy.GPS) {
                    if (sn.lldscaling.isEmpty() && sn.cdscaling.isEmpty() && (sn.jdscaling == null || sn.jdscaling.isEmpty())) {
                        for (int k = 0; k < K; k++) {
                            MatrixCell value = ((Map<Object, Map<Object, MatrixCell>>) PH).get(sn.stations.get(i)).get(sn.jobclasses.get(k));
                            if (!value.isEmpty()) {
                                double mean = Map_mean.map_mean(value) / S.get(i);
                                double UNarv_ik = 0.0;
                                for (int idx = 0; idx < wset.length(); idx++) {
                                    UNarv_ik += probSysState.get(idx) * arvRates[(int) wset.get(idx)][isf][k];
                                }
                                UNarv_ik = UNarv_ik * mean;
                                double UNdep_ik = TN.get(i, k) * mean;
                                UN.set(i, k, canDropClass[k] ? UNdep_ik : Maths.max(UNarv_ik, UNdep_ik));
                            }
                        }
                    } else {
                        // see _kb/06-solver-catalog.md for rationale
                        double ceff = S.get(i);
                        if (sn.lldscaling != null && !sn.lldscaling.isEmpty()
                                && i < sn.lldscaling.getNumRows()) {
                            for (int lldCol = 0; lldCol < sn.lldscaling.getNumCols(); lldCol++) {
                                ceff = Math.max(ceff, sn.lldscaling.get(i, lldCol));
                            }
                        }
                        for (int k = 0; k < K; k++) {
                            UN.set(i, k, 0);
                        }
                        for (int index = 0; index < wset.getNumCols(); index++) {
                            int st = (int) wset.get(index);
                            int c0 = (int) istSpaceShift.get(i);
                            int c1 = c0 + sn.space.get(sn.stateful.get(i)).getNumCols();
                            State.StateMarginalStatistics tm = ToMarginal.toMarginal(sn, ind,
                                    Matrix.extract(StateSpaceWork, st, st + 1, c0, c1), null, null, null, null, null);
                            Matrix ni = tm.ni;
                            Matrix nir = tm.nir;
                            boolean checkNi = true;
                            for (int niIndex = 0; niIndex < ni.length(); niIndex++) {
                                if (ni.get(niIndex) <= 0) checkNi = false;
                            }
                            if (checkNi) {
                                double lldnow = 1.0;
                                if (sn.lldscaling != null && !sn.lldscaling.isEmpty()
                                        && i < sn.lldscaling.getNumRows()) {
                                    double nitot = 0.0;
                                    for (int niIndex = 0; niIndex < ni.length(); niIndex++) {
                                        nitot += ni.get(niIndex);
                                    }
                                    int col = (int) Math.min(Math.max(nitot, 1), sn.lldscaling.getNumCols());
                                    lldnow = sn.lldscaling.get(i, col - 1);
                                }
                                for (int k = 0; k < K; k++) {
                                    double v = probSysState.get(st) * nir.get(k) * sn.schedparam.get(i, k);
                                    Matrix dividend = nir.mult(sn.schedparam.getRow(i).transpose());
                                    double addSum = 0.0;
                                    for (int divIndex = 0; divIndex < dividend.length(); divIndex++) {
                                        addSum += v / dividend.get(divIndex) * lldnow / ceff;
                                    }
                                    UN.set(i, k, addSum + UN.get(i, k));
                                }
                            }
                        }
                    }
                } else if (schedStrategy == SchedStrategy.PAS || schedStrategy == SchedStrategy.OI) {
                    for (int k = 0; k < K; k++) {
                        UN.set(i, k, 0);
                    }
                    for (int index = 0; index < wset.length(); index++) {
                        int st = (int) wset.get(index);
                        int c0 = (int) istSpaceShift.get(i);
                        int c1 = c0 + sn.space.get(sn.stateful.get(i)).getNumCols();
                        State.StateMarginalStatistics tm = ToMarginal.toMarginal(sn, ind,
                                Matrix.extract(StateSpaceWork, st, st + 1, c0, c1), null, null, null, null, null);
                        Matrix sir = tm.sir;
                        for (int k = 0; k < K; k++) {
                            UN.set(i, k, UN.get(i, k) + probSysState.get(st) * sir.get(k) / S.get(i));
                        }
                    }
                } else {
                    if ((sn.lldscaling == null || sn.lldscaling.isEmpty()) && (sn.cdscaling == null || sn.cdscaling.isEmpty()) && (sn.jdscaling == null || sn.jdscaling.isEmpty())) {
                        for (int k = 0; k < K; k++) {
                            MatrixCell value = ((Map<Object, Map<Object, MatrixCell>>) PH).get(sn.stations.get(i)).get(sn.jobclasses.get(k));
                            if (!value.isEmpty()) {
                                double mean = Map_mean.map_mean(value);
                                double UNarv_ik = 0.0;
                                for (int idx = 0; idx < wset.length(); idx++) {
                                    UNarv_ik += probSysState.get(idx) * arvRates[(int) wset.get(idx)][isf][k];
                                }
                                UNarv_ik = UNarv_ik * mean / S.get(i);
                                double UNdep_ik = TN.get(i, k) * mean / S.get(i);
                                UN.set(i, k, canDropClass[k] ? UNdep_ik : Maths.max(UNarv_ik, UNdep_ik));
                            }
                        }
                    } else {
                        for (int k = 0; k < K; k++) {
                            UN.set(i, k, 0);
                        }
                        for (int index = 0; index < wset.length(); index++) {
                            int st = (int) wset.get(index);
                            int c0 = (int) istSpaceShift.get(i);
                            int c1 = c0 + sn.space.get(sn.stateful.get(i)).getNumCols();
                            State.StateMarginalStatistics tm = ToMarginal.toMarginal(sn, ind,
                                    Matrix.extract(StateSpaceWork, st, st + 1, c0, c1), null, null, null, null, null);
                            Matrix ni = tm.ni;
                            Matrix sir = tm.sir;
                            boolean checkNir = true;
                            for (int niIndex = 0; niIndex < ni.length(); niIndex++) {
                                if (ni.get(niIndex) <= 0) checkNir = false;
                            }
                            if (checkNir) {
                                for (int k = 0; k < K; k++) {
                                    UN.set(i, k, UN.get(i, k) + probSysState.get(st) * sir.get(k) / S.get(i));
                                }
                            }
                        }
                    }
                }
                // see _kb/06-solver-catalog.md for rationale
                boolean anySignalLossy = false;
                for (int r = 0; r < K; r++) {
                    anySignalLossy = anySignalLossy || signalLossy[r];
                }
                if (anySignalLossy && schedStrategy != SchedStrategy.INF
                        && sn.lldscaling.isEmpty() && sn.cdscaling.isEmpty() && (sn.jdscaling == null || sn.jdscaling.isEmpty())) {
                    double[] UNb = CtmcSignalBusy.busyFraction(sn, ind, i, schedStrategy, S.get(i),
                            StateSpaceWork, istSpaceShift, wset, probSysState, K);
                    for (int k = 0; k < K; k++) {
                        if (signalLossy[k]) {
                            UN.set(i, k, UNb[k]);
                        }
                    }
                }
            }
        }

        // see _kb/06-solver-catalog.md for rationale
        if (sn.cdscaling != null && !sn.cdscaling.isEmpty()) {
            boolean allFinite = true;
            for (int k = 0; k < K; k++) {
                if (Double.isInfinite(sn.njobs.get(k))) { allFinite = false; break; }
            }
            if (allFinite) {
                for (int ist = 0; ist < M; ist++) {
                    Station stat = sn.stations.get(ist);
                    SerializableFunction<Matrix, Matrix> beta = sn.cdscaling.get(stat);
                    if (beta == null) continue;
                    Matrix peakVec = sn.cdscalingpeak != null ? sn.cdscalingpeak.get(stat) : null;
                    for (int k = 0; k < K; k++) {
                        double rate = sn.rates.get(ist, k);
                        double bmax = (peakVec != null) ? peakVec.get(0, k) : 1.0;
                        if (Double.isFinite(rate) && rate > 0 && bmax > 0) {
                            UN.set(ist, k, TN.get(ist, k) / rate / bmax);
                        } else {
                            UN.set(ist, k, 0.0);
                        }
                    }
                }
            }
        }

        // joint-dependence utilization normalization (U = T*S/peak using the
        // declared sn.jdscalingpeak), mirroring the class-dependence block.
        if (sn.jdscaling != null && !sn.jdscaling.isEmpty()) {
            boolean allFinite = true;
            for (int k = 0; k < K; k++) {
                if (Double.isInfinite(sn.njobs.get(k))) { allFinite = false; break; }
            }
            if (allFinite) {
                for (int ist = 0; ist < M; ist++) {
                    Station stat = sn.stations.get(ist);
                    SerializableFunction<Matrix, Matrix> eta = sn.jdscaling.get(stat);
                    if (eta == null) continue;
                    Matrix peakVec = sn.jdscalingpeak != null ? sn.jdscalingpeak.get(stat) : null;
                    for (int k = 0; k < K; k++) {
                        double rate = sn.rates.get(ist, k);
                        double bmax = (peakVec != null) ? peakVec.get(0, k) : 1.0;
                        if (Double.isFinite(rate) && rate > 0 && bmax > 0) {
                            UN.set(ist, k, TN.get(ist, k) / rate / bmax);
                        } else {
                            UN.set(ist, k, 0.0);
                        }
                    }
                }
            }
        }

        for (int k = 0; k < K; k++) {
            for (int i = 0; i < M; i++) {
                if (TN.get(i, k) > 0) {
                    RN.set(i, k, QN.get(i, k) / TN.get(i, k));
                } else {
                    RN.set(i, k, 0);
                }
            }
            CN.set(k, NK.get(k) / XN.get(k));
        }
        QN.setNaNToZero();
        CN.setNaNToZero();
        RN.setNaNToZero();
        UN.setNaNToZero();
        XN.setNaNToZero();
        TN.setNaNToZero();

        Result res = new Result();
        res.QN = QN; res.UN = UN; res.RN = RN; res.TN = TN; res.CN = CN; res.XN = XN;
        return res;
    }
}
