/*
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
package jline.inference.api;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.events.EventHandler;
import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator;

import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Station;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Infer_fluid_ps_rt_likelihood {
    private Infer_fluid_ps_rt_likelihood() {}

    /**
     * Result of building the augmented fluid model for response time likelihood.
     */
    public static final class FluidPsRtResult {
        public final FirstOrderDifferentialEquations ode;
        public final Matrix qIndices;
        public final Matrix augPhases;
        public final int stateSize;
        public final int refIdx;
        public final int newK;

        public FluidPsRtResult(FirstOrderDifferentialEquations ode,
                                Matrix qIndices,
                                Matrix augPhases,
                                int stateSize,
                                int refIdx,
                                int newK) {
            this.ode = ode;
            this.qIndices = qIndices;
            this.augPhases = augPhases;
            this.stateSize = stateSize;
            this.refIdx = refIdx;
            this.newK = newK;
        }
    }

    public static FluidPsRtResult infer_fluid_ps_rt_likelihood(NetworkStruct sn, int taggedClass) {
        final int M = sn.nstations;
        final int K = sn.nclasses;
        int N = sn.nclosedjobs;
        final int Kc = K + 1;

        // Find reference station (PS queue) and get server counts
        int refIdxLocal = -1;
        final int[] S = new int[M];
        for (int i = 0; i < M; i++) {
            Station station = sn.stations.get(i);
            if (sn.sched.get(station) != SchedStrategy.INF) {
                refIdxLocal = i;
            }
            S[i] = (int) sn.nservers.get(i, 0);
            if (S[i] == Integer.MAX_VALUE || S[i] < 0) {
                S[i] = N;
            }
        }
        final int refIdx = refIdxLocal;

        final double[][][] new_mu = new double[M][Kc][];
        final double[][][] new_pi = new double[M][Kc][];
        for (int j = 0; j < M; j++) {
            Station station = sn.stations.get(j);
            for (int k = 0; k < K; k++) {
                JobClass jobclass = sn.jobclasses.get(k);
                Matrix muMatrix = sn.mu.get(station) != null ? sn.mu.get(station).get(jobclass) : null;
                Matrix piMatrix = sn.phi.get(station) != null ? sn.phi.get(station).get(jobclass) : null;
                if (muMatrix != null && muMatrix.getNumRows() > 0) {
                    double[] muArr = new double[muMatrix.getNumRows()];
                    for (int r = 0; r < muMatrix.getNumRows(); r++) {
                        muArr[r] = muMatrix.get(r, 0);
                    }
                    double[] piArr;
                    if (piMatrix != null && piMatrix.getNumRows() > 0) {
                        piArr = new double[piMatrix.getNumRows()];
                        for (int r = 0; r < piMatrix.getNumRows(); r++) {
                            piArr[r] = piMatrix.get(r, 0);
                        }
                    } else {
                        piArr = new double[muMatrix.getNumRows()];
                        for (int r = 0; r < piArr.length; r++) piArr[r] = 1.0;
                    }
                    if (muArr.length == 1 && Double.isNaN(muArr[0])) {
                        new_mu[j][k] = null;
                        new_pi[j][k] = null;
                    } else {
                        new_mu[j][k] = muArr;
                        new_pi[j][k] = piArr;
                    }
                }
            }
            // Tagged class Kc-1: copy from taggedClass
            JobClass tcJobclass = sn.jobclasses.get(taggedClass);
            Matrix tcMu = sn.mu.get(station) != null ? sn.mu.get(station).get(tcJobclass) : null;
            Matrix tcPi = sn.phi.get(station) != null ? sn.phi.get(station).get(tcJobclass) : null;
            if (tcMu != null && tcMu.getNumRows() > 0) {
                double[] muArr = new double[tcMu.getNumRows()];
                for (int r = 0; r < tcMu.getNumRows(); r++) {
                    muArr[r] = tcMu.get(r, 0);
                }
                double[] piArr;
                if (tcPi != null && tcPi.getNumRows() > 0) {
                    piArr = new double[tcPi.getNumRows()];
                    for (int r = 0; r < tcPi.getNumRows(); r++) {
                        piArr[r] = tcPi.get(r, 0);
                    }
                } else {
                    piArr = new double[tcMu.getNumRows()];
                    for (int r = 0; r < piArr.length; r++) piArr[r] = 1.0;
                }
                if (muArr.length == 1 && Double.isNaN(muArr[0])) {
                    new_mu[j][Kc - 1] = null;
                    new_pi[j][Kc - 1] = null;
                } else {
                    new_mu[j][Kc - 1] = muArr;
                    new_pi[j][Kc - 1] = piArr;
                }
            }
        }

        final Matrix augPhases = new Matrix(M, Kc);
        for (int i = 0; i < M; i++) {
            for (int c = 0; c < Kc; c++) {
                double[] mu = new_mu[i][c];
                if (mu != null) {
                    augPhases.set(i, c, (double) mu.length);
                }
            }
        }

        final Matrix qIndices = new Matrix(M, Kc);
        int idx = 0;
        for (int i = 0; i < M; i++) {
            for (int c = 0; c < Kc; c++) {
                qIndices.set(i, c, (double) idx);
                int phases = (int) augPhases.get(i, c);
                if (phases > 0) {
                    idx += phases;
                }
            }
        }
        final int totalPhases = idx;

        // Build expanded routing table
        Matrix rt = sn.rt;
        final double[][] newRt = new double[M * Kc][M * Kc];

        for (int l = 0; l < K; l++) {
            for (int m = 0; m < K; m++) {
                for (int i = 0; i < M; i++) {
                    for (int j = 0; j < M; j++) {
                        newRt[i * Kc + l][j * Kc + m] = rt.get(i * K + l, j * K + m);
                    }
                }
            }
        }

        for (int i = 0; i < M; i++) {
            for (int j = 0; j < M; j++) {
                newRt[i * Kc + (Kc - 1)][j * Kc + (Kc - 1)] = rt.get(i * K + taggedClass, j * K + taggedClass);
            }
        }

        int chainIdx = -1;
        if (sn.chains != null) {
            for (int ch = 0; ch < sn.chains.getNumRows(); ch++) {
                if (sn.chains.get(ch, taggedClass) == 1.0) {
                    chainIdx = ch;
                    break;
                }
            }
        }
        if (chainIdx >= 0) {
            List<Integer> classesInChain = new ArrayList<Integer>();
            for (int c = 0; c < K; c++) {
                if (sn.chains.get(chainIdx, c) == 1.0) {
                    classesInChain.add(Integer.valueOf(c));
                }
            }
            for (Integer lBox : classesInChain) {
                int l = lBox.intValue();
                for (int j = 0; j < M; j++) {
                    newRt[refIdx * Kc + (Kc - 1)][j * Kc + l] = rt.get(refIdx * K + taggedClass, j * K + l);
                }
            }
        }
        for (int j = 0; j < M; j++) {
            newRt[refIdx * Kc + (Kc - 1)][j * Kc + (Kc - 1)] = 0.0;
        }

        final int stateSize = totalPhases;
        final NetworkStruct snFinal = sn;
        FirstOrderDifferentialEquations ode = new FirstOrderDifferentialEquations() {
            @Override
            public int getDimension() {
                return stateSize;
            }

            @Override
            public void computeDerivatives(double t, double[] y, double[] yDot) {
                for (int i = 0; i < yDot.length; i++) yDot[i] = 0.0;

                for (int i = 0; i < M; i++) {
                    Station station = snFinal.stations.get(i);
                    boolean isPS = snFinal.sched.get(station) != SchedStrategy.INF;

                    double nAtI = 0.0;
                    for (int c = 0; c < Kc; c++) {
                        if ((int) augPhases.get(i, c) > 0) {
                            nAtI += Math.max(0.0, y[(int) qIndices.get(i, c)]);
                        }
                    }

                    double capacity = (double) S[i];

                    for (int c = 0; c < Kc; c++) {
                        int phases = (int) augPhases.get(i, c);
                        if (phases > 0) {
                            int qi = (int) qIndices.get(i, c);
                            double muVal = (new_mu[i][c] != null) ? new_mu[i][c][0] : 0.0;
                            double phiVal = (new_pi[i][c] != null) ? new_pi[i][c][0] : 1.0;

                            if (nAtI > 0 && muVal > 0) {
                                double share = y[qi] / Math.max(nAtI, 1e-10);
                                double effectiveRate;
                                if (isPS) {
                                    effectiveRate = muVal * Math.min(nAtI, capacity) * share;
                                } else {
                                    effectiveRate = muVal * Math.max(0.0, y[qi]);
                                }

                                yDot[qi] -= effectiveRate * phiVal;

                                for (int j = 0; j < M; j++) {
                                    for (int d = 0; d < Kc; d++) {
                                        double pcd = newRt[i * Kc + c][j * Kc + d];
                                        if (pcd > 0 && (int) augPhases.get(j, d) > 0) {
                                            int qj = (int) qIndices.get(j, d);
                                            yDot[qj] += effectiveRate * phiVal * pcd;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        };

        return new FluidPsRtResult(ode, qIndices, augPhases, stateSize, refIdx, Kc);
    }

    public static double infer_fluid_ps_rt_solve(
            FluidPsRtResult result,
            Matrix y0Levels,
            double Rsampled,
            int taggedClass) {
        final double newFluid = 1.0;
        int M = y0Levels.getNumRows();
        int K = result.newK - 1;

        double[] y0 = new double[result.stateSize];
        for (int i = 0; i < M; i++) {
            for (int k = 0; k < K; k++) {
                if ((int) result.augPhases.get(i, k) > 0) {
                    int idx = (int) result.qIndices.get(i, k);
                    y0[idx] = y0Levels.get(i, k);
                }
            }
        }

        final int refTaggedOrigIdx = (int) result.qIndices.get(result.refIdx, taggedClass);
        final int refTaggedNewIdx = (int) result.qIndices.get(result.refIdx, K);
        y0[refTaggedOrigIdx] -= newFluid;
        y0[refTaggedNewIdx] = newFluid;

        for (int i = 0; i < y0.length; i++) y0[i] = Math.max(0.0, y0[i]);

        try {
            DormandPrince853Integrator integrator = new DormandPrince853Integrator(1e-10, Rsampled, 1e-8, 1e-5);

            integrator.addEventHandler(new EventHandler() {
                @Override
                public void init(double t0, double[] y0, double t) {}

                @Override
                public double g(double t, double[] y) {
                    return y[refTaggedNewIdx];
                }

                @Override
                public Action eventOccurred(double t, double[] y, boolean increasing) {
                    return Action.STOP;
                }

                @Override
                public void resetState(double t, double[] y) {}
            }, 1e-6, 1e-8, 100);

            double[] yFinal = y0.clone();
            double tFinal = integrator.integrate(result.ode, 0.0, y0, Rsampled, yFinal);

            if (Rsampled <= tFinal + 1e-10) {
                double[] lastRates = new double[result.stateSize];
                result.ode.computeDerivatives(tFinal, yFinal, lastRates);
                return Math.max(0.0, -lastRates[refTaggedNewIdx] / newFluid);
            }
        } catch (Exception e) {
            // ODE solver failed
        }
        return 0.0;
    }
}
