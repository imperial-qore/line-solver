/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */


package jline.solvers.fluid.handlers;

import jline.GlobalConstants;
import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.nodes.Station;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

import java.util.Map;
import java.util.Objects;

import static jline.api.mam.Map_pieKt.map_pie;
import static jline.io.InputOutputKt.line_error;
import static jline.io.InputOutputKt.mfilename;
import static jline.lang.constant.SchedStrategy.DPS;
import static jline.util.Maths.softmin;
import static org.apache.commons.math3.util.FastMath.min;

public class PassageTimeODE implements FirstOrderDifferentialEquations {
    // Softmin method implemented - provides smooth approximation of min function for ODE stability
    private final NetworkStruct sn;
    private final Map<Station, Map<JobClass, Matrix>> mu;
    private final Map<Station, Map<JobClass, Matrix>> phi;
    private final Map<Station, Map<JobClass, MatrixCell>> proc;
    private final Matrix rt;
    private final Matrix nservers;
    private final SolverOptions options;
    private final int numDimensions;

    // Precomputed structures (computed once in constructor, reused in computeDerivatives)
    private final boolean[][] cachedEnabled;
    private final Matrix cachedQIndices;
    private final Matrix cachedKic;
    private final Matrix cachedW;
    // Precomputed for the default (closing) method only
    private final Matrix cachedAllJumps;
    private final Matrix cachedRateBase;
    private final Matrix cachedEventIdx;

    public PassageTimeODE(
            NetworkStruct sn,
            Map<Station, Map<JobClass, Matrix>> mu,
            Map<Station, Map<JobClass, Matrix>> phi,
            Map<Station, Map<JobClass, MatrixCell>> proc,
            Matrix rt,
            Matrix S,
            SolverOptions options,
            int numDimensions) {
        this.sn = sn;
        this.mu = mu;
        this.phi = phi;
        this.proc = proc;
        this.rt = rt;
        this.nservers = S;
        this.options = options;
        this.numDimensions = numDimensions;

        // Precompute enabled, qIndices, Kic, w (these don't change between ODE steps)
        int M = S.length();
        int K = mu.get(sn.stations.get(0)).size();

        this.cachedEnabled = new boolean[M][K];
        this.cachedQIndices = new Matrix(M, K);
        this.cachedKic = new Matrix(M, K);
        this.cachedW = new Matrix(M, K);
        int cumSum = 0;

        for (int i = 0; i < M; i++) {
            Station station = sn.stations.get(i);
            for (int c = 0; c < K; c++) {
                JobClass jobClass = sn.jobclasses.get(c);
                cachedW.set(i, c, 1);
                cachedEnabled[i][c] = false;
                int numPhases = 0;

                Matrix muMatrix = mu.get(station).get(jobClass);
                if (muMatrix != null) {
                    int numNans = 0;
                    for (int row = 0; row < muMatrix.getNumRows(); row++) {
                        for (int col = 0; col < muMatrix.getNumCols(); col++) {
                            if (Double.isNaN(muMatrix.get(row, col))) {
                                numNans++;
                            }
                        }
                    }
                    if ((numNans != muMatrix.getNumElements()) && !muMatrix.isEmpty()) {
                        numPhases = muMatrix.length();
                        cachedEnabled[i][c] = true;
                    }
                }

                cachedQIndices.set(i, c, cumSum);
                cachedKic.set(i, c, numPhases);
                cumSum += numPhases;
            }

            if (sn.sched.get(station) == DPS) {
                for (int k = 0; k < K; k++) {
                    cachedW.set(i, k, sn.schedparam.get(i, k));
                }
                double sumWI = cachedW.sumRows(i);
                if (sumWI > 0) {
                    for (int k = 0; k < K; k++) {
                        cachedW.set(i, k, cachedW.get(i, k) / sumWI);
                    }
                }
            }
        }

        // Precompute allJumps, rateBase, eventIdx for the default (closing) method
        if (!Objects.equals(options.method, "statedep") && !Objects.equals(options.method, "softmin")) {
            Matrix tmpAllJumps = calculateJumps(cachedEnabled, cachedQIndices, cachedKic);
            Matrix tmpRateBase = new Matrix(tmpAllJumps.getNumCols(), 1);
            Matrix tmpEventIdx = new Matrix(tmpAllJumps.getNumCols(), 1);
            calculateRateBaseAndEventIdxs(cachedEnabled, cachedQIndices, cachedKic, tmpRateBase, tmpEventIdx);

            // Apply immediate elimination if configured (matching MATLAB solver_fluid_odes.m)
            if (options.config != null && options.config.hide_immediate) {
                ImmediateElimination.EliminationResult result =
                        ImmediateElimination.eliminateImmediate(tmpAllJumps, tmpRateBase, tmpEventIdx, sn, options);
                this.cachedAllJumps = result.allJumpsReduced;
                this.cachedRateBase = result.rateBaseReduced;
                this.cachedEventIdx = result.eventIdxReduced;
            } else {
                this.cachedAllJumps = tmpAllJumps;
                this.cachedRateBase = tmpRateBase;
                this.cachedEventIdx = tmpEventIdx;
            }
        } else {
            this.cachedAllJumps = null;
            this.cachedRateBase = null;
            this.cachedEventIdx = null;
        }
    }

    public PassageTimeODE(
            NetworkStruct sn,
            Map<Station, Map<JobClass, Matrix>> mu,
            Map<Station, Map<JobClass, Matrix>> phi,
            Map<Station, Map<JobClass, MatrixCell>> proc,
            Matrix rt,
            Matrix S,
            SolverOptions options) {
        this(sn, mu, phi, proc, rt, S, options, options.init_sol.length());
    }

    private Matrix calculateJumps(boolean[][] enabled, Matrix qIndices, Matrix Kic) {

        int M = sn.nstations; // Number of stations
        int K = mu.get(sn.stations.get(0)).size(); // Number of classes
        int jumpsRows = (int) Kic.elementSum();
        Matrix jumps =
                new Matrix(jumpsRows, 0); // Returns state changes triggered by all the events

        for (int i = 0; i < M; i++) { // state changes from departures in service phases 2
            for (int c = 0; c < K; c++) {
                if (enabled[i][c]) {
                    int xic = (int) qIndices.get(i, c); //  index of x_ic
                    for (int j = 0; j < M; j++) {
                        for (int l = 0; l < K; l++) {
                            if (rt.get(i * K + c, j * K + l) > 0) {
                                int xjl = (int) qIndices.get(j, l); // index of x_jl
                                for (int ki = 0; ki < Kic.get(i, c); ki++) { // job can leave from any phase in i
                                    for (int kj = 0; kj < Kic.get(j, l); kj++) { // job can start from any phase in j
                                        setNextJump(jumps, xic + ki, xjl + kj);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        for (int i = 0; i < M; i++) { // state changes: "next service phase" transition
            for (int c = 0; c < K; c++) {
                if (enabled[i][c]) {
                    int xic = (int) qIndices.get(i, c);
                    for (int ki = 0; ki < Kic.get(i, c) - 1; ki++) {
                        for (int kip = 0; kip < Kic.get(i, c); kip++) {
                            if (ki != kip) {
                                setNextJump(jumps, xic + ki, xic + kip);
                            }
                        }
                    }
                }
            }
        }

        return jumps;
    }

    private void calculateRateBaseAndEventIdxs(
            boolean[][] enabled,
            Matrix qIndices,
            Matrix Kic,
            Matrix rateBase,
            Matrix eventIdx) {

        int M = sn.nstations; // Number of stations
        int K = mu.get(sn.stations.get(0)).size(); // Number of classes
        int rateIdx = 0;

        // State changes from departures in service phases 2...
        for (int i = 0; i < M; i++) {
            for (int c = 0; c < K; c++) {
                if (enabled[i][c]) {
                    for (int j = 0; j < M; j++) {
                        for (int l = 0; l < K; l++) {
                            Matrix pie;
                            if (proc.get(sn.stations.get(j)).get(sn.jobclasses.get(l)).isEmpty()) {
                                pie = new Matrix(1, 1, 1);
                                pie.set(0, 0, 1);
                            } else {
                                Matrix D0 = proc.get(sn.stations.get(j)).get(sn.jobclasses.get(l)).get(0);
                                Matrix D1 = proc.get(sn.stations.get(j)).get(sn.jobclasses.get(l)).get(1);

                                pie = map_pie(D0, D1);
                            }
                            if (rt.get(i * K + c, j * K + l) > 0) {
                                for (int kicIdx = 0; kicIdx < Kic.get(i, c); kicIdx++) {
                                    for (int kjl = 0; kjl < Kic.get(j, l); kjl++) {
                                        rateBase.set(
                                                rateIdx,
                                                0,
                                                phi.get(sn.stations.get(i)).get(sn.jobclasses.get(c)).get(kicIdx, 0)
                                                        * mu.get(sn.stations.get(i)).get(sn.jobclasses.get(c)).get(kicIdx, 0)
                                                        * rt.get(i * K + c, j * K + l)
                                                        * pie.get(0, kjl));
                                        eventIdx.set(rateIdx, 0, qIndices.get(i, c) + kicIdx);
                                        rateIdx++;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        // State changes from "next service phase" transition in phases 2...
        for (int i = 0; i < M; i++) {
            for (int c = 0; c < K; c++) {
                if (enabled[i][c]) {
                    for (int kicIdx = 0; kicIdx < Kic.get(i, c) - 1; kicIdx++) {
                        for (int kicp = 0; kicp < Kic.get(i, c); kicp++) {
                            if (kicp != kicIdx) {
                                rateBase.set(
                                        rateIdx,
                                        0,
                                        proc.get(sn.stations.get(i))
                                                .get(sn.jobclasses.get(c))
                                                .get(0)
                                                .get(kicIdx, kicp));
                                eventIdx.set(rateIdx, 0, qIndices.get(i, c) + kicIdx);
                                rateIdx++;
                            }
                        }
                    }
                }
            }
        }
    }

    private Matrix calculatedxdtClosingMethod(
            double[] x,
            Matrix w,
            boolean[][] enabled,
            Matrix qIndices,
            Matrix Kic,
            Matrix allJumps,
            Matrix rateBase,
            Matrix eventIdx) {

        int M = sn.nstations; // Number of stations
        int K = mu.get(sn.stations.get(0)).size(); // Number of classes

        // Basic vector valid for INF and PS case min(ni, nservers(i)) = ni
        Matrix rates = new Matrix(x.length, 1);
        for (int i = 0; i < x.length; i++) {
            rates.set(i, 0, x[i]);
        }

        // Declare variables outside switch to avoid scope issues
        int idxIni, idxEnd;
        double ni;
        
        for (int i = 0; i < M; i++) {
            switch (sn.sched.get(sn.stations.get(i))) {
                case INF:
                    break;
                case EXT:
                    // This is treated by a delay except that we require mass conservation in the local
                    // population
                    for (int k = 0; k < K; k++) {
                        idxIni = (int) qIndices.get(i, k);
                        idxEnd = (int) qIndices.get(i, k) + (int) Kic.get(i, k);
                        if (enabled[i][k]) {
                            // Keep total mass 1 into the source for all classes at all
                            // times, not needed for idxIni+1:idxEnd as rates is initialized equal to x
                            double tmpSum = 0;
                            for (int idx = idxIni + 1; idx < idxEnd; idx++) {
                                tmpSum += x[idx];
                            }
                            rates.set(idxIni, 0, 1 - tmpSum);
                        }
                    }
                    break;
                case PS:
                case FCFS:
                    idxIni = (int) qIndices.get(i, 0);
                    idxEnd = (int) qIndices.get(i, K - 1) + (int) Kic.get(i, K - 1);
                    ni = 0;
                    for (int idx = idxIni; idx < idxEnd; idx++) {
                        ni += x[idx];
                    }
                    if (ni > nservers.get(i, 0)) { // case min = ni handled by rates = x
                        for (int idx = idxIni; idx < idxEnd; idx++) {
                            rates.set(idx, 0, x[idx] / ni * nservers.get(i, 0));
                        }
                    }
                    break;
                case DPS:
                    double sumWI = w.sumRows(i);
                    for (int col = 0; col < K; col++) {
                        w.set(i, col, w.get(i, col) / sumWI);
                    }
                    sumWI = w.sumRows(i);
                    ni = sumWI / K;

                    for (int k = 0; k < K; k++) {
                        idxIni = (int) qIndices.get(i, k);
                        idxEnd = (int) qIndices.get(i, k) + (int) Kic.get(i, k);
                        if (enabled[i][k]) {
                            double tmpSum = 0;
                            for (int idx = idxIni; idx < idxEnd; idx++) {
                                tmpSum += x[idx] * w.get(i, k);
                            }
                            ni += tmpSum;
                        }
                    }

                    for (int k = 0; k < K; k++) {
                        idxIni = (int) qIndices.get(i, k);
                        idxEnd = (int) qIndices.get(i, k) + (int) Kic.get(i, k);
                        if (enabled[i][k]) {
                            for (int idx = idxIni; idx < idxEnd; idx++) {
                                // Not needed for idxIni+1:idxEnd as rates is initialised equal to x
                                rates.set(idx, 0, w.get(i, k) * x[idx] / ni * nservers.get(i, 0));
                            }
                        }
                    }
            }
        }

        int numEventIndices = eventIdx.getNumRows(); // Use getNumRows(), not length() which returns max(rows,cols)
        Matrix newRates = new Matrix(numEventIndices, 1);
        for (int i = 0; i < numEventIndices; i++) {
            newRates.set(i, 0, rates.get((int) eventIdx.get(i, 0), 0));
        }
        newRates.elementMult(rateBase, newRates);
        
        
        return allJumps.mult(newRates, null);
    }

    private Matrix calculatedxdtStateDepMethod(
            double[] x, boolean[][] enabled, Matrix qIndices, Matrix Kic, Matrix w) {

        int M = sn.nstations; // Number of stations
        int K = mu.get(sn.stations.get(0)).size(); // Number of classes
        Matrix dxdt = new Matrix(x.length, 1);

        // Declare variables outside switch to avoid scope issues
        int idxIni, idxEnd;
        double ni;
        
        for (int i = 0; i < M; i++) {
            switch (sn.sched.get(sn.stations.get(i))) {
                case INF:
                    // Phase changes
                    for (int c = 0; c < K; c++) {
                        if (enabled[i][c]) {
                            int xic = (int) qIndices.get(i, c);
                            for (int kicIdx = 0; kicIdx < Kic.get(i, c) - 1; kicIdx++) {
                                for (int kic_p = 0; kic_p < Kic.get(i, c); kic_p++) {
                                    if (kicIdx != kic_p) {
                                        double rate =
                                                proc.get(sn.stations.get(i))
                                                        .get(sn.jobclasses.get(c))
                                                        .get(0)
                                                        .get(kicIdx, kic_p);
                                        dxdt.set(xic + kicIdx, 0, dxdt.get(xic + kicIdx, 0) - (x[xic + kicIdx] * rate));
                                        dxdt.set(xic + kic_p, 0, dxdt.get(xic + kic_p, 0) + (x[xic + kicIdx] * rate));
                                    }
                                }
                            }
                        }
                    }
                    // Service completions
                    for (int c = 0; c < K; c++) {
                        if (enabled[i][c]) {
                            int xic = (int) qIndices.get(i, c);
                            for (int j = 0; j < M; j++) {
                                for (int l = 0; l < K; l++) {
                                    int xjl = (int) qIndices.get(j, l);
                                    if (enabled[j][l]) {
                                        Matrix D0 = proc.get(sn.stations.get(j)).get(sn.jobclasses.get(l)).get(0);
                                        Matrix D1 = proc.get(sn.stations.get(j)).get(sn.jobclasses.get(l)).get(1);
                                        Matrix pie =
                                                map_pie(D0, D1);
                                        if (rt.get(i * K + c, j * K + l) > 0) {
                                            for (int kicIdx = 0; kicIdx < Kic.get(i, c); kicIdx++) {
                                                for (int kjl = 0; kjl < Kic.get(j, l); kjl++) {
                                                    if (j != i) {
                                                        double rate =
                                                                phi.get(sn.stations.get(i)).get(sn.jobclasses.get(c)).get(kicIdx, 0)
                                                                        * mu.get(sn.stations.get(i))
                                                                        .get(sn.jobclasses.get(c))
                                                                        .get(kicIdx, 0)
                                                                        * rt.get(i * K + c, j * K + l)
                                                                        * pie.get(0, kjl);
                                                        dxdt.set(
                                                                xic + kicIdx,
                                                                0,
                                                                dxdt.get(xic + kicIdx, 0) - (x[xic + kicIdx] * rate));
                                                        dxdt.set(
                                                                xjl + kjl, 0, dxdt.get(xjl + kjl, 0) + (x[xic + kicIdx] * rate));
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    break;

                case EXT:
                    // TODO: open models for state-dep
                    line_error(mfilename(new Object[]{}),
                            "State dependent ODE method does not support open models. Try with default method.");
                    break;

                case PS:
                    idxIni = (int) qIndices.get(i, 0);
                    idxEnd = (int) qIndices.get(i, K - 1) + (int) Kic.get(i, K - 1);
                    ni = 0;
                    for (int idx = idxIni; idx < idxEnd; idx++) {
                        ni += x[idx];
                    }
                    // Phase changes
                    for (int c = 0; c < K; c++) {
                        if (enabled[i][c]) {
                            int xic = (int) qIndices.get(i, c);
                            for (int kicIdx = 0; kicIdx < Kic.get(i, c) - 1; kicIdx++) {
                                for (int kic_p = 0; kic_p < Kic.get(i, c); kic_p++) {
                                    if (kicIdx != kic_p) {
                                        double rate =
                                                proc.get(sn.stations.get(i))
                                                        .get(sn.jobclasses.get(c))
                                                        .get(0)
                                                        .get(kicIdx, kic_p);
                                        if (ni > sn.nservers.get(i, 0)) {
                                            dxdt.set(
                                                    xic + kicIdx,
                                                    0,
                                                    dxdt.get(xic + kicIdx, 0)
                                                            - (x[xic + kicIdx] * rate * sn.nservers.get(i, 0) / ni));
                                            dxdt.set(
                                                    xic + kic_p,
                                                    0,
                                                    dxdt.get(xic + kic_p, 0)
                                                            + (x[xic + kicIdx] * rate * sn.nservers.get(i, 0) / ni));
                                        } else {
                                            dxdt.set(
                                                    xic + kicIdx, 0, dxdt.get(xic + kicIdx, 0) - (x[xic + kicIdx] * rate));
                                            dxdt.set(xic + kic_p, 0, dxdt.get(xic + kic_p, 0) + (x[xic + kicIdx] * rate));
                                        }
                                    }
                                }
                            }
                        }
                    }
                    // Service completions
                    for (int c = 0; c < K; c++) {
                        if (enabled[i][c]) {
                            int xic = (int) qIndices.get(i, c);
                            for (int j = 0; j < M; j++) {
                                for (int l = 0; l < K; l++) {
                                    int xjl = (int) qIndices.get(j, l);
                                    if (enabled[j][l]) {
                                        Matrix D0 = proc.get(sn.stations.get(j)).get(sn.jobclasses.get(l)).get(0);
                                        Matrix D1 = proc.get(sn.stations.get(j)).get(sn.jobclasses.get(l)).get(1);
                                        Matrix pie = map_pie(D0, D1);

                                        if (rt.get(i * K + c, j * K + l) > 0) {
                                            for (int kicIdx = 0; kicIdx < Kic.get(i, c); kicIdx++) {
                                                for (int kjl = 0; kjl < Kic.get(j, l); kjl++) {
                                                    double rate =
                                                            phi.get(sn.stations.get(i)).get(sn.jobclasses.get(c)).get(kicIdx, 0)
                                                                    * mu.get(sn.stations.get(i))
                                                                    .get(sn.jobclasses.get(c))
                                                                    .get(kicIdx, 0)
                                                                    * rt.get(i * K + c, j * K + l)
                                                                    * pie.get(0, kjl);
                                                    if (ni > sn.nservers.get(i, 0)) {
                                                        rate = 1/ni * sn.nservers.get(i, 0) * rate;
                                                    }
                                                    dxdt.set(
                                                            xic + kicIdx,
                                                            0,
                                                            dxdt.get(xic + kicIdx, 0) - (x[xic + kicIdx] * rate));
                                                    dxdt.set(xjl + kjl, 0, dxdt.get(xjl + kjl, 0) + (x[xic + kicIdx] * rate));
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    break;

                case FCFS:
                    idxIni = (int) qIndices.get(i, 0);
                    idxEnd = (int) qIndices.get(i, K - 1) + (int) Kic.get(i, K - 1);
                    ni = 0;
                    for (int idx = idxIni; idx < idxEnd; idx++) {
                        ni += x[idx];
                    }
                    double wni = GlobalConstants.FineTol;
                    for (int c = 0; c < K; c++) {
                        for (int kicIdx = 0; kicIdx < Kic.get(i, c); kicIdx++) {
                            if (enabled[i][c]) {
                                int xic = (int) qIndices.get(i, c);
                                w.set(
                                        c,
                                        kicIdx,
                                        -1
                                                / proc.get(sn.stations.get(i))
                                                .get(sn.jobclasses.get(c))
                                                .get(0)
                                                .get(kicIdx, kicIdx));
                                wni += w.get(c, kicIdx) * x[xic + kicIdx];
                            }
                        }
                    }
                    // Phase changes
                    for (int c = 0; c < K; c++) {
                        if (enabled[i][c]) {
                            int xic = (int) qIndices.get(i, c);
                            for (int kicIdx = 0; kicIdx < Kic.get(i, c) - 1; kicIdx++) {
                                for (int kic_p = 0; kic_p < Kic.get(i, c); kic_p++) {
                                    if (kicIdx != kic_p) {
                                        double rate =
                                                proc.get(sn.stations.get(i))
                                                        .get(sn.jobclasses.get(c))
                                                        .get(0)
                                                        .get(kicIdx, kic_p)
                                                        * min(ni, sn.nservers.get(i, 0))
                                                        * w.get(c, kicIdx)
                                                        / wni;
                                        dxdt.set(xic + kicIdx, 0, dxdt.get(xic + kicIdx, 0) - (x[xic + kicIdx] * rate));
                                        dxdt.set(xic + kic_p, 0, dxdt.get(xic + kic_p, 0) + (x[xic + kicIdx] * rate));
                                    }
                                }
                            }
                        }
                    }
                    // Service completions
                    for (int c = 0; c < K; c++) {
                        if (enabled[i][c]) {
                            int xic = (int) qIndices.get(i, c);
                            for (int j = 0; j < M; j++) {
                                for (int l = 0; l < K; l++) {
                                    int xjl = (int) qIndices.get(j, l);
                                    if (enabled[j][l]) {
                                        Matrix D0 = proc.get(sn.stations.get(j)).get(sn.jobclasses.get(l)).get(0);
                                        Matrix D1 = proc.get(sn.stations.get(j)).get(sn.jobclasses.get(l)).get(1);
                                        Matrix pie = map_pie(D0, D1);

                                        if (rt.get(i * K + c, j * K + l) > 0) {
                                            for (int kicIdx = 0; kicIdx < Kic.get(i, c); kicIdx++) {
                                                for (int kjl = 0; kjl < Kic.get(j, l); kjl++) {
                                                    double rate =
                                                            phi.get(sn.stations.get(i)).get(sn.jobclasses.get(c)).get(kicIdx, 0)
                                                                    * mu.get(sn.stations.get(i))
                                                                    .get(sn.jobclasses.get(c))
                                                                    .get(kicIdx, 0)
                                                                    * rt.get(i * K + c, j * K + l)
                                                                    * pie.get(0, kjl)
                                                                    * min(ni, sn.nservers.get(i, 0))
                                                                    * w.get(c, kicIdx)
                                                                    / wni;
                                                    dxdt.set(
                                                            xic + kicIdx,
                                                            0,
                                                            dxdt.get(xic + kicIdx, 0) - (x[xic + kicIdx] * rate));
                                                    dxdt.set(xjl + kjl, 0, dxdt.get(xjl + kjl, 0) + (x[xic + kicIdx] * rate));
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    break;

                case DPS:
                    double sumWI = w.sumRows(i);
                    for (int col = 0; col < w.getNumCols(); col++) {
                        w.set(i, col, w.get(i, col) / sumWI);
                    }
                    idxIni = (int) qIndices.get(i, 0);
                    idxEnd = (int) qIndices.get(i, K - 1) + (int) Kic.get(i, K - 1);
                    wni = 0;
                    for (int idx = idxIni; idx < idxEnd; idx++) {
                        wni += x[idx];
                    }
                    // Phase changes
                    for (int c = 0; c < K; c++) {
                        if (enabled[i][c]) {
                            int xic = (int) qIndices.get(i, c);
                            for (int kicIdx = 0; kicIdx < Kic.get(i, c) - 1; kicIdx++) {
                                for (int kic_p = 0; kic_p < Kic.get(i, c); kic_p++) {
                                    if (kicIdx != kic_p) {
                                        double rate =
                                                proc.get(sn.stations.get(i))
                                                        .get(sn.jobclasses.get(c))
                                                        .get(0)
                                                        .get(kicIdx, kic_p);
                                        if (wni > sn.nservers.get(i, 0)) {
                                            dxdt.set(
                                                    xic + kicIdx,
                                                    0,
                                                    dxdt.get(xic + kicIdx, 0)
                                                            - (x[xic + kicIdx]
                                                            * rate
                                                            * sn.nservers.get(i, 0)
                                                            * w.get(c, kicIdx)
                                                            / wni));
                                            dxdt.set(
                                                    xic + kic_p,
                                                    0,
                                                    dxdt.get(xic + kic_p, 0)
                                                            + (x[xic + kicIdx]
                                                            * rate
                                                            * sn.nservers.get(i, 0)
                                                            * w.get(c, kicIdx)
                                                            / wni));
                                        } else {
                                            dxdt.set(
                                                    xic + kicIdx, 0, dxdt.get(xic + kicIdx, 0) - (x[xic + kicIdx] * rate));
                                            dxdt.set(xic + kic_p, 0, dxdt.get(xic + kic_p, 0) + (x[xic + kicIdx] * rate));
                                        }
                                    }
                                }
                            }
                        }
                    }
                    // Service completions
                    for (int c = 0; c < K; c++) {
                        if (enabled[i][c]) {
                            int xic = (int) qIndices.get(i, c);
                            for (int j = 0; j < M; j++) {
                                for (int l = 0; l < K; l++) {
                                    int xjl = (int) qIndices.get(j, l);
                                    if (enabled[j][l]) {
                                        Matrix D0 = proc.get(sn.stations.get(j)).get(sn.jobclasses.get(l)).get(0);
                                        Matrix D1 = proc.get(sn.stations.get(j)).get(sn.jobclasses.get(l)).get(1);
                                        Matrix pie = map_pie(D0, D1);

                                        if (rt.get(i * K + c, j * K + l) > 0) {
                                            for (int kicIdx = 0; kicIdx < Kic.get(i, c); kicIdx++) {
                                                for (int kjl = 0; kjl < Kic.get(j, l); kjl++) {
                                                    double rate =
                                                            phi.get(sn.stations.get(i)).get(sn.jobclasses.get(c)).get(kicIdx, 0)
                                                                    * mu.get(sn.stations.get(i))
                                                                    .get(sn.jobclasses.get(c))
                                                                    .get(kicIdx, 0)
                                                                    * rt.get(i * K + c, j * K + l)
                                                                    * pie.get(0, kjl);
                                                    if (wni > sn.nservers.get(i, 0)) {
                                                        rate *= sn.nservers.get(i, 0) * w.get(c, kicIdx) / wni;
                                                    }
                                                    dxdt.set(
                                                            xic + kicIdx,
                                                            0,
                                                            dxdt.get(xic + kicIdx, 0) - (x[xic + kicIdx] * rate));
                                                    dxdt.set(xjl + kjl, 0, dxdt.get(xjl + kjl, 0) + (x[xic + kicIdx] * rate));
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
            }
        }

        return dxdt;
    }

    private Matrix calculatedxdtSoftminMethod(
            double[] x, boolean[][] enabled, Matrix qIndices, Matrix Kic, Matrix w) {

        int M = sn.nstations; // Number of stations
        int K = mu.get(sn.stations.get(0)).size(); // Number of classes
        Matrix dxdt = new Matrix(x.length, 1);
        double alpha = 20.0; // Softmin smoothing parameter (as used in MATLAB implementation)
        
        // Declare variables outside switch to avoid scope issues
        int idxIni, idxEnd;
        double ni;

        for (int i = 0; i < M; i++) {
            switch (sn.sched.get(sn.stations.get(i))) {
                case INF:
                    // Phase changes
                    for (int c = 0; c < K; c++) {
                        if (enabled[i][c]) {
                            int xic = (int) qIndices.get(i, c);
                            for (int kicIdx = 0; kicIdx < Kic.get(i, c) - 1; kicIdx++) {
                                for (int kic_p = 0; kic_p < Kic.get(i, c); kic_p++) {
                                    if (kicIdx != kic_p) {
                                        double rate =
                                                proc.get(sn.stations.get(i))
                                                        .get(sn.jobclasses.get(c))
                                                        .get(0)
                                                        .get(kicIdx, kic_p);
                                        dxdt.set(xic + kicIdx, 0, dxdt.get(xic + kicIdx, 0) - (x[xic + kicIdx] * rate));
                                        dxdt.set(xic + kic_p, 0, dxdt.get(xic + kic_p, 0) + (x[xic + kicIdx] * rate));
                                    }
                                }
                            }
                        }
                    }
                    // Service completions
                    for (int c = 0; c < K; c++) {
                        if (enabled[i][c]) {
                            int xic = (int) qIndices.get(i, c);
                            for (int j = 0; j < M; j++) {
                                for (int l = 0; l < K; l++) {
                                    int xjl = (int) qIndices.get(j, l);
                                    if (enabled[j][l]) {
                                        Matrix D0 = proc.get(sn.stations.get(j)).get(sn.jobclasses.get(l)).get(0);
                                        Matrix D1 = proc.get(sn.stations.get(j)).get(sn.jobclasses.get(l)).get(1);
                                        Matrix pie = map_pie(D0, D1);
                                        if (rt.get(i * K + c, j * K + l) > 0) {
                                            for (int kicIdx = 0; kicIdx < Kic.get(i, c); kicIdx++) {
                                                for (int kjl = 0; kjl < Kic.get(j, l); kjl++) {
                                                    if (j != i) {
                                                        double rate =
                                                                phi.get(sn.stations.get(i)).get(sn.jobclasses.get(c)).get(kicIdx, 0)
                                                                        * mu.get(sn.stations.get(i))
                                                                        .get(sn.jobclasses.get(c))
                                                                        .get(kicIdx, 0)
                                                                        * rt.get(i * K + c, j * K + l)
                                                                        * pie.get(0, kjl);
                                                        dxdt.set(
                                                                xic + kicIdx,
                                                                0,
                                                                dxdt.get(xic + kicIdx, 0) - (x[xic + kicIdx] * rate));
                                                        dxdt.set(
                                                                xjl + kjl, 0, dxdt.get(xjl + kjl, 0) + (x[xic + kicIdx] * rate));
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    break;

                case EXT:
                    // TODO: open models for state-dep
                    line_error(mfilename(new Object[]{}),
                            "Softmin ODE method does not support open models. Try with default method.");
                    break;

                case PS:
                    idxIni = (int) qIndices.get(i, 0);
                    idxEnd = (int) qIndices.get(i, K - 1) + (int) Kic.get(i, K - 1);
                    ni = 0;
                    for (int idx = idxIni; idx < idxEnd; idx++) {
                        ni += x[idx];
                    }
                    // Phase changes
                    for (int c = 0; c < K; c++) {
                        if (enabled[i][c]) {
                            int xic = (int) qIndices.get(i, c);
                            for (int kicIdx = 0; kicIdx < Kic.get(i, c) - 1; kicIdx++) {
                                for (int kic_p = 0; kic_p < Kic.get(i, c); kic_p++) {
                                    if (kicIdx != kic_p) {
                                        double rate =
                                                proc.get(sn.stations.get(i))
                                                        .get(sn.jobclasses.get(c))
                                                        .get(0)
                                                        .get(kicIdx, kic_p);
                                        // Use softmin instead of min for smooth approximation
                                        double softminFactor = softmin(ni, sn.nservers.get(i, 0), alpha) / ni;
                                        if (ni > 0) {
                                            dxdt.set(
                                                    xic + kicIdx,
                                                    0,
                                                    dxdt.get(xic + kicIdx, 0) - (x[xic + kicIdx] * rate * softminFactor));
                                            dxdt.set(
                                                    xic + kic_p,
                                                    0,
                                                    dxdt.get(xic + kic_p, 0) + (x[xic + kicIdx] * rate * softminFactor));
                                        }
                                    }
                                }
                            }
                        }
                    }
                    // Service completions
                    for (int c = 0; c < K; c++) {
                        if (enabled[i][c]) {
                            int xic = (int) qIndices.get(i, c);
                            for (int j = 0; j < M; j++) {
                                for (int l = 0; l < K; l++) {
                                    int xjl = (int) qIndices.get(j, l);
                                    if (enabled[j][l]) {
                                        Matrix D0 = proc.get(sn.stations.get(j)).get(sn.jobclasses.get(l)).get(0);
                                        Matrix D1 = proc.get(sn.stations.get(j)).get(sn.jobclasses.get(l)).get(1);
                                        Matrix pie = map_pie(D0, D1);

                                        if (rt.get(i * K + c, j * K + l) > 0) {
                                            for (int kicIdx = 0; kicIdx < Kic.get(i, c); kicIdx++) {
                                                for (int kjl = 0; kjl < Kic.get(j, l); kjl++) {
                                                    double rate =
                                                            phi.get(sn.stations.get(i)).get(sn.jobclasses.get(c)).get(kicIdx, 0)
                                                                    * mu.get(sn.stations.get(i))
                                                                    .get(sn.jobclasses.get(c))
                                                                    .get(kicIdx, 0)
                                                                    * rt.get(i * K + c, j * K + l)
                                                                    * pie.get(0, kjl);
                                                    // Use softmin instead of min for capacity constraint
                                                    rate *= softmin(ni, sn.nservers.get(i, 0), alpha) / ni;
                                                    dxdt.set(
                                                            xic + kicIdx,
                                                            0,
                                                            dxdt.get(xic + kicIdx, 0) - (x[xic + kicIdx] * rate));
                                                    dxdt.set(xjl + kjl, 0, dxdt.get(xjl + kjl, 0) + (x[xic + kicIdx] * rate));
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    break;

                case FCFS:
                    idxIni = (int) qIndices.get(i, 0);
                    idxEnd = (int) qIndices.get(i, K - 1) + (int) Kic.get(i, K - 1);
                    ni = 0;
                    for (int idx = idxIni; idx < idxEnd; idx++) {
                        ni += x[idx];
                    }
                    double wni = GlobalConstants.FineTol;
                    for (int c = 0; c < K; c++) {
                        for (int kicIdx = 0; kicIdx < Kic.get(i, c); kicIdx++) {
                            if (enabled[i][c]) {
                                int xic = (int) qIndices.get(i, c);
                                w.set(
                                        c,
                                        kicIdx,
                                        -1
                                                / proc.get(sn.stations.get(i))
                                                .get(sn.jobclasses.get(c))
                                                .get(0)
                                                .get(kicIdx, kicIdx));
                                wni += w.get(c, kicIdx) * x[xic + kicIdx];
                            }
                        }
                    }
                    // Phase changes
                    for (int c = 0; c < K; c++) {
                        if (enabled[i][c]) {
                            int xic = (int) qIndices.get(i, c);
                            for (int kicIdx = 0; kicIdx < Kic.get(i, c) - 1; kicIdx++) {
                                for (int kic_p = 0; kic_p < Kic.get(i, c); kic_p++) {
                                    if (kicIdx != kic_p) {
                                        double rate =
                                                proc.get(sn.stations.get(i))
                                                        .get(sn.jobclasses.get(c))
                                                        .get(0)
                                                        .get(kicIdx, kic_p)
                                                        * softmin(ni, sn.nservers.get(i, 0), alpha) // Use softmin here
                                                        * w.get(c, kicIdx)
                                                        / wni;
                                        dxdt.set(xic + kicIdx, 0, dxdt.get(xic + kicIdx, 0) - (x[xic + kicIdx] * rate));
                                        dxdt.set(xic + kic_p, 0, dxdt.get(xic + kic_p, 0) + (x[xic + kicIdx] * rate));
                                    }
                                }
                            }
                        }
                    }
                    // Service completions
                    for (int c = 0; c < K; c++) {
                        if (enabled[i][c]) {
                            int xic = (int) qIndices.get(i, c);
                            for (int j = 0; j < M; j++) {
                                for (int l = 0; l < K; l++) {
                                    int xjl = (int) qIndices.get(j, l);
                                    if (enabled[j][l]) {
                                        Matrix D0 = proc.get(sn.stations.get(j)).get(sn.jobclasses.get(l)).get(0);
                                        Matrix D1 = proc.get(sn.stations.get(j)).get(sn.jobclasses.get(l)).get(1);
                                        Matrix pie = map_pie(D0, D1);

                                        if (rt.get(i * K + c, j * K + l) > 0) {
                                            for (int kicIdx = 0; kicIdx < Kic.get(i, c); kicIdx++) {
                                                for (int kjl = 0; kjl < Kic.get(j, l); kjl++) {
                                                    double rate =
                                                            phi.get(sn.stations.get(i)).get(sn.jobclasses.get(c)).get(kicIdx, 0)
                                                                    * mu.get(sn.stations.get(i))
                                                                    .get(sn.jobclasses.get(c))
                                                                    .get(kicIdx, 0)
                                                                    * rt.get(i * K + c, j * K + l)
                                                                    * pie.get(0, kjl)
                                                                    * softmin(ni, sn.nservers.get(i, 0), alpha) // Use softmin here
                                                                    * w.get(c, kicIdx)
                                                                    / wni;
                                                    dxdt.set(
                                                            xic + kicIdx,
                                                            0,
                                                            dxdt.get(xic + kicIdx, 0) - (x[xic + kicIdx] * rate));
                                                    dxdt.set(xjl + kjl, 0, dxdt.get(xjl + kjl, 0) + (x[xic + kicIdx] * rate));
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    break;

                case DPS:
                    double sumWI = w.sumRows(i);
                    for (int col = 0; col < w.getNumCols(); col++) {
                        w.set(i, col, w.get(i, col) / sumWI);
                    }
                    idxIni = (int) qIndices.get(i, 0);
                    idxEnd = (int) qIndices.get(i, K - 1) + (int) Kic.get(i, K - 1);
                    wni = 0;
                    for (int idx = idxIni; idx < idxEnd; idx++) {
                        wni += x[idx];
                    }
                    // Phase changes
                    for (int c = 0; c < K; c++) {
                        if (enabled[i][c]) {
                            int xic = (int) qIndices.get(i, c);
                            for (int kicIdx = 0; kicIdx < Kic.get(i, c) - 1; kicIdx++) {
                                for (int kic_p = 0; kic_p < Kic.get(i, c); kic_p++) {
                                    if (kicIdx != kic_p) {
                                        double rate =
                                                proc.get(sn.stations.get(i))
                                                        .get(sn.jobclasses.get(c))
                                                        .get(0)
                                                        .get(kicIdx, kic_p);
                                        // Use softmin for DPS scheduling
                                        double softminFactor = softmin(wni, sn.nservers.get(i, 0), alpha) * w.get(c, kicIdx) / wni;
                                        dxdt.set(
                                                xic + kicIdx,
                                                0,
                                                dxdt.get(xic + kicIdx, 0) - (x[xic + kicIdx] * rate * softminFactor));
                                        dxdt.set(
                                                xic + kic_p,
                                                0,
                                                dxdt.get(xic + kic_p, 0) + (x[xic + kicIdx] * rate * softminFactor));
                                    }
                                }
                            }
                        }
                    }
                    // Service completions
                    for (int c = 0; c < K; c++) {
                        if (enabled[i][c]) {
                            int xic = (int) qIndices.get(i, c);
                            for (int j = 0; j < M; j++) {
                                for (int l = 0; l < K; l++) {
                                    int xjl = (int) qIndices.get(j, l);
                                    if (enabled[j][l]) {
                                        Matrix D0 = proc.get(sn.stations.get(j)).get(sn.jobclasses.get(l)).get(0);
                                        Matrix D1 = proc.get(sn.stations.get(j)).get(sn.jobclasses.get(l)).get(1);
                                        Matrix pie = map_pie(D0, D1);

                                        if (rt.get(i * K + c, j * K + l) > 0) {
                                            for (int kicIdx = 0; kicIdx < Kic.get(i, c); kicIdx++) {
                                                for (int kjl = 0; kjl < Kic.get(j, l); kjl++) {
                                                    double rate =
                                                            phi.get(sn.stations.get(i)).get(sn.jobclasses.get(c)).get(kicIdx, 0)
                                                                    * mu.get(sn.stations.get(i))
                                                                    .get(sn.jobclasses.get(c))
                                                                    .get(kicIdx, 0)
                                                                    * rt.get(i * K + c, j * K + l)
                                                                    * pie.get(0, kjl);
                                                    // Use softmin for capacity constraint
                                                    rate *= softmin(wni, sn.nservers.get(i, 0), alpha) * w.get(c, kicIdx) / wni;
                                                    dxdt.set(
                                                            xic + kicIdx,
                                                            0,
                                                            dxdt.get(xic + kicIdx, 0) - (x[xic + kicIdx] * rate));
                                                    dxdt.set(xjl + kjl, 0, dxdt.get(xjl + kjl, 0) + (x[xic + kicIdx] * rate));
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    break;
                    
                default:
                    // For other scheduling strategies, fall back to the state-dependent method
                    line_error(mfilename(new Object[]{}),
                            "Softmin ODE method does not support scheduling strategy: " + sn.sched.get(sn.stations.get(i)) + 
                            ". Try with statedep method.");
                    break;
            }
        }

        return dxdt;
    }

    @Override
    public void computeDerivatives(double t, double[] x, double[] dxdt)
            throws MaxCountExceededException, DimensionMismatchException {

        // Clamp state to non-negative before computing derivatives.
        // Matches MATLAB's ode15s NonNegative option: prevents negative fluid
        // levels from corrupting rate computations in the ODE dynamics.
        for (int idx = 0; idx < x.length; idx++) {
            if (x[idx] < 0) {
                x[idx] = 0;
            }
        }

        // Use precomputed structures (w is copied since some methods modify it in-place)
        Matrix w = cachedW.copy();

        Matrix dxdtTmp;
        if (Objects.equals(options.method, "statedep")) {
            dxdtTmp = calculatedxdtStateDepMethod(x, cachedEnabled, cachedQIndices, cachedKic, w);
        } else if (Objects.equals(options.method, "softmin")) {
            dxdtTmp = calculatedxdtSoftminMethod(x, cachedEnabled, cachedQIndices, cachedKic, w);
        } else {
            // Use precomputed (and potentially immediate-eliminated) allJumps/rateBase/eventIdx
            dxdtTmp =
                    calculatedxdtClosingMethod(x, w, cachedEnabled, cachedQIndices, cachedKic,
                            cachedAllJumps, cachedRateBase, cachedEventIdx);
        }

        for (int i = 0; i < dxdt.length; i++) {
            dxdt[i] = dxdtTmp.get(i);
        }
    }

    @Override
    public int getDimension() {
        return numDimensions;
    }

    private void setNextJump(Matrix jumps, int completionIdx, int startIdx) {

        int jumpsCols = jumps.getNumCols();
        jumps.expandMatrix(jumps.getNumRows(), jumpsCols + 1, jumps.getNumElements() + 2);
        jumps.set(completionIdx, jumpsCols, -1); // type c in stat i completes service
        jumps.set(startIdx, jumpsCols, 1); // type c job starts in stat j
    }
}
