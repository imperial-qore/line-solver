/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.fluid;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.constant.ProcessType;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Source;
import jline.lang.nodes.Station;
import jline.lang.processes.Distribution;
import jline.lang.processes.NHPP;
import jline.solvers.fluid.handlers.FluidRateMultiplier;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * Non-homogeneous Poisson (NHPP) support for the fluid solver.
 *
 * <p>An NHPP stores its process representation as the {breakpoints, rates,
 * cyclic} schedule rather than a {D0, D1} MAP, and its {@code pie} is NaN. The
 * fluid model treats such a source as a single-phase exponential at its nominal
 * (time-average) rate, which {@code sn.mu} already carries; the time-varying
 * intensity lambda(t) is applied separately by the per-event rate multiplier
 * ({@link FluidRateMultiplier}). Mirrors the NHPP handling of MATLAB
 * {@code solver_fluid.m} and {@code SolverFLD/getTranAvg.m}.</p>
 */
public final class FluidNhpp {

    private FluidNhpp() {
    }

    /**
     * Returns the process map to hand to the closing-rate builder, with every
     * NHPP station-class replaced by the equivalent one-phase exponential MAP
     * {D0, D1} = {-lam, lam} at its nominal rate {@code mu(i,k)(0)}, so that
     * {@code map_pie} sees a valid Markovian representation. The corresponding
     * {@code phi} entry is set to 1.
     *
     * <p>The input maps are not modified except for {@code phi}, which the
     * callers own as a local working copy; when no NHPP is present the original
     * {@code sn.proc} reference is returned unchanged.</p>
     *
     * @param sn  network structure (read for procid and the station/class lists)
     * @param mu  per-(station,class) phase rate vectors of the fluid ODE
     * @param phi per-(station,class) completion probability vectors (updated in place)
     * @return the process map with the NHPP substitutions applied
     */
    public static Map<Station, Map<JobClass, MatrixCell>> substituteNhppProc(
            NetworkStruct sn,
            Map<Station, Map<JobClass, Matrix>> mu,
            Map<Station, Map<JobClass, Matrix>> phi) {

        if (sn.procid == null) {
            return sn.proc;
        }
        Map<Station, Map<JobClass, MatrixCell>> proc = null;
        for (int i = 0; i < sn.nstations; i++) {
            Station station = sn.stations.get(i);
            Map<JobClass, ProcessType> procidI = sn.procid.get(station);
            if (procidI == null) {
                continue;
            }
            for (int k = 0; k < sn.nclasses; k++) {
                JobClass jobClass = sn.jobclasses.get(k);
                if (procidI.get(jobClass) != ProcessType.NHPP) {
                    continue;
                }
                Matrix muIK = (mu.get(station) == null) ? null : mu.get(station).get(jobClass);
                if (muIK == null || muIK.isEmpty() || Double.isNaN(muIK.get(0, 0))) {
                    continue;
                }
                double lam = muIK.get(0, 0);
                if (proc == null) {
                    proc = new HashMap<Station, Map<JobClass, MatrixCell>>();
                    for (int ii = 0; ii < sn.nstations; ii++) {
                        Station st = sn.stations.get(ii);
                        Map<JobClass, MatrixCell> src = sn.proc.get(st);
                        proc.put(st, (src == null) ? new HashMap<JobClass, MatrixCell>()
                                : new HashMap<JobClass, MatrixCell>(src));
                    }
                }
                Matrix d0 = new Matrix(1, 1, 1);
                d0.set(0, 0, -lam);
                Matrix d1 = new Matrix(1, 1, 1);
                d1.set(0, 0, lam);
                proc.get(station).put(jobClass, new MatrixCell(d0, d1));

                if (phi != null && phi.get(station) != null) {
                    Matrix phiOne = new Matrix(1, 1, 1);
                    phiOne.set(0, 0, 1.0);
                    phi.get(station).put(jobClass, phiOne);
                }
            }
        }
        return (proc == null) ? sn.proc : proc;
    }

    /**
     * Builds the NHPP schedule entries for every EXT/source station carrying a
     * non-homogeneous arrival process. Returns an empty list when there is
     * none. Station and class indices are in the sn index space.
     *
     * @param sn network structure whose stations are inspected
     * @return the schedule entries to inject in {@code options.config.nhpp_sched}
     */
    public static List<FluidRateMultiplier.NhppEntry> detectNhppSchedule(NetworkStruct sn) {
        List<FluidRateMultiplier.NhppEntry> sched = new ArrayList<FluidRateMultiplier.NhppEntry>();
        if (sn == null || sn.stations == null) {
            return sched;
        }
        for (int i = 0; i < sn.nstations; i++) {
            Station station = sn.stations.get(i);
            if (sn.sched.get(station) != SchedStrategy.EXT || !(station instanceof Source)) {
                continue;
            }
            Source source = (Source) station;
            for (int c = 0; c < sn.nclasses; c++) {
                JobClass jobClass = sn.jobclasses.get(c);
                Distribution distr = source.getArrivalDistribution(jobClass);
                if (distr instanceof NHPP) {
                    sched.add(new FluidRateMultiplier.NhppEntry(i, c, (NHPP) distr));
                }
            }
        }
        return sched;
    }
}
