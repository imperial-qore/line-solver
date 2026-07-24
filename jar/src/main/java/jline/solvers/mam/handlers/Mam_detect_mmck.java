/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.mam.handlers;

import java.util.ArrayList;
import java.util.List;

import jline.lang.NetworkStruct;
import jline.lang.constant.ProcessType;
import jline.util.matrix.MatrixCell;

/**
 * Decides if a station matches the M/M/c/K assumptions.
 *
 * <p>Port of matlab/src/solvers/MAM/mam_detect_mmck.m.</p>
 */
public final class Mam_detect_mmck {
    private Mam_detect_mmck() {}

    /**
     * Outcome of the M/M/c/K detection.
     */
    public static final class Result {
        /** true if the M/M/c/K closed form is exact at this station. */
        public final boolean isMmck;
        /** the shared service rate, NaN when isMmck is false. */
        public final double muRate;

        Result(boolean isMmck, double muRate) {
            this.isMmck = isMmck;
            this.muRate = muRate;
        }
    }

    /**
     * Returns isMmck=true (and the shared service rate) only when the aggregated
     * arrival MMAP is single-phase (Poisson superposition), every active class has
     * Exp service at the station, and all active classes share the same rate.
     */
    public static Result mam_detect_mmck(NetworkStruct sn, int ist, int K, MatrixCell mmapNode) {
        // arrivals must be a single-phase Poisson superposition
        if (mmapNode == null || mmapNode.isEmpty() || mmapNode.get(0).getNumRows() != 1) {
            return new Result(false, Double.NaN);
        }

        List<Double> muVals = new ArrayList<Double>();
        for (int k = 0; k < K; k++) {
            ProcessType pt = null;
            if (sn.procid != null && sn.procid.get(sn.stations.get(ist)) != null) {
                pt = sn.procid.get(sn.stations.get(ist)).get(sn.jobclasses.get(k));
            }
            double rate = sn.rates.get(ist, k);
            if (pt != ProcessType.EXP) {
                // Disabled (NaN rate) classes are skipped; everything else must be Exp
                if (Double.isNaN(rate)) {
                    continue;
                }
                return new Result(false, Double.NaN);
            }
            if (Double.isNaN(rate) || rate <= 0) {
                continue; // no inflow for this class
            }
            muVals.add(Double.valueOf(rate));
        }

        if (muVals.isEmpty()) {
            return new Result(false, Double.NaN);
        }

        double mx = muVals.get(0);
        double mn = muVals.get(0);
        for (int i = 1; i < muVals.size(); i++) {
            double v = muVals.get(i).doubleValue();
            if (v > mx) mx = v;
            if (v < mn) mn = v;
        }
        if (mx - mn > 1e-9 * Math.max(1.0, mx)) {
            return new Result(false, Double.NaN); // per-class service rates differ
        }

        return new Result(true, muVals.get(0).doubleValue());
    }
}
