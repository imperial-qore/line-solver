/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.sn;

import jline.lang.NetworkStruct;

/**
 * Physical buffer size of a station, in jobs, the one in service included.
 *
 * <p>Port of {@code matlab/src/api/sn/sn_get_buffer_size.m}. Kendall's K: the
 * total occupancy bound of station IST, taken as the tighter of the station
 * capacity {@code sn.cap(ist)} and the per-class capacities
 * {@code sn.classcap(ist,:)}. Returns +Inf when the station is unbounded. Both
 * fields are populated by refreshCapacity, which already folds setCapacity,
 * setClassCapacity, a finite orbit and the closed-chain population into them,
 * so this is the single place that decides whether a buffer BINDS.
 *
 * <p>Only a buffer that can actually BIND is reported. refreshCapacity derives a
 * FINITE classcap, the chain population, for EVERY closed model, so a plain
 * finiteness test reports a buffer at every station of every closed model; a
 * capacity at least as large as the total population can never refuse a job and
 * is returned as +Inf. sum(njobs) is +Inf as soon as one class is open, so any
 * finite capacity reachable by an open class binds.
 *
 * <p>The population counted is the one that can REACH this station, not the
 * model's total: using every class made a MIXED model read as finite-buffered at
 * stations an open class never visits.
 */
public final class SnGetBufferSize {

    private SnGetBufferSize() {
    }

    /**
     * @param sn  network structure
     * @param ist station index
     * @return buffer size in jobs, +Inf if unbounded or non-binding
     */
    public static double snGetBufferSize(NetworkStruct sn, int ist) {
        double N = Double.POSITIVE_INFINITY;
        if (sn.cap != null && sn.cap.length() > ist) {
            double c = unbounded(sn.cap.get(ist));
            if (c >= 0) {
                N = Math.min(N, c);
            }
        }
        boolean[] reach = null;
        int nclasses = sn.nclasses;
        if (sn.classcap != null && sn.classcap.getNumRows() > ist) {
            double sumccap = 0.0;
            boolean any = false;
            reach = new boolean[sn.classcap.getNumCols()];
            for (int r = 0; r < sn.classcap.getNumCols(); r++) {
                double v = unbounded(sn.classcap.get(ist, r));
                // a zero marks a class that is not served here
                reach[r] = v > 0;
                if (v > 0) {
                    sumccap += v;
                    any = true;
                }
            }
            if (any) {
                N = Math.min(N, sumccap);
            }
        }
        if (sn.njobs != null && sn.njobs.length() > 0) {
            double reachableJobs = 0.0;
            boolean anyReach = false;
            for (int r = 0; r < nclasses && r < sn.njobs.length(); r++) {
                boolean inReach = (reach == null) || (r < reach.length && reach[r]);
                if (inReach) {
                    anyReach = true;
                    reachableJobs += sn.njobs.get(r);
                }
            }
            if (!anyReach) {
                reachableJobs = 0.0;
            }
            if (N >= reachableJobs) {
                // declared but unreachable: the buffer can never refuse a job
                N = Double.POSITIVE_INFINITY;
            }
        }
        return N;
    }

    /**
     * MATLAB writes an unbounded capacity as Inf; the JAR writes
     * Integer.MAX_VALUE into cap/classcap/chaincap. Fold the sentinel so the
     * predicate reads the same on both sides.
     */
    private static double unbounded(double v) {
        return (v >= (double) jline.GlobalConstants.MaxInt) ? Double.POSITIVE_INFINITY : v;
    }
}
