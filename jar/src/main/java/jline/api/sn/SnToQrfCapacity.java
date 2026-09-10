/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.sn;

import jline.lang.NetworkStruct;

/**
 * Per-station occupancy bound F(i) for the QRF bounds.
 *
 * <p>Port of {@code matlab/src/api/sn/sn_to_qrf_capacity.m}.
 *
 * <p>The QRF bounds index every marginal by 0..F(i), so F is an occupancy bound
 * rather than a declared capacity: the station's buffer where that buffer
 * BINDS, and the population N everywhere else, since no queue of a closed model
 * can hold more than N jobs.
 *
 * <p>Binding is decided by {@link SnGetBufferSize}, the single place in LINE
 * that makes that call: refreshCapacity derives a finite classcap (the chain
 * population) at every station of every closed model, so a plain finiteness
 * test on {@code sn.cap} would report a buffer at every station.
 *
 * <p>Both QRF blocking bounds need this. {@code qrf.bas} needs it beside the
 * blocking tables {@link SnToQrfBlocking} derives; {@code qrf.rsrd} needs it
 * ALONE, since its PBB constraint reads only which queues can be full and it
 * carries no blocking tables at all.
 *
 * @since LINE 3.0
 */
public final class SnToQrfCapacity {

    private SnToQrfCapacity() {
    }

    /** Occupancy bounds, which of them bind, and why not when they cannot be built. */
    public static final class Result {
        /** (nstations) occupancy bound of each station, in jobs. */
        public final int[] F;
        /** (nstations) true where the buffer can refuse a job. */
        public final boolean[] binding;
        /** Empty on success, otherwise why F is not defined for this model. */
        public final String msg;

        Result(int[] F, boolean[] binding, String msg) {
            this.F = F;
            this.binding = binding;
            this.msg = msg;
        }
    }

    /**
     * @param sn network structure
     * @return the occupancy bounds, or a Result carrying a non-empty msg
     */
    public static Result snToQrfCapacity(NetworkStruct sn) {
        int M = sn.nstations;
        int[] F = new int[M];
        boolean[] binding = new boolean[M];

        double Nd = 0.0;
        if (sn.njobs != null) {
            Nd = sn.njobs.elementSum();
        }
        if (!Double.isFinite(Nd) || Nd < 1) {
            return new Result(F, binding,
                    "the QRF bounds need a closed model with a finite population.");
        }
        int N = (int) Math.round(Nd);

        for (int i = 0; i < M; i++) {
            double b = SnGetBufferSize.snGetBufferSize(sn, i);
            binding[i] = Double.isFinite(b);
            F[i] = (!Double.isFinite(b) || b > N) ? N : (int) Math.round(b);
            if (F[i] < 1) {
                return new Result(F, binding, String.format(
                        "station %d has capacity %d: the QRF bounds need every queue to be able "
                                + "to hold at least one job.", i + 1, F[i]));
            }
        }
        return new Result(F, binding, "");
    }
}
