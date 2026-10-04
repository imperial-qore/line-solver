/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.sn;

import java.util.Map;
import java.util.function.DoubleUnaryOperator;

import jline.api.mam.Map_pie;
import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.constant.ProcessType;
import jline.lang.nodes.Station;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * The arrival rate of a station-class pair AS A FUNCTION OF TIME.
 *
 * <p>The time-varying analyses (Mt/G/inf, the modified offered load, the
 * Gt/Mt/st+GI fluid queue) consume lambda(t) itself, not a mean rate: their
 * whole content is the LAG between when work arrives and when it is felt, and a
 * time-averaged rate has no lag. LINE carries a time-varying arrival as an NHPP
 * or a MAPt, whose {@code sn.proc} slot is a piecewise-constant schedule, so
 * lambda(t) is read off the segment in force at t.
 *
 * <p>For any other process the rate is constant and the handle returns it, which
 * is what lets a caller ask for the time-varying analysis of a stationary model
 * and get the stationary answer rather than an error.
 *
 * <p>Java twin of {@code matlab/src/api/sn/sn_arrival_rate_fun.m}.
 *
 * @since LINE 3.1.0
 */
public final class SnArrivalRateFun {

    private SnArrivalRateFun() {}

    /** lambda(t), with whether it actually varies and the cycle length. */
    public static final class RateFun {
        /** The rate as a function of time. */
        public final DoubleUnaryOperator lambda;
        /** Whether the rate depends on t at all. */
        public final boolean timeVarying;
        /** The cycle length when the schedule is cyclic, infinite otherwise. */
        public final double period;

        RateFun(DoubleUnaryOperator lambda, boolean timeVarying, double period) {
            this.lambda = lambda;
            this.timeVarying = timeVarying;
            this.period = period;
        }
    }

    /**
     * Build lambda(t) for station {@code ist}, class {@code r}.
     *
     * @param sn  the network struct
     * @param ist station index
     * @param r   class index
     * @return the rate function
     */
    public static RateFun snArrivalRateFun(NetworkStruct sn, int ist, int r) {
        Station station = sn.stations.get(ist);
        JobClass jobClass = sn.jobclasses.get(r);
        final double rate = sn.rates.get(ist, r);
        ProcessType ty = null;
        if (sn.procid != null && sn.procid.get(station) != null) {
            ty = sn.procid.get(station).get(jobClass);
        }
        boolean isNhpp = ty == ProcessType.NHPP;
        boolean isSched = isNhpp || ty == ProcessType.MAPT || ty == ProcessType.PHT;
        if (!isSched) {
            return new RateFun(new DoubleUnaryOperator() {
                @Override
                public double applyAsDouble(double t) {
                    return rate;
                }
            }, false, Double.POSITIVE_INFINITY);
        }

        Map<JobClass, MatrixCell> procMap = sn.proc.get(station);
        MatrixCell slot = procMap == null ? null : procMap.get(jobClass);
        if (slot == null || slot.size() < 3) {
            throw new RuntimeException("snArrivalRateFun: the schedule slot of sn.proc is malformed");
        }
        final double[] bp = rowOf(slot.get(0));
        final double[] segRate;
        final boolean cyclic;
        if (isNhpp) {
            // An NHPP slot is {breakpoints, rates, cyclic}: the rates ARE
            // lambda(t), one per interval, so there is no MAP pair to reduce.
            segRate = rowOf(slot.get(1));
            cyclic = slot.get(2).get(0) != 0;
        } else {
            // A MAPt/PHt slot is flat, [breakpoints, A_1..A_n, B_1..B_n, cyclic].
            int n = (slot.size() - 2) / 2;
            segRate = new double[n];
            for (int k = 0; k < n; k++) {
                Matrix a = slot.get(1 + k);
                Matrix b = slot.get(1 + n + k);
                Matrix D0;
                Matrix D1;
                if (ty == ProcessType.MAPT) {
                    D0 = a;
                    D1 = b;
                } else {
                    // PHt: a is the alpha row, b the sub-generator S, and the
                    // equivalent MAP pair is (S, s*alpha) with s = -S e.
                    D0 = b;
                    D1 = new Matrix(b.getNumRows(), b.getNumCols());
                    for (int i = 0; i < b.getNumRows(); i++) {
                        double exit = 0;
                        for (int j = 0; j < b.getNumCols(); j++) {
                            exit -= b.get(i, j);
                        }
                        for (int j = 0; j < b.getNumCols(); j++) {
                            D1.set(i, j, exit * a.get(0, j));
                        }
                    }
                }
                if (D0.getNumRows() == 1) {
                    segRate[k] = D1.get(0, 0);
                } else {
                    // The arrival rate of a segment is pie_k D1_k e, the
                    // stationary throughput of that segment's own MAP.
                    Matrix pie = Map_pie.map_pie(D0, D1);
                    double v = 0;
                    for (int i = 0; i < D1.getNumRows(); i++) {
                        for (int j = 0; j < D1.getNumCols(); j++) {
                            v += pie.get(i) * D1.get(i, j);
                        }
                    }
                    segRate[k] = v;
                }
            }
            cyclic = slot.get(slot.size() - 1).get(0) != 0;
        }

        boolean varying = false;
        for (int k = 1; k < segRate.length; k++) {
            if (Math.abs(segRate[k] - segRate[0]) > 1e-12) {
                varying = true;
                break;
            }
        }
        final boolean cyc = cyclic;
        DoubleUnaryOperator lam = new DoubleUnaryOperator() {
            @Override
            public double applyAsDouble(double t) {
                double u = t;
                if (cyc && bp[bp.length - 1] > bp[0]) {
                    double span = bp[bp.length - 1] - bp[0];
                    u = bp[0] + ((t - bp[0]) % span + span) % span;
                }
                // Segment k is in force on [bp[k], bp[k+1]). Before the first
                // breakpoint the first segment holds and after the last the last
                // one does, so a caller integrating over an infinite past (the
                // Mt/G/inf convolution) gets a defined rate everywhere.
                int k = 0;
                for (int i = 0; i + 1 < bp.length; i++) {
                    if (u >= bp[i]) {
                        k = i;
                    }
                }
                if (k >= segRate.length) {
                    k = segRate.length - 1;
                }
                return segRate[k];
            }
        };
        double period = cyclic ? bp[bp.length - 1] - bp[0] : Double.POSITIVE_INFINITY;
        return new RateFun(lam, varying, period);
    }

    private static double[] rowOf(Matrix m) {
        int n = m.getNumRows() * m.getNumCols();
        double[] out = new double[n];
        int idx = 0;
        for (int i = 0; i < m.getNumRows(); i++) {
            for (int j = 0; j < m.getNumCols(); j++) {
                out[idx++] = m.get(i, j);
            }
        }
        return out;
    }
}
