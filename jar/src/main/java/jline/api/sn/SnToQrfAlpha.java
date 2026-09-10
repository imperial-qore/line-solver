/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.sn;

import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Station;
import jline.util.matrix.MatrixCell;

/**
 * Per-station load-dependent rate scaling alpha(i,n) for the QRF bounds.
 *
 * <p>Port of {@code matlab/src/api/sn/sn_to_qrf_alpha.m}.
 *
 * <p>The load-dependent QRF arms carry a scaling {@code alpha(i,n)} that
 * multiplies EVERY rate out of station i while it holds n jobs, completions mu
 * and background phase changes v alike -- see the q construction in
 * {@link jline.api.mapqn.Mapqn_qrf_noblo_mmi_ld}. That is exactly the rate law
 * of
 *
 * <pre>
 *   an infinite server        alpha(i,n) = n
 *   a c-server station        alpha(i,n) = min(n, c_i)
 *   limited load dependence   alpha(i,n) = sn.lldscaling(i,n)
 * </pre>
 *
 * so the three COMPOSE BY MULTIPLICATION and not one of them is an
 * approximation: the relaxed chain is the model's own, and the QRF answer keeps
 * whatever status it had on a single-server model.
 *
 * <p>WHERE IT STOPS BEING THE MODEL'S OWN IS PHASE-TYPE SERVICE AT A STATION
 * THAT SERVES SEVERAL JOBS AT ONCE. The QRF local state carries ONE phase per
 * station, a faithful description of one job in service and of nothing else:
 * min(n,c) jobs served in parallel each advance through a phase of their own,
 * and no scaling of a single-phase process reproduces that joint motion. A
 * multiserver or delay station must therefore be exponential -- scaling a PH
 * server by min(n,c) would answer a DIFFERENT chain, so the relaxation would
 * stop containing the model's stationary distribution and the number would
 * bound nothing. Limited load dependence at a SINGLE server is exempt and
 * admits PH freely: one job is in service whatever the rate.
 *
 * <p>THE UTILIZATION NORMALIZER IS THE DECLARED PEAK, NOT max(alpha). LINE
 * reports U = T*S/peak at every station whose rate scales with the population,
 * one convention shared by multiserver, lld and class dependence. {@code peak}
 * is therefore nservers(i) times the largest lld scaling the model can REACH,
 * and not max(alpha(i,:)): at c = 3 with N = 2 the reachable alpha peaks at 2
 * while the station still has three servers, and normalizing by 2 would report
 * a utilization the model never attains. Infinite at a delay, where LINE
 * reports U = QN instead.
 *
 * @since LINE 3.0
 */
public final class SnToQrfAlpha {

    private SnToQrfAlpha() {
    }

    /** The scaling, the utilization normalizer, and why they may not exist. */
    public static final class Result {
        /** (nstations x N) rate scaling at population n = 1..N. */
        public final double[][] alpha;
        /** Empty on success, otherwise why alpha is not defined for this model. */
        public final String msg;
        /** True when alpha is not identically 1, i.e. the model needs a load-dependent arm. */
        public final boolean ld;
        /** (nstations) utilization normalizer; infinite at a delay. */
        public final double[] peak;

        Result(double[][] alpha, String msg, boolean ld, double[] peak) {
            this.alpha = alpha;
            this.msg = msg;
            this.ld = ld;
            this.peak = peak;
        }
    }

    /**
     * @param sn network structure
     * @return the scaling, or a Result carrying a non-empty msg
     */
    public static Result snToQrfAlpha(NetworkStruct sn) {
        int M = sn.nstations;
        double[] peak = new double[M];
        for (int i = 0; i < M; i++) {
            peak[i] = 1.0;
        }

        double Nd = sn.njobs != null ? sn.njobs.elementSum() : 0.0;
        if (!Double.isFinite(Nd) || Nd < 1) {
            return new Result(ones(M, 1), "the QRF bounds need a closed model with a finite "
                    + "population.", false, peak);
        }
        int N = (int) Math.round(Nd);

        double[][] alpha = ones(M, N);
        int smax = (sn.lldscaling != null && sn.lldscaling.getNumRows() > 0)
                ? sn.lldscaling.getNumCols() : 0;
        boolean ld = false;

        for (int i = 0; i < M; i++) {
            double c = sn.nservers.get(i, 0);
            Station station = sn.stations.get(i);
            boolean isDelay = Double.isInfinite(c)
                    || sn.sched.get(station) == SchedStrategy.INF;
            boolean servesMany = isDelay || c > 1;
            int ki = phases(sn, i);
            if (servesMany && ki > 1) {
                // ld stays TRUE through the refusal: the model IS load dependent, and
                // the caller has to tell "no arm serves this" from "the arm you asked
                // for does not". Clearing it here would report the latter for both.
                return new Result(alpha, String.format(
                        "station %d serves %s jobs at once with %d-phase service, and the QRF "
                                + "local state carries one phase per station, which describes one "
                                + "job in service and no more. Give that station exponential "
                                + "service, or use a single-server model.",
                        i + 1, isDelay ? "unboundedly many" : "up to " + (int) c, ki),
                        true, peak);
            }
            double lldpeak = 1.0;
            for (int n = 1; n <= N; n++) {
                if (isDelay) {
                    alpha[i][n - 1] = n;
                } else if (c > 1) {
                    alpha[i][n - 1] = Math.min(n, c);
                }
                if (smax > 0) {
                    double s = sn.lldscaling.get(i, Math.min(n, smax) - 1);
                    lldpeak = Math.max(lldpeak, s);
                    alpha[i][n - 1] *= s;
                }
                if (alpha[i][n - 1] != 1.0) {
                    ld = true;
                }
            }
            peak[i] = isDelay ? Double.POSITIVE_INFINITY : c * lldpeak;
        }
        return new Result(alpha, "", ld, peak);
    }

    /**
     * Phases of station i's service process, read from {@code sn.proc} as the
     * QRF adapter reads them: that {D0, D1} pair is what sizes the local state,
     * so testing it keeps the refusal and the formulation on one quantity.
     */
    private static int phases(NetworkStruct sn, int i) {
        if (sn.proc == null || sn.stations == null || sn.jobclasses == null) {
            return 1;
        }
        Station station = sn.stations.get(i);
        if (sn.proc.get(station) == null) {
            return 1;
        }
        MatrixCell procI = sn.proc.get(station).get(sn.jobclasses.get(0));
        if (procI == null || procI.size() == 0 || procI.get(0) == null) {
            return 1;
        }
        return procI.get(0).getNumRows();
    }

    private static double[][] ones(int rows, int cols) {
        double[][] a = new double[rows][cols];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                a[i][j] = 1.0;
            }
        }
        return a;
    }
}
