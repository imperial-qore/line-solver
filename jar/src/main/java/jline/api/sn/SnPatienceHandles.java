/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.sn;

import java.util.Map;
import java.util.function.DoubleUnaryOperator;

import jline.api.mam.Map_cdf;
import jline.api.mam.Map_pdf;
import jline.api.qsys.Patience;
import jline.lang.NetworkStruct;
import jline.lang.constant.ImpatienceType;
import jline.lang.constant.ProcessType;
import jline.lang.JobClass;
import jline.lang.nodes.Station;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * Patience (time-to-abandon) handles derived from a {@link NetworkStruct}.
 *
 * <p>The abandonment solvers need the patience law as FUNCTIONS -- a
 * complementary cdf and a hazard rate -- not as moments, because that is what
 * the underlying theory consumes: Whitt's engineering solution reads the hazard
 * near the origin, and the fluid models integrate the ccdf. LINE stores the law
 * as a MAP/PH pair in {@code sn.impatienceProc}, from which both are available
 * in closed form.
 *
 * <p>Java twin of {@code matlab/src/api/sn/sn_patience_handles.m} and
 * {@code python/line_solver/api/sn/patience.py}.
 *
 * @since LINE 3.1.0
 */
public final class SnPatienceHandles {

    private SnPatienceHandles() {}

    /** The patience law of one station-class pair, in the forms the solvers consume. */
    public static final class Handles {
        /** F^c(t) = P(patience > t). */
        public final DoubleUnaryOperator ccdf;
        /** The patience density. */
        public final DoubleUnaryOperator pdf;
        /** The hazard rate h = f/(1-F). */
        public final DoubleUnaryOperator hazard;
        /** Mean patience, infinite when no rate is on record. */
        public final double mean;
        /** Whether the law is exponential, in which case the analysis is exact. */
        public final boolean isExponential;
        /** The abandonment rate 1/mean. */
        public final double rate;

        Handles(DoubleUnaryOperator ccdf, DoubleUnaryOperator pdf, DoubleUnaryOperator hazard,
                double mean, boolean isExponential, double rate) {
            this.ccdf = ccdf;
            this.pdf = pdf;
            this.hazard = hazard;
            this.mean = mean;
            this.isExponential = isExponential;
            this.rate = rate;
        }

        /** This law as the {@link Patience} argument of the abandonment solvers. */
        public Patience asPatience() {
            if (isExponential) {
                return Patience.exponential(rate);
            }
            return Patience.hazard(hazard);
        }
    }

    /**
     * Build the patience handles of station {@code ist}, class {@code r}.
     *
     * @param sn  the network struct
     * @param ist station index
     * @param r   class index
     * @return the handles, or {@code null} when the station-class pair has no
     *         reneging patience configured
     */
    public static Handles snPatienceHandles(NetworkStruct sn, int ist, int r) {
        if (sn == null || sn.impatienceClass == null) {
            return null;
        }
        if (sn.stations == null || ist < 0 || ist >= sn.stations.size()) {
            return null;
        }
        if (sn.jobclasses == null || r < 0 || r >= sn.jobclasses.size()) {
            return null;
        }
        Station station = sn.stations.get(ist);
        JobClass jobClass = sn.jobclasses.get(r);
        Map<JobClass, ImpatienceType> clsMap = sn.impatienceClass.get(station);
        if (clsMap == null || clsMap.get(jobClass) != ImpatienceType.RENEGING) {
            return null;
        }

        double rate = 0.0;
        if (sn.impatienceMu != null && sn.impatienceMu.get(station) != null) {
            Matrix muM = sn.impatienceMu.get(station).get(jobClass);
            if (muM != null && muM.length() > 0) {
                rate = muM.get(0);
            }
        }
        boolean isExp = false;
        if (sn.impatienceType != null && sn.impatienceType.get(station) != null) {
            isExp = sn.impatienceType.get(station).get(jobClass) == ProcessType.EXP;
        }

        MatrixCell pair = null;
        if (sn.impatienceProc != null && sn.impatienceProc.get(station) != null) {
            pair = sn.impatienceProc.get(station).get(jobClass);
        }

        if (pair == null || pair.size() < 2 || pair.get(0) == null || pair.get(1) == null) {
            if (rate <= 0) {
                return null;
            }
            // Only the rate is on record, so the law is exponential by construction.
            final double theta = rate;
            return new Handles(new DoubleUnaryOperator() {
                @Override
                public double applyAsDouble(double t) {
                    return Math.exp(-theta * t);
                }
            }, new DoubleUnaryOperator() {
                @Override
                public double applyAsDouble(double t) {
                    return theta * Math.exp(-theta * t);
                }
            }, new DoubleUnaryOperator() {
                @Override
                public double applyAsDouble(double t) {
                    return theta;
                }
            }, 1.0 / theta, true, theta);
        }

        final Matrix D0 = pair.get(0);
        final Matrix D1 = pair.get(1);
        final double asymptoticRate = rate;
        final DoubleUnaryOperator ccdf = new DoubleUnaryOperator() {
            @Override
            public double applyAsDouble(double t) {
                Matrix pt = new Matrix(1, 1);
                pt.set(0, 0, t);
                return 1.0 - Map_cdf.map_cdf(D0, D1, pt).get(0);
            }
        };
        final DoubleUnaryOperator pdf = new DoubleUnaryOperator() {
            @Override
            public double applyAsDouble(double t) {
                return Map_pdf.map_pdf(D0, D1, t);
            }
        };
        DoubleUnaryOperator hazard = new DoubleUnaryOperator() {
            @Override
            public double applyAsDouble(double t) {
                // h = f/(1-F). Past the point where the ccdf underflows the
                // hazard is the asymptotic decay rate, and returning that is
                // better conditioned than dividing two zeros.
                double c = ccdf.applyAsDouble(t);
                if (c <= 1e-300) {
                    return asymptoticRate > 0 ? asymptoticRate : 0.0;
                }
                return pdf.applyAsDouble(t) / c;
            }
        };
        double mean = rate > 0 ? 1.0 / rate : Double.POSITIVE_INFINITY;
        return new Handles(ccdf, pdf, hazard, mean, isExp, rate);
    }
}
