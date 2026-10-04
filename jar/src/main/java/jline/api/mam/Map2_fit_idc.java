/**
 * @file MAP(2) fitting from three moments and an index of dispersion
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.GlobalConstants;
import jline.io.Ret;

/**
 * Fits a second-order MAP matching the first three moments and the asymptotic index of
 * dispersion.
 *
 * A MAP(2) has a geometrically decaying autocorrelation, so its index of dispersion obeys
 * I = SCV + (SCV-1)*g2/(1-g2), as reported in Section 5.2.2 of Casale, Mi, Cherkasova and
 * Smirni, IEEE Trans. Soft. Eng. 37(5), 2011. The relation is inverted in closed form as
 * g2 = (I-SCV)/(I-1) and the decay rate is passed to map2_fit. Burstiness cannot be
 * represented when SCV &lt;= 1 or when I &lt; SCV, and the paper falls back to an
 * exponential in both cases. A third moment outside the feasible region is replaced by its
 * lower limit (3/2)*e2^2/e1, the largest heavy-tail decay a MAP(2) admits.
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
public final class Map2_fit_idc {
    private Map2_fit_idc() {}

    /** Exact fit of the four descriptors. */
    public static final int STATUS_EXACT = 0;
    /** Burstiness not representable, an exponential was returned. */
    public static final int STATUS_EXPONENTIAL = 1;
    /** Third moment clamped to its feasible lower limit. */
    public static final int STATUS_E3_CLAMPED = 2;
    /** Third moment selected automatically by map2_fit. */
    public static final int STATUS_E3_AUTO = 3;
    /** Fit failed, an exponential was returned. */
    public static final int STATUS_FAILED = 4;

    /**
     * Fits a MAP(2) to three moments and an index of dispersion.
     *
     * @param e1  mean inter-arrival time
     * @param e2  second moment of the inter-arrival times
     * @param e3  third moment of the inter-arrival times
     * @param idc asymptotic index of dispersion
     * @return the fitted MAP and the fallback taken, if any
     */
    public static Ret.mamMAPFitIdcReturn map2_fit_idc(double e1, double e2, double e3, double idc) {
        double scv = (e2 - e1 * e1) / (e1 * e1);
        Ret.mamMAPFitIdcReturn out = new Ret.mamMAPFitIdcReturn();

        if (scv <= 1 + GlobalConstants.FineTol || idc < scv) {
            out.MAP = Map_exponential.map_exponential(e1);
            out.status = STATUS_EXPONENTIAL;
            return out;
        }

        double g2 = (idc - scv) / (idc - 1);

        Ret.mamMAPFitReturn fit = Map2_fit.map2_fit(e1, e2, e3, g2);
        if (fit.error == 0) {
            out.MAP = fit.MAP;
            out.status = STATUS_EXACT;
            return out;
        }

        double e3min = (3.0 / 2 + 1e-6) * e2 * e2 / e1;
        if (e3 < e3min) {
            fit = Map2_fit.map2_fit(e1, e2, e3min, g2);
            if (fit.error == 0) {
                out.MAP = fit.MAP;
                out.status = STATUS_E3_CLAMPED;
                return out;
            }
        }

        fit = Map2_fit.map2_fit(e1, e2, -1, g2);
        if (fit.error == 0) {
            out.MAP = fit.MAP;
            out.status = STATUS_E3_AUTO;
            return out;
        }

        out.MAP = Map_exponential.map_exponential(e1);
        out.status = STATUS_FAILED;
        return out;
    }
}
