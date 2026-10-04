/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dph;

import java.util.Arrays;
import jline.lib.butools.dph.CanonicalFromDPH3.DPH3Representation;
import jline.lib.butools.dph.MGFromMoments.MGRepresentation;

public final class DPH3From5Moments {
    private DPH3From5Moments() {}

    /**
     * Returns an order-3 discrete phase-type distribution which has the same 5 moments as given.
     *
     * @param moms The moments to match (length 5)
     * @param prec Numerical precision, default value is 1e-14
     * @return The DPH3Representation containing alpha (initial probability vector) and A (transition probability matrix)
     *
     * Note: Raises an error if the moments are not feasible with a DPH(3).
     *       This procedure first calls 'MGFromMoments', then transforms it to DPH(3) by 'CanonicalFromDPH3'.
     */
    public static DPH3Representation dph3From5Moments(double[] moms, double prec) {
        if (moms.length < 5) {
            throw new IllegalArgumentException("DPH3From5Moments: At least 5 moments are required!");
        }

        // Get MG representation from first 5 moments
        MGRepresentation mgRep = MGFromMoments.mgFromMoments(Arrays.copyOfRange(moms, 0, 5));

        // Convert to canonical DPH(3) form
        return CanonicalFromDPH3.canonicalFromDPH3(mgRep.getAlpha(), mgRep.getA(), prec);
    }

    public static DPH3Representation dph3From5Moments(double[] moms) {
        return dph3From5Moments(moms, 1e-14);
    }
}
