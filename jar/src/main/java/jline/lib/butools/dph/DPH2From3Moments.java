/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dph;

import java.util.Arrays;
import jline.lib.butools.dph.CanonicalFromDPH2.DPH2Representation;
import jline.lib.butools.dph.MGFromMoments.MGRepresentation;

public final class DPH2From3Moments {
    private DPH2From3Moments() {}

    /**
     * Returns an order-2 discrete phase-type distribution which has the same 3 moments as given.
     *
     * @param moms The moments to match (length 3)
     * @return The DPH2Representation containing alpha (initial probability vector) and A (transition probability matrix)
     *
     * Note: Raises an error if the moments are not feasible with a DPH(2).
     *       This procedure first calls 'MGFromMoments', then transforms it to DPH(2) by 'CanonicalFromDPH2'.
     */
    public static DPH2Representation dph2From3Moments(double[] moms) {
        if (moms.length < 3) {
            throw new IllegalArgumentException("DPH2From3Moments: At least 3 moments are required!");
        }

        // Get MG representation from first 3 moments
        MGRepresentation mgRep = MGFromMoments.mgFromMoments(Arrays.copyOfRange(moms, 0, 3));

        // Convert to canonical DPH(2) form
        return CanonicalFromDPH2.canonicalFromDPH2(mgRep.getAlpha(), mgRep.getA());
    }
}
