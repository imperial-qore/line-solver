/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dph

/**
 * Returns an order-2 discrete phase-type distribution which has the same 3 moments as given.
 *
 * @param moms The moments to match (length 3)
 * @return The DPH2Representation containing alpha (initial probability vector) and A (transition probability matrix)
 *
 * Note: Raises an error if the moments are not feasible with a DPH(2).
 *       This procedure first calls 'MGFromMoments', then transforms it to DPH(2) by 'CanonicalFromDPH2'.
 */
fun dph2From3Moments(moms: DoubleArray): DPH2Representation {
    require(moms.size >= 3) { "DPH2From3Moments: At least 3 moments are required!" }

    // Get MG representation from first 3 moments
    val mgRep = mgFromMoments(moms.sliceArray(0..2))

    // Convert to canonical DPH(2) form
    return canonicalFromDPH2(mgRep.alpha, mgRep.A)
}
