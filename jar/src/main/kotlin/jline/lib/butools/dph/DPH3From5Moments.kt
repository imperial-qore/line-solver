/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dph

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
fun dph3From5Moments(moms: DoubleArray, prec: Double = 1e-14): DPH3Representation {
    require(moms.size >= 5) { "DPH3From5Moments: At least 5 moments are required!" }

    // Get MG representation from first 5 moments
    val mgRep = mgFromMoments(moms.sliceArray(0..4))

    // Convert to canonical DPH(3) form
    return canonicalFromDPH3(mgRep.alpha, mgRep.A, prec)
}
