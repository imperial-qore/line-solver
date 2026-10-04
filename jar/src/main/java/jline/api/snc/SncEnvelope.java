/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.snc;

/**
 * An MGF envelope evaluated at a Chernoff parameter.
 *
 * <p>The whole family passes envelopes as functions of theta rather than as
 * numbers, because every bound is an infimum over theta and a composition
 * (superposition, leftover service, concatenation, departure) has to be
 * re-evaluated at whatever theta the search asks for. Building the numbers
 * eagerly at one theta would defeat the optimization.</p>
 */
public interface SncEnvelope {
    /**
     * @param theta Chernoff parameter, theta &gt; 0
     * @return the pair {sigma(theta), rho(theta)}; a non-finite entry marks an
     *         infeasible theta and is discarded by {@link Snc_thetaopt}
     */
    double[] eval(double theta);
}
