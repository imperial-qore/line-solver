/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.pfqn.nc;

import jline.io.Ret;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;

/**
 * Importance-sampling (IS) estimate of the normalizing constant of a closed
 * LOAD-INDEPENDENT product-form queueing network with M single-server queues of
 * per-class demand L and an aggregated delay of think time Z.
 *
 * <p>This is the load-independent case of {@link Pfqn_ld_is} (capacities
 * mu_i(k) = 1), and the ordinary-network counterpart of the order-independent
 * {@link Pfqn_oi_is} and the pass-and-swap {@link Pfqn_pas_is}: all four are the
 * same sample-an-ordering estimator, differing only in the per-position factor
 * of each station's balance function. For a single-server queue that factor is
 * the demand of the class at that position, L(i,q_p); for the delay it is
 * Z(q_p)/p; for an OI/P&amp;S station it is the reciprocal rank rate
 * 1/mu_i(supp(q_1..q_p)).
 *
 * <p>Writing ell = sum(N), an ordering c of all ell jobs is drawn by placing a
 * uniformly random present class at each step (probability p(c) = product of the
 * reciprocal branching factors), and the sum over ALL ways of cutting c into
 * contiguous per-station segments is computed exactly by dynamic programming:
 * G(N) = E_{C~p}[ S(C)/p(C) ], which is unbiased for the exact constant of
 * {@link Pfqn_nc}.
 *
 * <p>Port of matlab/src/api/pfqn/pfqn_is.m.
 *
 * @see Pfqn_ld_is
 * @see Pfqn_oi_is
 * @see Pfqn_pas_is
 */
public final class Pfqn_is {
    private Pfqn_is() {}

    /**
     * Importance-sampling estimate of the load-independent normalizing constant.
     *
     * @param L       (M x R) per-class service demands at the M single-server queues.
     * @param N       (1 x R) closed population vector, finite.
     * @param Z       (1 x R) aggregated think time (delay) demand; null or zeros if none.
     * @param options solver options; uses options.samples and options.seed.
     * @return G and lG = log(G).
     */
    public static Ret.pfqnNc pfqn_is(Matrix L, Matrix N, Matrix Z, SolverOptions options) {
        return Pfqn_ld_is.pfqn_ld_is(L, N, Z, null, options);
    }
}
