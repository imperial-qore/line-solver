/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.pfqn;

/**
 * Sojourn-time moments at the processor-sharing station of the closed
 * terminal-driven system of Mitra and Morrison (1983).
 *
 * <p>Port of MATLAB pfqn_respt_ps_moments.m.</p>
 */
public class PfqnResptPsResult {
    /** Per-class mean sojourn times; NaN where no route applies. */
    public final double[] W;
    /** Per-class second moments of the sojourn time; NaN where no route applies. */
    public final double[] W2;
    /** Per-class route taken: "exact", "asymptotic", "unavailable" or "none". */
    public final String[] method;
    /** Leading coefficient of E[W^2] ~ c0 + c1/expansionParam; NaN on the exact route. */
    public final double[] c0;
    /** First correction of E[W^2] ~ c0 + c1/expansionParam; NaN on the exact route. */
    public final double[] c1;
    /** Per-class unutilized fraction 1 - sum_r lambda_r/q_r of the CPU in the open counterpart. */
    public final double[] alpha;
    /** Size of the exact state space that the tagged class would need. */
    public final double[] nstates;
    /** The large parameter Nexp = max_r Z(r)/S(r). */
    public final double expansionParam;

    public PfqnResptPsResult(double[] W, double[] W2, String[] method, double[] c0,
                             double[] c1, double[] alpha, double[] nstates,
                             double expansionParam) {
        this.W = W;
        this.W2 = W2;
        this.method = method;
        this.c0 = c0;
        this.c1 = c1;
        this.alpha = alpha;
        this.nstates = nstates;
        this.expansionParam = expansionParam;
    }
}
