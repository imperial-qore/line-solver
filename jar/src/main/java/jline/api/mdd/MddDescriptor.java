package jline.api.mdd;

import java.util.List;

/**
 * Kronecker rate descriptor of a structured model, the input of {@link Mdd_mcd}.
 *
 * <p>Built by {@link Mdd_descriptor} (count-plus-in-service-phase local states,
 * non-preemptive) or {@link Mdd_ps} (per-phase-count local states, shared
 * servers).</p>
 */
public class MddDescriptor {

    /** Number of levels, i.e. stations. */
    public int K;
    /** Closed population; the conservation law the level marginals must satisfy. */
    public int N;
    /** Local domain per level. */
    public int[] domain;
    /** Station service rates, 1/E[S]. */
    public double[] mu;
    /** Servers per station; Double.POSITIVE_INFINITY for a delay station. */
    public double[] servers;
    /** Station-to-station routing matrix. */
    public double[][] P;
    /** Phases per station, 1 when exponential. */
    public int[] nphases;
    /**
     * valuemap[i][idx] is the physical occupancy of station i in local state idx.
     *
     * <p>A level whose local state encodes more than a count (a station holding
     * both a population and a service phase) needs this map; without one the
     * index would be the quantity.</p>
     */
    public double[][] valuemap;
    /** Initial local index per station. */
    public int[] init;
    /** Successor function over local indices, for {@link Mdd_reachset}. */
    public MddNextState nextfun;
    /** The events of the descriptor. */
    public List<MddEvent> events;
    /**
     * Optional conservation law as weights' * QLen = value, overriding the
     * closed-population test. Null when the population N is the invariant.
     */
    public double[] invariantWeights;
    /** Value of the invariant when invariantWeights is set. */
    public double invariantValue;
}
