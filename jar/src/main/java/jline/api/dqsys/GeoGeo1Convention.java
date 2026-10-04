/**
 * @file Slot-boundary convention for discrete-time Geo/Geo/1 analysis
 *
 * @since LINE 3.1.0
 */
package jline.api.dqsys;

/**
 * Observation epoch for a discrete-time Geo/Geo/1 queue.
 *
 * <p>The two entries are not two different systems. They are the same system --
 * Daduna's LA-rule (all arrivals and departures occur at the end of their slot)
 * combined with the D/A-rule (when both fall on one epoch, the departure is
 * resolved first) -- observed at two different instants within the slot
 * boundary. The queue-length process is
 * {@code X(t+1) = X(t) - D(t) + A(t)}; {@link #LAS_DA} is the law of
 * {@code X}, taken after both events, and {@link #EAS} is the law of the
 * intermediate {@code Y(t) = X(t) - D(t)}, taken after the departure and before
 * the arrival.
 *
 * <p>Since the two epochs are one departure apart and the departure rate equals
 * the arrival rate {@code a} in steady state, the means differ by exactly
 * {@code a} in population and, by Little's law, by exactly one slot in sojourn
 * time. The queueing delay is the same under both.
 *
 * <p>Note that the regulation rule itself is not what separates them. Daduna
 * (2001), section 3.1, records that with ample waiting room the D/A-rule and the
 * A/D-rule yield the same stationary boundary distribution, because a service
 * occupying at least one whole slot prevents a job from departing in its own
 * arrival slot. That is also why the LDES engine's intra-slot event ordering
 * does not move the estimates; it fixes a convention that would otherwise be
 * emergent rather than modelled.
 *
 * <p><b>Why two constants suffice.</b> The literature also names a late-arrival
 * system with immediate access (LAS-IA), which raises the question of whether a
 * third constant is missing. It is not, because the access rule and the
 * observation epoch are not separately identifiable. Delaying by one slot the
 * first slot in which an arriving job may be served shifts every mean the same
 * way that stepping the observation epoch back by one does: verified by
 * slot-level simulation, the delayed-access system's mid-slot population equals
 * the immediate-access system's boundary population (0.69897 against 0.69897 at
 * a=0.3, s=0.6; 0.53272 against 0.53338 at a=0.2, s=0.5, within sampling
 * noise). So the reachable means form the one-parameter family
 * {@code E[T] = (1-a)/(s-a) + m}, {@code E[N] = a E[T]}, indexed by an integer
 * slot offset m, and LAS-IA and EAS name the same member of it. The two
 * constants below are {@code m = 0} and {@code m = -1}, the only two that carry
 * standard names; any further offset re-indexes the same family rather than
 * describing a new system.
 *
 * <p>Reference: H. Daduna, Queueing Networks with Discrete Time Scale, LNCS
 * 2046, Springer 2001, section 2.1 (figures 2.2 and 2.3, theorem 2.3,
 * corollary 2.7).
 */
public enum GeoGeo1Convention {

    /**
     * Late arrival system with delayed access: the state at the slot boundary,
     * after both the departure and the arrival of the slot. A job never departs
     * in its own arrival slot and the service time is geometric on
     * {@code {1,2,...}} with mean {@code 1/s}.
     *
     * <p>This is the convention realized by the LDES engine when a Source with a
     * {@code Geometric(a)} interarrival time feeds an FCFS Queue with a
     * {@code Geometric(s)} service time, because both are supported on
     * {@code {1,2,...}}.
     *
     * <p>Matches Daduna corollary 2.7 term for term, including the tail ratio
     * {@code bq/(cp)} in his notation ({@code b = a}, {@code p = s}).
     */
    LAS_DA,

    /**
     * Early arrival system: the intermediate state, after the departure of the
     * slot and before its arrival. Equivalently, the sojourn is measured net of
     * the arrival slot, so the effective service time is geometric on
     * {@code {0,1,...}} with mean {@code (1-s)/s}.
     *
     * <p>This is the process observed at {@code t+1/2} in Daduna's figure 2.3.
     * Its empty probability is {@code 1-r}, NOT {@code 1-rho}: a job can arrive
     * and be cleared within one slot, so the empty-at-epoch probability and the
     * busy-slot fraction come apart even though the utilization is still
     * {@code rho} in both.
     */
    EAS
}
