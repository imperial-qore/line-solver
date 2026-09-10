/**
 * @file Signals that the SJN waiting time equation has no solution at a population
 *
 * SJN starves long jobs as the station saturates and the arrival theorem then fails: once the work
 * brought by the jobs no longer than the tagged one saturates the server, the conditional waiting
 * time equation has no solution. Capping the utilization is not an option in the population
 * lattice, the cap rescaling the profile that the next population step reads back, so the
 * correction compounds and the recursion oscillates. Pfqn_amvasjn solves the same profile by a
 * fixed point and does cap, the iteration being self-consistent.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.mva;

public class SjnStarvationException extends RuntimeException {
    private static final long serialVersionUID = 1L;

    public SjnStarvationException(String message) {
        super(message);
    }
}
