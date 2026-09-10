/**
 * Markovian Arrival Process feasibility tolerance constants.
 *
 * Provides standard tolerance values for numerical feasibility checks in MAP algorithms.
 * Ensures consistent numerical accuracy across different MAP validation procedures.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

public final class Map_feastol {
    private Map_feastol() {}

    /**
     * Returns the feasibility tolerance EXPONENT for MAPs.
     *
     * This is the exponent k of the toolbox feasibility tolerance 10^-k, so the
     * tolerance itself is 10^-8. Callers either raise it, as
     * Math.pow(10.0, -map_feastol()), or compare a tolerance magnitude against it
     * directly, as in {@link Map_isfeasible}. It is NOT the tolerance:
     * map_feastol() returns 8, not 1e-8.
     *
     * @return the exponent k of the feasibility tolerance 10^-k
     */
    public static int map_feastol() {
        return 8;
    }
}
