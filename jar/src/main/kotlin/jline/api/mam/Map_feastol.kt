/**
 * @file Markovian Arrival Process feasibility tolerance constants
 * 
 * Provides standard tolerance values for numerical feasibility checks in MAP algorithms.
 * Ensures consistent numerical accuracy across different MAP validation procedures.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

/**
 * Returns the feasibility tolerance for MAPs.
 *
 *
 * This method provides a standard tolerance level used in feasibility checks for MAPs.
 * It might be used to determine the numerical accuracy required in calculations involving MAPs.
 *
 * @return the feasibility tolerance value
 */
fun map_feastol(): Int {
    return 8
}
/**
 * MAP feastol algorithms
 */
@Suppress("unused")
class MapFeastolAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}