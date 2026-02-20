/**
 * @file Markovian Arrival Process large map threshold constants
 * 
 * Provides size thresholds for determining when MAP algorithms should switch to
 * optimized implementations for large-scale problems. Used for performance optimization.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

/**
 * Returns the maximum size considered for large MAPs.
 *
 *
 * This method provides a standard size threshold for determining what constitutes a large MAP.
 * This threshold might be used in performance optimizations or to decide when to switch to more efficient algorithms.
 *
 * @return the size threshold for large MAPs
 */
fun map_largemap(): Int {
    return 100
}
/**
 * MAP largemap algorithms
 */
@Suppress("unused")
class MapLargemapAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}