/**
 * @file Phase-type distribution reindexing utilities
 * 
 * Reindexes phase-type distribution maps for network models using integer station and class indices.
 * Essential utility for integrating MAM algorithms with LINE network modeling infrastructure.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.lang.JobClass
import jline.lang.NetworkStruct
import jline.lang.nodes.Station
import jline.util.matrix.MatrixCell

/**
 * Reindexes phase-type (PH) distributions for a network model based on station and job class indices.
 *
 * This method converts PH distributions from Station/JobClass object indexing to integer-based indexing.
 * The new index system corresponds to the internal indexing used in a NetworkStruct object,
 * which describes the structure of the network, including the number of stations and classes.
 *
 * @param sn the NetworkStruct object containing the network structure, indexed PH distributions (sn.proc), and internal indexing information
 * @return a reindexed map where the key is an integer station index and the value is another map with integer job class indices and MatrixCell values
 */
fun ph_reindex(sn: NetworkStruct): Map<Int, MutableMap<Int, MatrixCell?>> {
    val result: MutableMap<Int, MutableMap<Int, MatrixCell?>> = HashMap()

    for (i in 0..<sn.nstations) {
        for (j in 0..<sn.nclasses) {
            if (j == 0) {
                result[i] = HashMap()
            }
            result[i]!![j] = sn.proc[sn.stations[i]]!![sn.jobclasses[j]]
        }
    }
    return result
}
/**
 * Ph Reindex algorithms
 */
@Suppress("unused")
class PhReindex {
    companion object {
        // Class documentation marker for Dokka
    }
}