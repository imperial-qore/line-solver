/**
 * @file Markovian Arrival Process superposition operations
 * 
 * Creates superposition of MAP processes using Kronecker product techniques.
 * Fundamental for modeling independent arrival stream combinations in queueing networks.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.MatrixCell

/**
 * Creates a superposition of two Markovian Arrival Processes (MAPs) to form a new MAP.
 *
 *
 * This function takes two input MAPs and constructs a superposed MAP by performing the Kronecker product
 * of their corresponding matrices. The resulting MAP represents the combined behavior of the two input processes.
 * After constructing the superposed MAP, it is normalized to ensure it represents a valid stochastic process.
 *
 *
 * @param MAPa The first Markovian Arrival Process stored in a MatrixCell.
 * @param MAPb The second Markovian Arrival Process stored in a MatrixCell.
 * @return A MatrixCell representing the superposed MAP, formed by combining the input MAPs.
 */
fun map_super(MAPa: MatrixCell, MAPb: MatrixCell): MatrixCell {
    val sup = MatrixCell()
    val sizeA = MAPa.size() - 2
    val sizeB = MAPb.size() - 2
    val lengthA = MAPa[0].length()
    val lengthB = MAPb[0].length()

    sup[0] = MAPa[0].krons(MAPb[0])
    sup[1] = MAPa[1].krons(MAPb[1])

    return map_normalize(sup)
}
/**
 * MAP super algorithms
 */
@Suppress("unused")
class MapSuperAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}