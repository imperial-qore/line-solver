/**
 * @file Marked Markovian Arrival Process superposition operations
 * 
 * Combines multiple MMAP processes into superposed multiclass arrival streams.
 * Fundamental for modeling aggregated traffic sources in complex queueing networks.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Combines two MMAPs into one superposed MMAP.
 *
 *
 * This method creates a superposed MMAP from two given MMAPs (MMAPa and MMAPb). Depending on the option specified:
 * - "default": Treats the MMAPs as unrelated, and each class in MMAPa and MMAPb is considered distinct in the resulting superposed MMAP.
 * - "match": Assumes that MMAPa and MMAPb have the same number of classes and maps corresponding classes from both into the same class in the superposed MMAP.
 *
 * @param MMAPa the first MMAP
 * @param MMAPb the second MMAP
 * @param opt   "default" for unrelated MMAPs, "match" for matching classes
 * @return a new MMAP representing the superposition of the input MMAPs
 */
fun mmap_super(MMAPa: MatrixCell, MMAPb: MatrixCell, opt: String): MatrixCell? {
    val sup = MatrixCell()
    if (opt == "default") {
        val K1 = MMAPa.size() - 2
        val K2 = MMAPb.size() - 2
        val n1 = MMAPa[0].length()
        val n2 = MMAPb[0].length()

        sup[0] = MMAPa[0].krons(MMAPb[0])
        sup[1] = MMAPa[1].krons(MMAPb[1])

        for (i in 0..<K1) {
            sup[2 + i] = MMAPa[2 + i].krons(Matrix(n2, n2, 0))
        }

        for (j in 0..<K2) {
            val a = Matrix(n1, n1, 0)
            sup[2 + K1 + j] = a.krons(MMAPb[2 + j])
        }
    } else if (opt == "match") {
        val K1 = MMAPa.size()
        val K2 = MMAPb.size()

        if (K1 != K2) {
            throw RuntimeException("class matching failed: MMAPs have different number of classes")
        }

        for (i in 0..<K1) {
            sup[i] = MMAPa[i].krons(MMAPa[i])
        }
    } else {
        throw RuntimeException("unrecognized option")
    }

    return mmap_normalize(sup)
}

/**
 * Combines two MMAPs into one superposed MMAP using the default option.
 *
 * @param MMAPa the first MMAP
 * @param MMAPb the second MMAP
 * @return a new MMAP representing the superposition of the input MMAPs
 */

fun mmap_super(MMAPa: MatrixCell, MMAPb: MatrixCell): MatrixCell? {
    return mmap_super(MMAPa, MMAPb, "default")
}

/**
 * Combines a list of MMAPs into one superposed MMAP.
 *
 *
 * This method removes null entries from the list and then superposes the remaining MMAPs.
 *
 * @param MMAPa the list of MMAPs to be combined
 * @return a new MMAP representing the superposition of the input MMAPs
 */

fun mmap_super(MMAPa: MatrixCell): MatrixCell? {
    MMAPa.removeNull()
    var SUP: MatrixCell? = MatrixCell()
    SUP!![0] = MMAPa[0]
    for (i in 1..<MMAPa.size()) {
        val MMAPb = MatrixCell()
        MMAPb[0] = MMAPa[i]
        SUP = mmap_super(SUP!!, MMAPb)
    }
    return SUP
}
/**
 * MMAP super algorithms
 */
@Suppress("unused")
class MmapSuperAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}