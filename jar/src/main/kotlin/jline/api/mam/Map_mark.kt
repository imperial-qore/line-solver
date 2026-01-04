/**
 * @file Markovian Arrival Process marking for multiclass processes
 * 
 * Creates Marked MAP (MMAP) representations by adding class labels to MAP arrivals.
 * Essential for modeling multiclass arrival streams in queueing networks and service systems.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.io.line_warning
import jline.io.mfilename
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Creates a Marked Markovian Arrival Process (MMAP) by marking a given Markovian Arrival Process (MAP)
 * with additional phases based on specified marking probabilities.
 *
 *
 * This function constructs a new MMAP by extending the original MAP with additional phases, where each phase corresponds to
 * a marking according to the provided probability matrix. If the marking probabilities do not sum to 1, they are normalized.
 * The original MAP matrices are retained, and new matrices corresponding to the additional marked phases are generated and scaled
 * according to the specified probabilities.
 *
 *
 * @param MAP  The original Markovian Arrival Process stored in a MatrixCell.
 * @param prob A matrix containing the marking probabilities. Each column represents a different marking, and the sum of the
 * elements is expected to be 1. If not, the probabilities are normalized.
 * @return A MatrixCell representing the Marked Markovian Arrival Process (MMAP), which includes the original MAP matrices
 * and the additional matrices for the marked phases.
 */
fun map_mark(MAP: MatrixCell, prob: Matrix): MatrixCell {
    val prob_local = prob.copy()
    if ((prob_local.elementSum() < 1 - 1e-6) || (prob_local.elementSum() > 1 + 1e-6))
        line_warning(mfilename(object : Any() {}),
        "Input marking probabilities do not sum to 1. Normalizing.")
    prob_local.scaleEq(prob_local.elementSum())
    val R = prob_local.numCols
    val mmap = MatrixCell(2 + R)
    mmap[0] = MAP[0].copy()
    mmap[1] = MAP[1].copy()
    for (r in 0..<R) {
        mmap[2 + r] = Matrix.createLike(MAP[1])
        val a = MAP[1].copy()
        a.scaleEq(prob_local[r])
        mmap[2 + r] = a
    }
    return mmap
}
/**
 * MAP mark algorithms
 */
@Suppress("unused")
class MapMarkAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}