/**
 * @file Quasi-Birth-Death process BMAP/BMAP/1 queue analysis
 *
 * Constructs QBD transition blocks for batch arrival and service systems.
 * Handles batch arrivals with specified probability distribution.
 *
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Result containing QBD matrices for BMAP/BMAP/1 queue
 *
 * @param A0 Local transition block (Kronecker sum of arrival D0 and service D0)
 * @param A_1 Downward transition block (service completions)
 * @param A1 List of upward transition blocks for each batch size
 * @param B0 Boundary local block (Kronecker sum of arrival D0 and identity)
 * @param B1 List of boundary upward blocks for each batch size
 */
data class QbdBmapResult(
    val A0: Matrix,
    val A_1: Matrix,
    val A1: List<Matrix>,
    val B0: Matrix,
    val B1: List<Matrix>
)

/**
 * Set up QBD matrices for BMAP/BMAP/1 queue analysis
 *
 * Constructs QBD transition blocks matching MATLAB qbd_bmapbmap1.m:
 *   A1{b} = kron(MAPa{2}*pbatcha(b), eye(ns))
 *   A0 = krons(MAPa{1}, MAPs{1})
 *   A_1 = kron(eye(na), MAPs{2})
 *   B0 = krons(MAPa{1}, eye(ns))
 *   B1{b} = kron(MAPa{2}*pbatcha(b), eye(ns))
 *
 * @param MAPa Arrival MAP (Batch Markovian Arrival Process)
 * @param pbatcha Batch probabilities for arrivals
 * @param MAPs Service MAP (Batch Markovian Service Process)
 * @return QbdBmapResult containing the QBD matrices
 */
fun qbd_bmapbmap1(MAPa: MatrixCell, pbatcha: Matrix, MAPs: MatrixCell): QbdBmapResult {
    val na = MAPa[0].numRows
    val ns = MAPs[0].numRows
    val maxbatch = pbatcha.length()

    val IA = Matrix.eye(na)
    val IS = Matrix.eye(ns)

    // Upward transition blocks for each batch size
    // A1{b} = kron(MAPa{2}*pbatcha(b), eye(ns))
    val A1 = ArrayList<Matrix>()
    for (b in 0 until maxbatch) {
        val A1_b = MAPa[1].scale(pbatcha[b]).kron(IS)
        A1.add(A1_b)
    }

    // Local transitions: A0 = krons(MAPa{1}, MAPs{1}) = kron(D0_arr, I_ns) + kron(I_na, D0_srv)
    val A0 = MAPa[0].kron(IS).add(IA.kron(MAPs[0]))

    // Downward transitions: A_1 = kron(eye(na), MAPs{2})
    val A_1 = IA.kron(MAPs[1])

    // Boundary local: B0 = krons(MAPa{1}, eye(ns)) = kron(D0_arr, I_ns) + kron(I_na, I_ns)
    val B0 = MAPa[0].kron(IS).add(Matrix.eye(na * ns))

    // Boundary upward blocks: B1{b} = kron(MAPa{2}*pbatcha(b), eye(ns))
    val B1 = ArrayList<Matrix>()
    for (b in 0 until maxbatch) {
        val B1_b = MAPa[1].scale(pbatcha[b]).kron(IS)
        B1.add(B1_b)
    }

    return QbdBmapResult(A0, A_1, A1, B0, B1)
}

/**
 * QBD bmapbmap1 algorithms
 */
@Suppress("unused")
class QbdBmapbmap1Algo {
    companion object {
        // Class documentation marker for Dokka
    }
}
