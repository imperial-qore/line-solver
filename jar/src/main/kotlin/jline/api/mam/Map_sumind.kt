/**
 * @file Markovian Arrival Process independent summation operations
 * 
 * Computes MAP representations of sums of independent MAP processes for modeling
 * parallel arrival streams. Used for load aggregation and multi-source system analysis.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Computes the Markovian Arrival Process (MAP) representing the sum of `n` independent MAPs.
 *
 *
 * This method constructs a new MAP with transition matrices `D0` and `D1` by summing `n` independent MAPs,
 * specified by the input array `MAPs`. The result is a larger MAP with the combined state space of the original MAPs,
 * where transitions are structured to account for the independent nature of the original MAPs.
 *
 * @param MAPs an array of MatrixCells, each containing the transition matrices `D0` and `D1` of an independent MAP
 * @return a MatrixCell containing the transition matrices `D0` and `D1` of the resulting summed MAP
 */
fun map_sumind(MAPs: Array<MatrixCell>): MatrixCell {
    val n = MAPs.size
    val order = Matrix(1, n)
    for (i in 0..<n) {
        order[i] = MAPs[i][0].numRows.toDouble()
    }
    val D0 = Matrix.zeros(order.elementSum().toInt(), order.elementSum().toInt())
    val D1 = Matrix.zeros(order.elementSum().toInt(), order.elementSum().toInt())
    var curpos = 0
    for (i in 0..<n) {
        D0.setSliceEq(curpos, curpos + order[i].toInt(), curpos, curpos + order[i].toInt(), MAPs[i][0])
        if (i < (n - 1)) {
            D0.setSliceEq(curpos,
                curpos + order[i].toInt(),
                curpos + order[i].toInt(),
                curpos + order[i].toInt() + order[i + 1].toInt(),
                MAPs[i][1].mult(Matrix.ones(order[i].toInt(), 1)).mult(map_pie(MAPs[i + 1])))
        } else {
            D1.setSliceEq(curpos,
                curpos + order[i].toInt(),
                0,
                order[i].toInt(),
                MAPs[i][1].mult(Matrix.ones(order[i].toInt(), 1)).mult(map_pie(MAPs[i])))
        }
        curpos += order[i].toInt()
    }
    return MatrixCell(D0, D1)
}


//        Matrix D0 = new Matrix("[-2,1;0,-3]");
//        Matrix D1 = new Matrix("[1,0;2,1]");
//        Matrix C0 = new Matrix("[-11,11,0;0,-22,22;0,0,-33]");
//        Matrix C1 = new Matrix("[0,0,0;0,0,0;33,0,0]");
//        Matrix B0 = new Matrix("[-111,111,0;0,-222,111;0,0,-333]");
//        Matrix B1 = new Matrix("[0,0,0;111,0,0;333,0,0]");
//        MatrixCell D = new MatrixCell(D0,D1);
//        MatrixCell C = new MatrixCell(C0,C1);
//        MatrixCell B = new MatrixCell(B0,B1);
//        MatrixCell[] MAPs = new MatrixCell[3];
//        MAPs[0]=D;
//        MAPs[1]=B;
//        MAPs[2]=C;
//        MatrixCell SMAP = map_sumind(MAPs);
//        SMAP.get(0).print();
//        SMAP.get(1).print();
//    }
/**
 * MAP sumind algorithms
 */
@Suppress("unused")
class MapSumindAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}