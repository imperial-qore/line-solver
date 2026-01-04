package jline.api.mc

import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath
import java.util.*
import java.util.stream.Collectors
import java.util.stream.IntStream

/**
 * Returns the stochastic complement of a DTMC
 *
 * @param P Transition matrix of the DTMC
 * @param I Indexes of states to be kept in the stochastic complement
 * @return Transition matrix of the stochastic complement
 * @throws RuntimeException If the transition matrix is null
 */
fun dtmc_stochcomp(P: Matrix, I: MutableList<Int?>?): Matrix {
    //Note that in this function, List is used instead of Matrix for performance consideration
    var I = I

    val lengthP = FastMath.max(P.numCols, P.numRows)
    if (I == null || I.size == 0) {
        I = ArrayList()
        for (i in 0..<FastMath.ceil(lengthP / 2.0).toInt()) I.add(i)
    }

    val Ic = IntStream.rangeClosed(0, lengthP - 1).boxed().collect(Collectors.toList())
    Ic.removeAll(I)

    val P11 = Matrix(I.size, I.size)
    val P12 = Matrix(I.size, Ic.size)
    val P21 = Matrix(Ic.size, I.size)
    val P22 = Matrix(Ic.size, Ic.size)

    for (colIdx in 0..<P.numCols) {
        for (rowIdx in 0..<P.numRows) {
            val value = P[rowIdx, colIdx]
            if (value > 0.0 || java.lang.Double.isNaN(value)) {
                if (I.contains(colIdx)) {
                    if (I.contains(rowIdx)) {
                        P11[I.indexOf(rowIdx), I.indexOf(colIdx)] = value
                    } else {
                        P21[Ic.indexOf(rowIdx), I.indexOf(colIdx)] = value
                    }
                } else {
                    if (I.contains(rowIdx)) {
                        P12[I.indexOf(rowIdx), Ic.indexOf(colIdx)] = value
                    } else {
                        P22[Ic.indexOf(rowIdx), Ic.indexOf(colIdx)] = value
                    }
                }
            }
        }
    }

    val values = DoubleArray(Ic.size)
    Arrays.fill(values, 1.0)
    val S2 = Matrix.diagMatrix(null, values, 0, values.size).sub(P22)

    // S=P11+P12*(S2 \ P21);
    val s2_p21 = Matrix(S2.numRows, P21.numCols)
    
    // Check if S2 is singular and handle it
    if (S2.numRows == 0 || S2.numCols == 0) {
        // If S2 is empty, return P11 directly
        return P11
    }
    
    Matrix.solve(S2, P21, s2_p21)
    val S = P11.add(1.0, P12.mult(s2_p21, null))
    return S
}
/**
 * DTMC stochcomp algorithms
 */
@Suppress("unused")
class DtmcStochcompAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}