/**
 * @file Markovian Arrival Process maximum operation
 * 
 * Computes MAP representation of maximum inter-arrival times from independent MAP processes.
 * Used for analyzing synchronization and worst-case timing in parallel arrival systems.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Computes the MAP that represents the maximum of two independent MAPs.
 *
 * Given two independent MAPs A and B, this function computes the MAP that represents
 * the maximum of the two inter-arrival times. The resulting MAP describes the process
 * where each arrival corresponds to the maximum of the two independent arrival times.
 *
 * @param A The first MAP stored in a MatrixCell, containing (D0, D1) matrices
 * @param B The second MAP stored in a MatrixCell, containing (D0, D1) matrices
 * @return A MatrixCell representing the MAP of the maximum process
 */
fun map_max(A: MatrixCell, B: MatrixCell): MatrixCell {
    val na = A[0].numRows
    val nb = B[0].numRows
    
    // Compute exit rates: a = -A{1} * ones(na,1) and b = -B{1} * ones(nb,1)
    val onesA = Matrix.ones(na, 1)
    val onesB = Matrix.ones(nb, 1)
    val a = A[0].scale(-1.0).mult(onesA)
    val b = B[0].scale(-1.0).mult(onesB)
    
    // Construct the M0 matrix (generator matrix)
    val M0 = Matrix(na * nb + na + nb, na * nb + na + nb)
    
    // First block row: [krons(A{1},B{1}), kron(a,eye(na)), kron(b,eye(nb))]
    val kronsAB = kroneckerSum(A[0], B[0])
    val eyeNA = Matrix.eye(na)
    val eyeNB = Matrix.eye(nb)
    val kronA = kroneckerProduct(a, eyeNA)
    val kronB = kroneckerProduct(b, eyeNB)
    
    // Set the first block row
    for (i in 0 until na * nb) {
        for (j in 0 until na * nb) {
            M0.set(i, j, kronsAB.get(i, j))
        }
        for (j in 0 until na) {
            M0.set(i, na * nb + j, kronA.get(i, j))
        }
        for (j in 0 until nb) {
            M0.set(i, na * nb + na + j, kronB.get(i, j))
        }
    }
    
    // Second block row: [zeros(nb,na*nb), B{1}, zeros(nb,nb)]
    for (i in 0 until nb) {
        for (j in 0 until na * nb) {
            M0.set(na * nb + i, j, 0.0)
        }
        for (j in 0 until nb) {
            M0.set(na * nb + i, na * nb + j, B[0].get(i, j))
        }
        for (j in 0 until nb) {
            M0.set(na * nb + i, na * nb + na + j, 0.0)
        }
    }
    
    // Third block row: [zeros(na,na*nb), zeros(na,nb), A{1}]
    for (i in 0 until na) {
        for (j in 0 until na * nb) {
            M0.set(na * nb + nb + i, j, 0.0)
        }
        for (j in 0 until nb) {
            M0.set(na * nb + nb + i, na * nb + j, 0.0)
        }
        for (j in 0 until na) {
            M0.set(na * nb + nb + i, na * nb + na + j, A[0].get(i, j))
        }
    }
    
    // Compute the pie vector: [kron(map_pie(A),map_pie(B)), zeros(1,na), zeros(1,nb)]
    val pieA = map_pie(A)
    val pieB = map_pie(B)
    val pieKron = kroneckerProduct(pieA, pieB)
    val pie = Matrix(1, na * nb + na + nb)
    
    // Set the pie vector
    for (j in 0 until na * nb) {
        pie.set(0, j, pieKron.get(0, j))
    }
    for (j in 0 until na) {
        pie.set(0, na * nb + j, 0.0)
    }
    for (j in 0 until nb) {
        pie.set(0, na * nb + na + j, 0.0)
    }
    
    // Compute M1 = pie' * ones(1, size(M0,1))
    val pieTrans = pie.transpose()
    val onesRow = Matrix.ones(1, M0.numRows)
    val M1 = pieTrans.mult(onesRow)
    
    // Return the MAP
    val result = MatrixCell()
    result[0] = M0
    result[1] = M1
    return result
}

/**
 * Computes the Kronecker sum of two matrices: A ⊕ B = A \otimes I + I \otimes B
 */
private fun kroneckerSum(A: Matrix, B: Matrix): Matrix {
    val eyeA = Matrix.eye(A.numRows)
    val eyeB = Matrix.eye(B.numRows)
    val kronA = kroneckerProduct(A, eyeB)
    val kronB = kroneckerProduct(eyeA, B)
    return kronA.add(1.0, kronB)
}

/**
 * Computes the Kronecker product of two matrices.
 */
private fun kroneckerProduct(A: Matrix, B: Matrix): Matrix {
    val result = Matrix(A.numRows * B.numRows, A.numCols * B.numCols)
    
    for (i in 0 until A.numRows) {
        for (j in 0 until A.numCols) {
            val aij = A.get(i, j)
            for (k in 0 until B.numRows) {
                for (l in 0 until B.numCols) {
                    result.set(i * B.numRows + k, j * B.numCols + l, aij * B.get(k, l))
                }
            }
        }
    }
    
    return result
}
/**
 * MAP max algorithms
 */
@Suppress("unused")
class MapMaxAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}