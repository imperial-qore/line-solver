package jline.lib.smc

import jline.util.matrix.Matrix
import kotlin.math.abs
import kotlin.math.max
import kotlin.math.min

/**
 * QBD_RAP_ParsePara checks the validity of the input matrices A0, A1 and A2
 * for a QBD with RAP components. The evaluated conditions are necessary but not
 * sufficient.
 */
fun QBD_RAP_ParsePara(A0: Matrix, A1: Matrix, A2: Matrix) {
    // Check dimensions
    if (A0.numRows != A0.numCols) {
        throw IllegalArgumentException("A0 is not a square matrix")
    }
    if (A1.numRows != A1.numCols) {
        throw IllegalArgumentException("A1 is not a square matrix")
    }
    if (A2.numRows != A2.numCols) {
        throw IllegalArgumentException("A2 is not a square matrix")
    }
    if (A0.numRows != A1.numRows) {
        throw IllegalArgumentException("The matrices A0 and A1 do not have the same dimension")
    }
    if (A0.numRows != A2.numRows) {
        throw IllegalArgumentException("The matrices A0 and A2 do not have the same dimension")
    }
    
    // Check zero row sum
    val sum_matrix = A0.add(1.0, A1).add(1.0, A2)
    val row_sums = sum_matrix.sumRows()
    val max_row_sum = row_sums.elementMax()
    val min_row_sum = row_sums.elementMin()
    
    if (max_row_sum > 1e-14 || min_row_sum < -1e-14) {
        throw IllegalArgumentException("The matrix A0+A1+A2 must have zero row sum")
    }
    
    // Note: Eigenvalue checks would require implementing Matrix.eig() 
    // For now, we skip these advanced validation checks
    // In production, these would verify the matrix properties needed for QBD convergence
}