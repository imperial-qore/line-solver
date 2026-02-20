/**
 * @file CTMC stochastic complementarity analysis
 *
 * Implements stochastic complementarity analysis for CTMCs to identify strongly
 * connected components and analyze their communication structure. Essential for
 * decomposing large Markov chains and understanding long-run behavior.
 *
 * @since LINE 3.0
 */
package jline.api.mc

import jline.solvers.ctmc.SolverCTMC
import jline.util.matrix.Matrix
import kotlin.math.ceil
import org.ejml.data.DMatrixRMaj
import org.ejml.data.DMatrixSparseCSC
import org.ejml.ops.DConvertMatrixStruct
import org.ejml.dense.row.factory.LinearSolverFactory_DDRM
import org.ejml.dense.row.CommonOps_DDRM

fun ctmc_stochcomp(Q: Matrix, I_list: MutableList<Double?>): SolverCTMC.StochCompResult {

    var I = Matrix(I_list)
    if (I_list.isEmpty()) {
        // Match MATLAB: I = 1:ceil(length(Q)/2) -> 0-indexed: 0:ceil(n/2)-1
        val halfSize = ceil(Q.getNumCols().toDouble() / 2).toInt()
        val defaultList: MutableList<Double?> = ArrayList<Double?>()
        for (idx in 0..<halfSize) {
            defaultList.add(idx.toDouble())
        }
        I = Matrix(defaultList)
    }

    // Compute Ic = setdiff(0..n-1, I) using HashSet for O(n) instead of O(n*|I|)
    val iSet = HashSet<Int>(I.getNumRows() * 2)
    for (iRow in 0..<I.getNumRows()) {
        iSet.add(I.get(iRow).toInt())
    }
    val diff_values: MutableList<Double?> = ArrayList<Double?>()
    for (checkValue in 0..<Q.getNumCols()) {
        if (!iSet.contains(checkValue)) {
            diff_values.add(checkValue.toDouble())
        }
    }
    val Ic = Matrix(diff_values)
    val Q11: Matrix
    val Q12: Matrix
    val Q21: Matrix?
    val Q22: Matrix

    Q11 = Q.getSubMatrix(I, I)
    Q12 = Q.getSubMatrix(I, Ic)
    Q21 = Q.getSubMatrix(Ic, I)
    Q22 = Q.getSubMatrix(Ic, Ic)

    // Solve (-Q22) * T = Q21 using dense LU factorization for performance.
    // EJML's sparse LU lacks fill-reducing ordering, making it much slower than
    // MATLAB's UMFPACK for matrices with significant fill-in. Dense LU is predictable
    // and fast for the typical stochcomp sizes (< 2000).
    // Convert Q22 directly to dense and negate there (avoids slow sparse copy).
    val n = Q22.getNumRows()
    val nrhs = Q21!!.getNumCols()

    val denseNegQ22 = DMatrixRMaj(n, n)
    DConvertMatrixStruct.convert(Q22.getData() as DMatrixSparseCSC, denseNegQ22)
    CommonOps_DDRM.scale(-1.0, denseNegQ22)

    // LU factorize (reusable for per-event loop)
    val luSolver = LinearSolverFactory_DDRM.lu(n)
    val luOk = luSolver.setA(denseNegQ22)

    var T: Matrix
    if (luOk) {
        // Dense solve for all columns of Q21
        val denseQ21 = DMatrixRMaj(n, nrhs)
        DConvertMatrixStruct.convert(Q21.getData() as DMatrixSparseCSC, denseQ21)
        val denseT = DMatrixRMaj(n, nrhs)
        luSolver.solve(denseQ21, denseT)

        // Convert back to sparse
        val sparseT = DConvertMatrixStruct.convert(denseT, null as DMatrixSparseCSC?, 1e-15)
        T = Matrix(sparseT as org.ejml.data.DMatrix)
    } else {
        // LU failed (singular matrix), fall back to pseudoinverse
        val negQ22sparse = Q22.neg()
        val negQ22pinv = negQ22sparse.pinv()
        T = negQ22pinv.mult(Q21)
    }
    T = Q12.mult(T)
    val S = Q11.add(1.0, T)
    val result = SolverCTMC.StochCompResult(S, Q11, Q12, Q21, Q22, T)
    // Store the pre-factored LU solver for reuse in per-event stochcomp
    if (luOk) {
        result.denseLUSolver = luSolver
    }
    return result
}
/**
 * CTMC stochcomp algorithms
 */
@Suppress("unused")
class CtmcStochcompAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}
