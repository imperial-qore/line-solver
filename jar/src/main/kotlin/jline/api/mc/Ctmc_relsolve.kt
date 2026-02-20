package jline.api.mc

import jline.util.matrix.Matrix

/**
 * Equilibrium distribution of a continuous-time Markov chain re-normalized with respect to 
 * the probability of a reference state.
 * 
 * @param Q Infinitesimal generator matrix of the CTMC
 * @param refstate Reference state (default: 0 for 0-indexed)
 * @param options Optional parameters (not implemented yet)
 * @return Array containing [p, Q, nConnComp, connComp] where:
 *         p: equilibrium distribution
 *         Q: processed generator matrix  
 *         nConnComp: number of connected components
 *         connComp: component assignments
 */
fun ctmc_relsolve(Q: Matrix, refstate: Int = 0, options: Map<String, Any>? = null): Array<Any> {
    var Qmat = Q.copy()
    
    // Size warning for large matrices (>6000 elements)
    if (Qmat.length() > 6000) {
        println("ctmc_relsolve: the order of Q is large (${Qmat.length()}).")
    }
    
    // Handle trivial case
    if (Qmat.length() == 1) {
        val p = Matrix(1, 1)
        p[0, 0] = 1.0
        return arrayOf(p, Qmat, 1, Matrix.ones(1, 1))
    }
    
    // Make proper infinitesimal generator
    Qmat = ctmc_makeinfgen(Qmat)
    val n = Qmat.length()
    
    // Check connectivity: B = abs(Q+Q')>0
    var B = Qmat.add(1.0, Qmat.transpose())
    B.absEq()
    for (colIdx in 0..<B.numCols) {
        val col1 = B.colIndexes[colIdx]
        val col2 = B.colIndexes[colIdx + 1]
        
        for (i in col1..<col2) {
            if (B.nonZeroValues[i] > 0) B.nonZeroValues[i] = 1.0
        }
    }
    
    // Find weakly connected components
    val sets = Matrix.weaklyConnect(B, null)
    val nConnComp = sets.size
    
    // Create component assignment matrix
    val connComp = Matrix(1, n)
    for ((compIdx, set) in sets.withIndex()) {
        for (nodeIdx in set) {
            connComp[0, nodeIdx] = compIdx.toDouble()
        }
    }
    
    // Handle reducible generator - solve each component recursively
    if (nConnComp > 1) {
        val p = Matrix(1, n)
        
        for ((c, set_c) in sets.withIndex()) {
            // Extract subgenerator Qc = Q(connComp==c, connComp==c)
            val Qc = Matrix(set_c.size, set_c.size)
            var Qc_row = 0
            for (q_row in set_c) {
                var Qc_col = 0
                for (q_col in set_c) {
                    Qc[Qc_row, Qc_col] = Qmat[q_row, q_col]
                    Qc_col++
                }
                Qc_row++
            }
            
            // Solve the component recursively
            val QcProper = ctmc_makeinfgen(Qc)
            val pc = ctmc_solve(QcProper)
            
            // Assign solution back to main probability vector
            var idx = 0
            for (i in set_c) {
                p[0, i] = pc[0, idx]
                idx++
            }
        }
        
        // Normalize
        val pSum = p.sumRows(0)
        p.divide(pSum, p, true)
        
        return arrayOf(p, Qmat, nConnComp, connComp)
    }
    
    // Handle all-zero generator
    if (Qmat.nonZeroLength == 0) {
        val p = Matrix(1, n)
        p.fill(1.0 / n)
        return arrayOf(p, Qmat, nConnComp, connComp)
    }
    
    // Main solution algorithm - similar to ctmc_solve but with reference state normalization
    val p = Matrix(1, n)
    val b = Matrix(n, 1)
    
    var nnzel = Matrix(1, n)
    for (i in 0..<n) nnzel[0, i] = i.toDouble()
    var Qnnz = Qmat.copy()
    var bnnz = b.copy()
    var Qnnz_1 = Qmat.copy()
    var bnnz_1 = bnnz.copy()
    
    var isReducible = false
    var goon = true
    
    while (goon) {
        // Find non-zero elements: nnzel = find(sum(abs(Qnnz),1)~=0 & sum(abs(Qnnz),2)'~=0)
        val Qnnz_abs = Qnnz.copy()
        Qnnz_abs.absEq()
        val Qnnz_abs_sum_col = Qnnz_abs.sumCols()
        val Qnnz_abs_sum_rows = Qnnz_abs.sumRows()
        val find_res = Matrix(1, Qnnz_abs_sum_col.numCols)
        
        for (i in 0..<Qnnz_abs_sum_col.numCols) {
            if (Qnnz_abs_sum_col[i] != 0.0 && Qnnz_abs_sum_rows[i] != 0.0) {
                find_res[0, i] = 1.0
            }
        }
        nnzel = find_res.find().transpose()
        
        if (nnzel.length() < n && !isReducible) {
            isReducible = true
        }
        
        // Extract non-zero submatrix: Qnnz = Qnnz(nnzel, nnzel)
        val new_Qnnz = Matrix(nnzel.numCols, nnzel.numCols)
        for (i in 0..<nnzel.numCols) {
            for (j in 0..<nnzel.numCols) {
                val matrixValue = Qnnz[nnzel[0, i].toInt(), nnzel[0, j].toInt()]
                if (matrixValue != 0.0) new_Qnnz[i, j] = matrixValue
            }
        }
        Qnnz = new_Qnnz
        
        // Extract corresponding b vector: bnnz = bnnz(nnzel)
        val new_bnnz = Matrix(nnzel.numCols, 1)
        for (i in 0..<nnzel.numCols) {
            new_bnnz[i, 0] = bnnz[nnzel[0, i].toInt(), 0]
        }
        bnnz = new_bnnz
        
        // Make proper infinitesimal generator
        Qnnz = ctmc_makeinfgen(Qnnz)
        
        // Check convergence
        if ((Qnnz.numCols * Qnnz.numRows == Qnnz_1.numCols * Qnnz_1.numRows) && 
            (bnnz.numCols * bnnz.numRows == bnnz_1.numCols * bnnz_1.numRows)) {
            goon = false
        } else {
            Qnnz_1 = Qnnz.copy()
            bnnz_1 = bnnz.copy()
            nnzel = Matrix(1, Qnnz.length())
            for (i in 0..<Qnnz.length()) nnzel[0, i] = i.toDouble()
        }
    }
    
    if (Qnnz.isEmpty) {
        p.fill(1.0 / n)
        return arrayOf(p, Qmat, nConnComp, connComp)
    }
    
    // Set up reference state normalization: Qnnz(:,end) = 0; Qnnz(refstate,end) = 1; bnnz(end) = 1
    val lastCol = Qnnz.numCols - 1
    val adjustedRefstate = if (refstate < Qnnz.numRows) refstate else 0

    for (i in 0..<Qnnz.numRows) {
        Qnnz[i, lastCol] = 0.0
    }
    Qnnz[adjustedRefstate, lastCol] = 1.0
    bnnz[lastCol, 0] = 1.0
    
    // Solve linear system: p(nnzel) = Qnnz' \ bnnz
    val solve_res = Matrix(0, 0)
    Matrix.solveSafe(Qnnz.transpose(), bnnz, solve_res)
    
    for (i in 0..<nnzel.numCols) {
        p[0, nnzel[0, i].toInt()] = solve_res[i, 0]
    }
    
    return arrayOf(p, Qmat, nConnComp, connComp)
}
/**
 * CTMC relsolve algorithms
 */
@Suppress("unused")
class CtmcRelsolveAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}