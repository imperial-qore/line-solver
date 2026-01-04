/**
 * @file Continuous-time Markov chain steady-state solver
 * 
 * Computes the steady-state probability distribution for CTMCs by solving the linear
 * system π·Q = 0 with normalization π·1 = 1. Handles reducible chains by decomposing
 * into strongly connected components and solving each component separately.
 *
 * @since LINE 3.0
 */
package jline.api.mc

import jline.GlobalConstants
import jline.util.matrix.Matrix

/**
 * Return the steady-state probability of a CTMC.
 *
 * @param Q Infinitesimal generator of the CTMC
 * @return Steady-state probability vector
 */

fun ctmc_solve(Q: Matrix): Matrix {
    var Qmat = Q.copy()
    if (Qmat.length() == 1) {
        val p = Matrix(1, 1, 1)
        p[0, 0] = 1
        return p
    }

    Qmat = ctmc_makeinfgen(Qmat)
    val n = Qmat.length()

    // B = abs(Q+Q')>0
    var B = Qmat.add(1.0, Qmat.transpose())
    B.absEq()
    for (colIdx in 0..<B.numCols) {
        val col1 = B.colIndexes[colIdx]
        val col2 = B.colIndexes[colIdx + 1]

        for (i in col1..<col2) {
            if (B.nonZeroValues[i] > GlobalConstants.Zero) B.nonZeroValues[i] = 1.0
        }
    }

    // [nConnComp, connComp] = weaklyconncomp(B);
    var sets = Matrix.weaklyConnect(B, null)
    if (sets.size > 1) {
        val p = Matrix(1, n)

        for (set_c in sets) {
            // Qc = Q(connComp==c,connComp==c);
            var Qc = Matrix(set_c.size, set_c.size)
            var Qc_row = 0
            var Qc_col = 0
            for (q_row in set_c) {
                for (q_col in set_c) {
                    Qc[Qc_row, Qc_col++] = Qmat[q_row, q_col]
                }
                Qc_row++
                Qc_col = 0
            }
            // Qc = ctmc_makeinfgen(Qc);
            Qc = ctmc_makeinfgen(Qc)
            // p(connComp==c) = ctmc_solve(Qc);
            val ctmc_solve_Qc = ctmc_solve(Qc)
            var idx = 0
            for (i in set_c) {
                p[0, i] = ctmc_solve_Qc[0, idx++]
            }
        }
        p.divide(p.sumRows(0), p, true)
        return p
    }

    if (Qmat.nonZeroLength == 0) {
        val p = Matrix(1, n)
        p.fill(1.0 / n)
        return p
    }

    val p = Matrix(1, n)
    val b = Matrix(n, 1)
    var nnzel = Matrix(1, n)
    for (i in 0..<n) nnzel[0, i] = i
    var Qnnz = Qmat.copy()
    var bnnz = b.copy()
    var Qnnz_1 = Qmat.copy()
    var bnnz_1 = bnnz.copy()

    var isReducible = false
    var goon = true
    while (goon) {
        // nnzel = find(sum(abs(Qnnz),1)~=0 & sum(abs(Qnnz),2)'~=0);
        val Qnnz_abs = Qnnz!!.copy()
        Qnnz_abs.absEq()
        val Qnnz_abs_sum_col = Qnnz_abs.sumCols()
        val Qnnz_abs_sum_rows = Qnnz_abs.sumRows()
        val find_res = Matrix(1, Qnnz_abs_sum_col.numCols)
        for (i in 0..<Qnnz_abs_sum_col.numCols) {
            if (Qnnz_abs_sum_col[i] != 0.0 && Qnnz_abs_sum_rows[i] != 0.0) find_res[0, i] = 1
        }
        nnzel = find_res.find().transpose()

        if (nnzel.length() < n && !isReducible) {
            isReducible = true
            // if (nargin > 1 && options.verbose == 2) % debug
            // fprintf(1,'ctmc_solve: the infinitesimal generator is reducible.\n');
            // end
        }

        // Qnnz = Qnnz(nnzel, nnzel);
        val new_Qnnz = Matrix(nnzel.numCols, nnzel.numCols)
        for (i in 0..<nnzel.numCols) {
            for (j in 0..<nnzel.numCols) {
                val matrixValue = Qnnz[nnzel[0, i].toInt(), nnzel[0, j].toInt()]
                if (matrixValue > 0) new_Qnnz[i, j] = matrixValue
            }
        }
        Qnnz = new_Qnnz

        // bnnz = bnnz(nnzel);
        val new_bnnz = Matrix(nnzel.numCols, 1)
        for (i in 0..<nnzel.numCols) {
            new_bnnz[i, 0] = bnnz[nnzel[0, i].toInt(), 0]
        }
        bnnz = new_bnnz

        // Qnnz = ctmc_makeinfgen(Qnnz);
        Qnnz = ctmc_makeinfgen(Qnnz)

        if ((Qnnz.numCols * Qnnz.numRows == Qnnz_1.numCols * Qnnz_1.numRows) && (bnnz.numCols * bnnz.numRows == bnnz_1.numCols * bnnz_1.numRows)) {
            goon = false
        } else {
            Qnnz_1 = Qnnz.copy()
            bnnz_1 = bnnz.copy()
            nnzel = Matrix(1, Qnnz.length())
            for (i in 0..<Qnnz.length()) nnzel[0, i] = i
        }
    }

    if ((Qnnz == null) || (Qnnz.isEmpty)) {
        p.fill(1.0 / n)
        return p
    }

    // Qnnz(:,end) = 1;
    for (i in 0..<Qnnz.numRows) Qnnz[i, Qnnz.numCols - 1] = 1.0

    // bnnz(end) = 1;
    bnnz[bnnz.numRows - 1, 0] = 1.0

    // p(nnzel)=Qnnz'\ bnnz;
    val solve_res = Matrix(bnnz.numRows, 1)
    // Use solveSafe to match MATLAB behavior - returns NaN for singular matrices instead of throwing
    Matrix.solveSafe(Qnnz.transpose(), bnnz, solve_res)
    for (i in 0..<nnzel.numCols) p[0, nnzel[0, i].toInt()] = solve_res[i, 0]

    if (p.hasNaN()) {
        // B = abs(Qnnz+Qnnz')>0;
        B = Qnnz.add(1.0, Qnnz.transpose())
        B.absEq()
        for (colIdx in 0..<B.numCols) {
            val col1 = B.colIndexes[colIdx]
            val col2 = B.colIndexes[colIdx + 1]

            for (i in col1..<col2) {
                if (B.nonZeroValues[i] > 0) B.nonZeroValues[i] = 1.0
            }
        }

        // [nConnComp, connComp] = weaklyconncomp(B);
        sets = Matrix.weaklyConnect(B, null)
        if (sets.size > 1) {
            // p(nnzel) = zeros(1,n);
            for (i in 0..<nnzel.numCols) p.remove(0, nnzel[0, i].toInt())

            for (set_c in sets) {
                // Qc = Q(connComp==c,connComp==c);
                var Qc = Matrix(set_c.size, set_c.size)
                var Qc_row = 0
                var Qc_col = 0
                for (q_row in set_c) {
                    for (q_col in set_c) {
                        Qc[Qc_row, Qc_col++] = Qmat[q_row, q_col]
                    }
                    Qc_row++
                    Qc_col = 0
                }
                Qc = ctmc_makeinfgen(Qc)
                // p(intersect(find(connComp==c),nnzel)) = ctmc_solve(Qc);
                val ctmc_solve_Qc = ctmc_solve(Qc)
                var idx = 0
                for (i in 0..<nnzel.numCols) {
                    val nnzelValue = nnzel[0, i].toInt()
                    if (set_c.contains(nnzelValue)) p[0, nnzelValue] = ctmc_solve_Qc[0, idx++]
                }
            }
            p.divide(p.sumRows(0), p, true)
            return p
        }
    }
    return p
}
/**
 * CTMC solve algorithms
 */
@Suppress("unused")
class CtmcSolveAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}