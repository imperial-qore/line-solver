package jline.api.mc

import jline.util.matrix.Matrix

/**
 * Multi-level aggregation method for CTMCs
 * Basic one-step implementation of standard multigrid method
 *
 * @param Q Infinitesimal generator matrix
 * @param MS List where MS[i] is the set of rows of Q in macrostate i
 * @param MSS List where MSS[i] is the set of rows of G in macromacrostate i
 * @return Triple of (p: steady-state probabilities, eps: NCD index, epsMAX: max acceptable eps)
 */
fun ctmc_multi(Q: Matrix, MS: List<List<Int>>, MSS: List<List<Int>>): Triple<Matrix, Double, Double> {
    val nMacroStates = MS.size
    val n = Q.numRows
    
    // Rearrange infinitesimal generator according to macrostates
    val v = mutableListOf<Int>()
    for (macroState in MS) {
        v.addAll(macroState)
    }
    
    // Create permuted matrix
    val Qperm = Matrix.zeros(n, n)
    for (i in 0 until n) {
        for (j in 0 until n) {
            Qperm.set(i, j, Q.get(v[i], v[j]))
        }
    }
    
    // Create decomposed matrix
    val Qdec = Qperm.copy()
    var procRows = 0
    for (i in 0 until nMacroStates) {
        val blockSize = MS[i].size
        // Zero out inter-block transitions
        if (procRows > 0) {
            for (row in procRows until procRows + blockSize) {
                for (col in 0 until procRows) {
                    Qdec.set(row, col, 0.0)
                }
            }
        }
        for (row in procRows until procRows + blockSize) {
            for (col in procRows + blockSize until n) {
                Qdec.set(row, col, 0.0)
            }
        }
        procRows += blockSize
    }
    
    // Make each diagonal block a valid infinitesimal generator
    val QdecNew = ctmc_makeinfgen(Qdec)
    
    // Apply randomization
    val q = 1.05 * Qperm.elementMaxAbs()
    val (P, _) = ctmc_randomization(Qperm, q)
    
    // Create block-diagonal matrix A
    val A = P.copy()
    procRows = 0
    for (i in 0 until nMacroStates) {
        val blockSize = MS[i].size
        if (procRows > 0) {
            for (row in procRows until procRows + blockSize) {
                for (col in 0 until procRows) {
                    A.set(row, col, 0.0)
                }
            }
        }
        for (row in procRows until procRows + blockSize) {
            for (col in procRows + blockSize until n) {
                A.set(row, col, 0.0)
            }
        }
        procRows += blockSize
    }
    
    // Compute NCD error index
    val B = P.sub(A)
    val eps = B.sumRows().elementMax()

    // Compute epsMAX (max acceptable NCD index)
    val epsMAX = computeEpsMax(A, MS)
    
    // Compute microprobabilities
    val pmicro = Matrix.zeros(n, 1)
    procRows = 0
    for (i in 0 until nMacroStates) {
        val blockSize = MS[i].size
        val Qmicrostate = Matrix.zeros(blockSize, blockSize)
        for (row in 0 until blockSize) {
            for (col in 0 until blockSize) {
                Qmicrostate.set(row, col, QdecNew.get(procRows + row, procRows + col))
            }
        }
        val microProb = ctmc_solve(Qmicrostate)
        for (j in 0 until blockSize) {
            pmicro.set(procRows + j, 0, microProb.get(j, 0))
        }
        procRows += blockSize
    }
    
    // Compute macroprobabilities transition matrix G
    val G = Matrix.zeros(nMacroStates, nMacroStates)
    procRows = 0
    for (i in 0 until nMacroStates) {
        val blockSizeI = MS[i].size
        var procCols = 0
        for (j in 0 until nMacroStates) {
            val blockSizeJ = MS[j].size
            if (i != j) {
                for (iState in 0 until blockSizeI) {
                    var sum = 0.0
                    for (jState in 0 until blockSizeJ) {
                        sum += P.get(procRows + iState, procCols + jState)
                    }
                    G.set(i, j, G.get(i, j) + pmicro.get(procRows + iState, 0) * sum)
                }
            }
            procCols += blockSizeJ
        }
        procRows += blockSizeI
    }
    
    // Make G stochastic
    for (i in 0 until nMacroStates) {
        var rowSum = 0.0
        for (j in 0 until nMacroStates) {
            if (i != j) {
                rowSum += G.get(i, j)
            }
        }
        G.set(i, i, 1.0 - rowSum)
    }
    
    // Convert G to infinitesimal generator for Courtois method
    val Ginf = G.copy()
    for (i in 0 until nMacroStates) {
        for (j in 0 until nMacroStates) {
            if (i != j) {
                Ginf.set(i, j, G.get(i, j))
            } else {
                Ginf.set(i, i, G.get(i, i) - 1.0)
            }
        }
    }
    
    // Use Courtois method on the macro level with MSS
    val (pMacro, _, _) = ctmc_courtois(Ginf, MSS)
    
    // Compute final probabilities
    val p = Matrix.zeros(n, 1)
    procRows = 0
    for (i in 0 until nMacroStates) {
        val blockSize = MS[i].size
        for (j in 0 until blockSize) {
            p.set(procRows + j, 0, pMacro.get(i, 0) * pmicro.get(procRows + j, 0))
        }
        procRows += blockSize
    }
    
    // Reorder back to original ordering
    val pout = Matrix.zeros(n, 1)
    for (i in 0 until n) {
        pout.set(v[i], 0, p.get(i, 0))
    }

    return Triple(pout, eps, epsMAX)
}
/**
 * CTMC multi algorithms
 */
@Suppress("unused")
class CtmcMultiAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}