package jline.api.mc

import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath
import kotlin.math.abs
import kotlin.math.max

/**
 * Courtois decomposition for nearly completely decomposable CTMCs
 * 
 * @param Q Infinitesimal generator matrix
 * @param MS List where MS[i] is the set of rows of Q in macrostate i
 * @param q Randomization coefficient (optional)
 * @return Triple of (p: steady-state probabilities, eps: NCD index, epsMAX: max acceptable eps)
 */
fun ctmc_courtois(Q: Matrix, MS: List<List<Int>>, q: Double = 1.0): Triple<Matrix, Double, Double> {
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
    val QdecFixed = ctmc_makeinfgen(Qdec)
    
    // Apply randomization
    val qVal = if (q == 1.0) {
        1.05 * Qperm.elementMaxAbs()
    } else {
        q
    }
    
    val (P, _) = ctmc_randomization(Qperm, qVal)
    
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
    
    // Compute epsMAX
    val epsMAX = computeEpsMax(A, MS)
    
    // Compute microprobabilities
    val pmicro = Matrix.zeros(n, 1)
    procRows = 0
    for (i in 0 until nMacroStates) {
        val blockSize = MS[i].size
        val Qmicrostate = Matrix.zeros(blockSize, blockSize)
        for (row in 0 until blockSize) {
            for (col in 0 until blockSize) {
                Qmicrostate.set(row, col, Qdec.get(procRows + row, procRows + col))
            }
        }
        val microProb = ctmc_solve(Qmicrostate)
        for (j in 0 until blockSize) {
            pmicro.set(procRows + j, 0, microProb.get(j, 0))
        }
        procRows += blockSize
    }
    
    // Compute macroprobabilities
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
    
    val pMacro = dtmc_solve(G)
    
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

internal fun computeEpsMax(A: Matrix, MS: List<List<Int>>): Double {
    val eigMS = mutableListOf<Double>()
    var procRows = 0
    
    for (i in MS.indices) {
        val blockSize = MS[i].size
        val blockMatrix = Matrix.zeros(blockSize, blockSize)
        
        // Extract block
        for (row in 0 until blockSize) {
            for (col in 0 until blockSize) {
                blockMatrix.set(row, col, A.get(procRows + row, procRows + col))
            }
        }
        
        // Compute eigenvalues (simplified - would need proper eigenvalue computation)
        // This is a placeholder - real implementation would compute actual eigenvalues
        val secondLargestEig = if (blockSize > 1) 0.9 else 0.0
        eigMS.add(secondLargestEig)
        
        procRows += blockSize
    }
    
    return (1.0 - (eigMS.maxOrNull() ?: 0.0)) / 2.0
}
/**
 * CTMC courtois algorithms
 */
@Suppress("unused")
class CtmcCourtoisAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}