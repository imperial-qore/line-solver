package jline.api.mc

import jline.util.matrix.Matrix

/**
 * Koury-McAllister-Stewart aggregation-disaggregation method for CTMCs
 *
 * @param Q Infinitesimal generator matrix
 * @param MS List where MS[i] is the set of rows of Q in macrostate i
 * @param numSteps Number of iterative steps
 * @return Triple of (p: steady-state probabilities, eps: NCD index, epsMAX: max acceptable eps)
 */
fun ctmc_kms(Q: Matrix, MS: List<List<Int>>, numSteps: Int): Triple<Matrix, Double, Double> {
    val nMacroStates = MS.size
    val nStates = Q.numRows

    // Start from Courtois decomposition solution
    val (pcourt, eps, epsMAX) = ctmc_courtois(Q, MS)
    
    // Apply randomization to get P matrix
    val qVal = 1.05 * Q.elementMaxAbs()
    val (P, _) = ctmc_randomization(Q, qVal)
    
    // Rearrange P according to macrostates
    val v = mutableListOf<Int>()
    for (macroState in MS) {
        v.addAll(macroState)
    }
    val Pperm = Matrix.zeros(nStates, nStates)
    for (i in 0 until nStates) {
        for (j in 0 until nStates) {
            Pperm.set(i, j, P.get(v[i], v[j]))
        }
    }
    
    // Initialize with Courtois solution
    var pn = pcourt.copy()
    
    // Main loop
    for (n in 0 until numSteps) {
        val pn_1 = pn.transpose()
        
        // Aggregation step
        val pcondn_1 = Matrix.zeros(1, nStates)
        var procRows = 0
        for (I in 0 until nMacroStates) {
            val blockSize = MS[I].size
            var sum = 0.0
            for (j in 0 until blockSize) {
                sum += pn_1.get(0, procRows + j)
            }
            if (sum > 0) {
                for (j in 0 until blockSize) {
                    pcondn_1.set(0, procRows + j, pn_1.get(0, procRows + j) / sum)
                }
            }
            procRows += blockSize
        }
        
        // Compute aggregated transition matrix G
        val G = Matrix.zeros(nMacroStates, nMacroStates)
        var procCols = 0
        for (I in 0 until nMacroStates) {
            procRows = 0
            for (J in 0 until nMacroStates) {
                var sum = 0.0
                for (j in 0 until MS[J].size) {
                    for (i in 0 until MS[I].size) {
                        sum += pcondn_1.get(0, procRows + j) * Pperm.get(procRows + j, procCols + i)
                    }
                }
                G.set(I, J, sum)
                procRows += MS[J].size
            }
            procCols += MS[I].size
        }
        
        val w = dtmc_solve(G.transpose())
        
        // Disaggregation step
        val z = pcondn_1.copy()
        val zn = Matrix.zeros(1, nStates)
        val L = Matrix.zeros(nStates, nStates)
        val D = Matrix.zeros(nStates, nStates)
        val U = Matrix.zeros(nStates, nStates)
        
        procRows = 0
        for (I in 0 until nMacroStates) {
            val blockSizeI = MS[I].size
            for (j in 0 until blockSizeI) {
                zn.set(0, procRows + j, w.get(I, 0) * z.get(0, procRows + j))
            }
            
            procCols = 0
            for (J in 0 until nMacroStates) {
                val blockSizeJ = MS[J].size
                when {
                    I > J -> {
                        // Fill L matrix
                        for (i in 0 until blockSizeI) {
                            for (j in 0 until blockSizeJ) {
                                L.set(procRows + i, procCols + j, Pperm.get(procRows + i, procCols + j))
                            }
                        }
                    }
                    I == J -> {
                        // Fill D matrix
                        for (i in 0 until blockSizeI) {
                            for (j in 0 until blockSizeJ) {
                                D.set(procRows + i, procCols + j, if (i == j) {
                                    1.0 - Pperm.get(procRows + i, procCols + j)
                                } else {
                                    -Pperm.get(procRows + i, procCols + j)
                                })
                            }
                        }
                    }
                    I < J -> {
                        // Fill U matrix
                        for (i in 0 until blockSizeI) {
                            for (j in 0 until blockSizeJ) {
                                U.set(procRows + i, procCols + j, Pperm.get(procRows + i, procCols + j))
                            }
                        }
                    }
                }
                procCols += blockSizeJ
            }
            procRows += blockSizeI
        }
        
        // Solve system using block decomposition
        val M = D.sub(U)
        val MPL = (nStates - 2) / 2
        
        // Extract blocks
        val A = Matrix.extract(M, 0, MPL + 1, 0, MPL + 1)
        val B = Matrix.extract(M, 0, MPL + 1, MPL + 1, nStates)
        val C = Matrix.extract(M, MPL + 1, nStates, MPL + 1, nStates)
        
        val invA = A.inv()
        val invC = C.inv()
        
        // Compute solution
        val solutionMatrix = Matrix.zeros(nStates, nStates)
        
        // Upper-left block: invA
        for (i in 0 until invA.numRows) {
            for (j in 0 until invA.numCols) {
                solutionMatrix.set(i, j, invA.get(i, j).toDouble())
            }
        }
        
        // Upper-right block: -invA*B*invC
        val temp = invA.mult(B).mult(invC).scale(-1.0)
        for (i in 0 until temp.numRows) {
            for (j in 0 until temp.numCols) {
                solutionMatrix.set(i, MPL + 1 + j, temp.get(i, j).toDouble())
            }
        }
        
        // Lower-right block: invC
        for (i in 0 until invC.numRows) {
            for (j in 0 until invC.numCols) {
                solutionMatrix.set(MPL + 1 + i, MPL + 1 + j, invC.get(i, j).toDouble())
            }
        }
        
        pn = zn.mult(L).mult(solutionMatrix)
    }
    
    // Reorder back to original ordering
    val pout = Matrix.zeros(nStates, 1)
    for (i in 0 until nStates) {
        pout.set(v[i], 0, pn.get(0, i))
    }

    return Triple(pout, eps, epsMAX)
}
/**
 * CTMC kms algorithms
 */
@Suppress("unused")
class CtmcKmsAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}