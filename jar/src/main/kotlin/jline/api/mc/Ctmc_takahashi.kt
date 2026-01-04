package jline.api.mc

import jline.util.matrix.Matrix

/**
 * Takahashi's aggregation-disaggregation method for CTMCs
 *
 * @param Q Infinitesimal generator matrix
 * @param MS List where MS[i] is the set of rows of Q in macrostate i
 * @param numSteps Number of iterative steps
 * @return Triple of (p: steady-state probabilities, eps: NCD index, epsMAX: max acceptable eps)
 */
fun ctmc_takahashi(Q: Matrix, MS: List<List<Int>>, numSteps: Int): Triple<Matrix, Double, Double> {
    val nMacroStates = MS.size
    val nStates = Q.numRows

    // Start from Courtois decomposition solution
    val (pcourt, eps, epsMAX) = ctmc_courtois(Q, MS)
    val (P, _) = ctmc_randomization(Q)
    
    // Initialize with Courtois solution
    var pn = pcourt.copy()
    
    // Main loop
    for (n in 0 until numSteps) {
        val pn_1 = pn.copy()
        
        // Aggregation step
        val G = Matrix.zeros(nMacroStates, nMacroStates)
        var procRows = 0
        for (I in 0 until nMacroStates) {
            val blockSizeI = MS[I].size
            var S = 0.0
            for (i in 0 until blockSizeI) {
                S += pn_1.get(procRows + i, 0)
            }
            
            var procCols = 0
            for (J in 0 until nMacroStates) {
                val blockSizeJ = MS[J].size
                if (I != J) {
                    for (i in 0 until blockSizeI) {
                        for (j in 0 until blockSizeJ) {
                            if (S > 1e-14) {
                                G.set(I, J, G.get(I, J) + P.get(procRows + i, procCols + j) * pn_1.get(procRows + i, 0) / S)
                            }
                        }
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
        
        // Compute macroprobabilities
        val gamma = dtmc_solve(G)
        
        // Disaggregation step
        val GI = Matrix.zeros(nMacroStates, nStates)
        procRows = 0
        for (I in 0 until nMacroStates) {
            val blockSizeI = MS[I].size
            var S = 0.0
            for (i in 0 until blockSizeI) {
                S += pn_1.get(procRows + i, 0)
            }
            
            for (j in 0 until nStates) {
                var sum = 0.0
                for (i in 0 until blockSizeI) {
                    sum += P.get(procRows + i, j) * pn_1.get(procRows + i, 0)
                }
                if (S > 0) {
                    GI.set(I, j, sum / S)
                }
            }
            procRows += blockSizeI
        }
        
        // Solve for new probabilities
        procRows = 0
        for (I in 0 until nMacroStates) {
            val blockSize = MS[I].size
            val A = Matrix.eye(blockSize)
            val b = Matrix.zeros(blockSize, 1)
            
            // Set up linear system
            for (i in 0 until blockSize) {
                for (j in 0 until blockSize) {
                    A.set(i, j, A.get(i, j) - P.get(procRows + j, procRows + i))
                }
                
                for (K in 0 until nMacroStates) {
                    if (K != I) {
                        b.set(i, 0, b.get(i, 0) + gamma.get(K, 0) * GI.get(K, procRows + i))
                    }
                }
            }
            
            // Solve A * x = b
            val x = Matrix.zeros(blockSize, 1)
            Matrix.solve(A, b, x)
            for (i in 0 until blockSize) {
                pn.set(procRows + i, 0, x.get(i, 0))
            }
            
            procRows += blockSize
        }
        
        // Normalize
        var sum = 0.0
        for (i in 0 until nStates) {
            sum += pn.get(i, 0)
        }
        for (i in 0 until nStates) {
            pn.set(i, 0, pn.get(i, 0) / sum)
        }
    }

    return Triple(pn, eps, epsMAX)
}
/**
 * CTMC takahashi algorithms
 */
@Suppress("unused")
class CtmcTakahashiAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}