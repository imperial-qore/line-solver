/**
 * @file NSF_GHT - Non-Skip-Free Gail-Hantler-Taylor Algorithm
 *
 * Computes the minimal nonnegative solution G for Non-Skip-Free (NSF)
 * Markov chains using the Gail-Hantler-Taylor algorithm.
 *
 * For M/G/1-type Markov chains where multiple transitions can occur
 * per step (e.g., MAP/D/c queues where multiple arrivals are possible
 * in one service interval).
 *
 * Based on the SMC Solver implementation by Benny Van Houdt.
 *
 * @since LINE 3.1.0
 */
package jline.lib.smc

import jline.io.line_warning
import jline.util.matrix.Matrix
import kotlin.math.abs

/**
 * Options for NSF_GHT solver
 */
data class NSFGHTOptions(
    val maxNumIt: Int = 10000,
    val verbose: Int = 0,
    val firstBlockRow: Boolean = false
)

/**
 * Gail-Hantler-Taylor algorithm for Non-Skip-Free Markov chains.
 *
 * Computes the minimal nonnegative solution to the matrix equation:
 *   G = C0 + C1*G + C2*G^2 + C3*G^3 + ... + C_max*G^max
 *
 * where the Ci matrices are derived from reblocking the original A matrices.
 *
 * @param A Block matrix [A0 A1 A2 ... A_max] with m rows and m*(max+1) columns.
 *          Must be nonnegative with (A0+A1+...+A_max) irreducible and stochastic.
 * @param N Number of non-zero blocks below the main diagonal (the skip factor).
 *          For example, N=c for MAP/D/c where up to c customers can be served per interval.
 * @param options Solver options
 * @return G matrix of dimension (m*N) x (m*N)
 */
fun nsfGht(
    A: Matrix,
    N: Int,
    options: NSFGHTOptions = NSFGHTOptions()
): Matrix {
    val m = A.numRows
    val K = A.numCols / m - 1  // degree of A(z)

    require(N >= 1) { "N must be at least 1" }
    require(K >= 0) { "A must have at least one block" }

    var numit = 0
    var check = 1.0

    // Initialize G = A(:, 1:m*N) = [A0 A1 ... A_{N-1}]
    val initialCols = minOf(m * N, A.numCols)
    var G = if (initialCols > 0) {
        A.extractCols(0, initialCols)
    } else {
        Matrix.zeros(m, m * N)
    }
    // Ensure G has correct dimensions
    if (G.numCols < m * N) {
        val paddedG = Matrix.zeros(m, m * N)
        for (i in 0 until m) {
            for (j in 0 until G.numCols) {
                paddedG[i, j] = G[i, j]
            }
        }
        G = paddedG
    }

    // Iterative computation
    while (check > 1e-14 && numit < options.maxNumIt) {
        numit++
        val Gold = G.copy()
        var temp = G.copy()

        // G = A(:,1:m*N) + A(:,m*N+1:m*(N+1)) * G
        val AN = if ((N + 1) * m <= A.numCols) {
            A.extractCols(N * m, (N + 1) * m)
        } else {
            Matrix.zeros(m, m)
        }
        G = A.extractCols(0, m * N).add(AN.mult(G))

        // For j = N+1 to K
        for (j in (N + 1)..K) {
            // temp = [zeros(m) temp(:,1:(N-1)*m)] + temp(:,(N-1)*m+1:end) * Gold
            val newTemp = Matrix.zeros(m, m * N)

            // Copy temp(:,1:(N-1)*m) shifted right by m
            for (i in 0 until m) {
                for (col in 0 until (N - 1) * m) {
                    newTemp[i, col + m] = temp[i, col]
                }
            }

            // Add temp(:,(N-1)*m+1:end) * Gold
            val tempRight = temp.extractCols((N - 1) * m, N * m)
            val product = tempRight.mult(Gold)
            for (i in 0 until m) {
                for (col in 0 until m * N) {
                    newTemp[i, col] = newTemp[i, col] + product[i, col]
                }
            }
            temp = newTemp

            // G = G + A(:,j*m+1:(j+1)*m) * temp
            if ((j + 1) * m <= A.numCols) {
                val Aj = A.extractCols(j * m, (j + 1) * m)
                G = G.add(Aj.mult(temp))
            }
        }

        check = Matrix.infNorm(G.sub(Gold))
    }

    if (numit == options.maxNumIt && check > 1e-14) {
        line_warning("NSF_GHT", "Maximum Number of Iterations %d reached", numit)
    }

    // If not firstBlockRow, compute remaining block rows of G
    if (!options.firstBlockRow && N > 1) {
        // Expand G to full (m*N) x (m*N) matrix
        val fullG = Matrix.zeros(m * N, m * N)

        // Copy first block row
        for (i in 0 until m) {
            for (j in 0 until m * N) {
                fullG[i, j] = G[i, j]
            }
        }

        // Compute remaining block rows
        // G((j-1)*m+1:j*m,:) = [zeros(m) G((j-2)*m+1:(j-1)*m,1:(N-1)*m)] +
        //                      G((j-2)*m+1:(j-1)*m,(N-1)*m+1:end) * G(1:m,:)
        val G1 = G.copy()  // First block row

        for (j in 2..N) {
            val prevRow = Matrix(m, m * N)
            for (i in 0 until m) {
                for (c in 0 until m * N) {
                    prevRow[i, c] = fullG[(j - 2) * m + i, c]
                }
            }

            // [zeros(m) prevRow(:,1:(N-1)*m)]
            val shiftedPart = Matrix.zeros(m, m * N)
            for (i in 0 until m) {
                for (c in 0 until (N - 1) * m) {
                    shiftedPart[i, c + m] = prevRow[i, c]
                }
            }

            // prevRow(:,(N-1)*m+1:end) * G(1:m,:)
            val rightPart = prevRow.extractCols((N - 1) * m, N * m)
            val multPart = rightPart.mult(G1)

            // Sum
            for (i in 0 until m) {
                for (c in 0 until m * N) {
                    fullG[(j - 1) * m + i, c] = shiftedPart[i, c] + multPart[i, c]
                }
            }
        }

        G = fullG
    }

    return G
}

/**
 * Compute residual error for verification
 */
private fun computeResidual(A: Matrix, G: Matrix, N: Int, m: Int): Double {
    // Pad A if needed
    val K = A.numCols / m - 1
    val extraBlocks = N - 1 - (K % N)
    val paddedA = if (extraBlocks > 0) {
        val newA = Matrix.zeros(m, A.numCols + m * extraBlocks)
        for (i in 0 until m) {
            for (j in 0 until A.numCols) {
                newA[i, j] = A[i, j]
            }
        }
        newA
    } else {
        A
    }

    val Nb = paddedA.numCols / (m * N) - 1
    var Gcheck = paddedA.extractCols(Nb * N * m, paddedA.numCols)

    // First block row of G for multiplication
    val G1 = Matrix(m, m * N)
    for (i in 0 until m) {
        for (j in 0 until m * N) {
            G1[i, j] = G[i, j]
        }
    }

    for (j in (Nb - 1) downTo 0) {
        Gcheck = paddedA.extractCols(j * N * m, (j + 1) * N * m).add(Gcheck.mult(G))
    }

    return Matrix.infNorm(G1.sub(Gcheck))
}

/**
 * NSF_pi - Stationary vector of a Non-Skip-Free Markov chain
 *
 * Computes the stationary vector for an NSF chain using the computed G matrix.
 *
 * @param B Boundary block matrix (or null to use homogeneous case)
 * @param A Repeating block matrix [A0 A1 ... A_max]
 * @param G The G matrix computed by nsfGht
 * @param options Solver options
 * @return Stationary distribution vector
 */
fun nsfPi(
    B: Matrix?,
    A: Matrix,
    G: Matrix,
    options: NSFPiOptions = NSFPiOptions()
): Matrix {
    val m = A.numRows
    val K = A.numCols / m - 1
    val N = G.numCols / m

    // Ensure G has full block rows if firstBlockRow was used
    val fullG = if (options.firstBlockRow && G.numRows == m) {
        expandGBlockRows(G, N, m)
    } else {
        G
    }

    // Construct reblocked C matrices
    val extraBlocksC = N - 1 - ((K - 1) % N)
    val CWidth = m * (N - 1) + A.numCols + m * extraBlocksC
    val C = Matrix.zeros(N * m, CWidth)

    // Fill C with shifted copies of A
    for (row in 0 until m) {
        for (col in 0 until (K + 1) * m) {
            C[row, col] = A[row, col]
        }
    }
    for (i in 1 until N) {
        for (row in 0 until m) {
            for (col in 0 until (K + 1) * m) {
                C[i * m + row, i * m + col] = A[row, col]
            }
        }
    }

    if (B != null) {
        // Use recursive call with B
        val extraBlocksB = N - 1 - (K % N)
        val paddedB = if (extraBlocksB > 0) {
            val newB = Matrix.zeros(N * m, B.numCols + m * extraBlocksB)
            for (i in 0 until B.numRows) {
                for (j in 0 until B.numCols) {
                    newB[i, j] = B[i, j]
                }
            }
            newB
        } else {
            B
        }
        return mg1_pi(paddedB, C, MG1PiOptions(maxNumComp = options.maxNumComp / N))
    } else {
        // Homogeneous case
        // Compute FirstRowSumsG
        val firstRowSumsG = Matrix(m, N * m)
        for (i in 0 until m) {
            for (j in 0 until m) {
                firstRowSumsG[i, j] = fullG[i, j]
            }
        }
        for (block in 1 until N) {
            for (i in 0 until m) {
                for (j in 0 until m) {
                    firstRowSumsG[i, block * m + j] = firstRowSumsG[i, (block - 1) * m + j] + fullG[i, block * m + j]
                }
            }
        }

        // ghat = stat(FirstRowSumsG(:,(N-1)*m+1:end))
        val lastBlock = firstRowSumsG.extractCols((N - 1) * m, N * m)
        val ghat = stat(lastBlock)

        // pi0 = ghat * G(1:m,:)
        val G1Row = fullG.extractRows(0, m)
        var pi0 = ghat.mult(G1Row)

        // Normalize using first m components of mu1
        val g = ghat.mult(firstRowSumsG)

        // beta = Cmax*e + (Cmax + Cmax-1)*e + ... + (Cmax + ... + C1)*e
        val beta = Matrix.zeros(m * N, 1)
        var CSum = Matrix.zeros(m * N, m * N)
        val numCBlocks = CWidth / (m * N)

        for (i in (numCBlocks - 1) downTo 1) {
            for (r in 0 until m * N) {
                for (c in 0 until m * N) {
                    CSum[r, c] = CSum[r, c] + C[r, i * m * N + c]
                }
            }
            val ones = Matrix.ones(m * N, 1)
            val CSumOnes = CSum.mult(ones)
            for (r in 0 until m * N) {
                beta[r, 0] = beta[r, 0] + CSumOnes[r, 0]
            }
        }

        // Add C0 to CSum
        for (r in 0 until m * N) {
            for (c in 0 until m * N) {
                CSum[r, c] = CSum[r, c] + C[r, c]
            }
        }

        // temp = (eye(m*N) - CSum + (ones(m*N,1) - beta) * g)^(-1) * ones(N*m,1)
        val eye = Matrix.eye(m * N)
        val ones = Matrix.ones(m * N, 1)
        val onesMinusBeta = ones.sub(beta)
        val outerProduct = onesMinusBeta.mult(g)
        val matrixToInvert = eye.sub(CSum).add(outerProduct)
        val invMatrix = matrixToInvert.inv()
        var temp = invMatrix.mult(ones)

        // temp = ([eye(m) zeros(m,(N-1)*m)] - G(1:m,:) + ones(m,1) * g) * temp
        val eyeZeros = Matrix(m, N * m)
        for (i in 0 until m) {
            eyeZeros[i, i] = 1.0
        }
        val onesM = Matrix.ones(m, 1)
        val finalMatrix = eyeZeros.sub(G1Row).add(onesM.mult(g))
        temp = finalMatrix.mult(temp)

        // Normalize pi0
        val normFactor = ghat.mult(temp)[0, 0]
        pi0 = pi0.scale(1.0 / normFactor)

        // Compute higher levels using MG1_pi pattern
        val Bmatrix = Matrix(N * m, A.numCols)
        for (i in 0 until N) {
            for (r in 0 until m) {
                for (c in 0 until A.numCols) {
                    Bmatrix[i * m + r, c] = A[r, c]
                }
            }
        }

        val extraBlocksB = N - 1 - (K % N)
        val paddedB = if (extraBlocksB > 0) {
            val newB = Matrix.zeros(N * m, Bmatrix.numCols + m * extraBlocksB)
            for (i in 0 until Bmatrix.numRows) {
                for (j in 0 until Bmatrix.numCols) {
                    newB[i, j] = Bmatrix[i, j]
                }
            }
            newB
        } else {
            Bmatrix
        }

        return mg1_pi(paddedB, C, MG1PiOptions(
            maxNumComp = options.maxNumComp / N,
            verbose = options.verbose
        ))
    }
}

/**
 * Options for NSF_pi solver
 */
data class NSFPiOptions(
    val maxNumComp: Int = 1000,
    val verbose: Boolean = false,
    val firstBlockRow: Boolean = false
)

/**
 * Expand G from first block row to full matrix
 */
private fun expandGBlockRows(G: Matrix, N: Int, m: Int): Matrix {
    if (G.numRows >= N * m) {
        return G
    }

    val fullG = Matrix.zeros(m * N, m * N)

    // Copy first block row
    for (i in 0 until m) {
        for (j in 0 until m * N) {
            fullG[i, j] = G[i, j]
        }
    }

    // Compute remaining block rows
    for (j in 2..N) {
        val prevRow = Matrix(m, m * N)
        for (i in 0 until m) {
            for (c in 0 until m * N) {
                prevRow[i, c] = fullG[(j - 2) * m + i, c]
            }
        }

        val shiftedPart = Matrix.zeros(m, m * N)
        for (i in 0 until m) {
            for (c in 0 until (N - 1) * m) {
                shiftedPart[i, c + m] = prevRow[i, c]
            }
        }

        val rightPart = prevRow.extractCols((N - 1) * m, N * m)
        val G1 = fullG.extractRows(0, m)
        val multPart = rightPart.mult(G1)

        for (i in 0 until m) {
            for (c in 0 until m * N) {
                fullG[(j - 1) * m + i, c] = shiftedPart[i, c] + multPart[i, c]
            }
        }
    }

    return fullG
}
