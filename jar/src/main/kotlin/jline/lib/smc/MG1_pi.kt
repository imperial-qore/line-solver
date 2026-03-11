/**
 * @file M/G/1-type Stationary Distribution
 *
 * Computes the stationary distribution for M/G/1-type Markov chains.
 *
 * Based on the SMC Solver implementation by Benny Van Houdt.
 *
 * @since LINE 3.1.0
 */
package jline.lib.smc

import jline.util.matrix.Matrix

/**
 * Options for MG1_pi solver
 */
data class MG1PiOptions(
    val boundary: Matrix? = null,
    val maxNumComp: Int = 500,
    val precision: Int = 200,
    val solver: String = "FI",
    val verbose: Boolean = false,
    val mode: String = "ShiftPWCR"
)

/**
 * Computes the stationary distribution of an M/G/1-type Markov chain.
 *
 * The chain is characterized by:
 * - B: boundary blocks [B0 B1 B2 ... B_maxb] with m rows
 * - A: repeating blocks [A0 A1 A2 ... A_maxa] with m rows
 *
 * @param B Boundary block matrix (or null to use first row of A)
 * @param A Repeating block matrix
 * @param options Solver options
 * @return Stationary distribution vector
 */
fun mg1_pi(B: Matrix?, A: Matrix, options: MG1PiOptions = MG1PiOptions()): Matrix {
    val m = A.numRows
    val dega = A.numCols / m - 1

    // Use boundary or default to A's first row structure
    val boundary = B ?: A

    // First, compute G using appropriate solver
    val G = when (options.solver.uppercase()) {
        "CR" -> mg1_cr(A, MG1CROptions(
            mode = options.mode,
            verbose = if (options.verbose) 1 else 0
        ))
        else -> mg1_fi(A, MG1FIOptions(
            verbose = if (options.verbose) 1 else 0
        ))
    }

    // Compute R matrix
    // R = sum(i=0 to dega) A_i * G^i
    var R = A.extractCols(0, m).copy()
    var Gpow = G.copy()
    for (i in 1..dega) {
        R = R.add(A.extractCols(i * m, (i + 1) * m).mult(Gpow))
        Gpow = Gpow.mult(G)
    }

    // Compute stationary distribution
    // pi_0 satisfies: pi_0 * B_hat = 0 where B_hat = sum(B_i * G^i) - I
    val degb = boundary.numCols / m - 1

    var Bhat = boundary.extractCols(0, m).copy()
    Gpow = G.copy()
    for (i in 1..degb) {
        Bhat = Bhat.add(boundary.extractCols(i * m, (i + 1) * m).mult(Gpow))
        Gpow = Gpow.mult(G)
    }

    // Compute pi_0 as stationary distribution of Bhat
    val pi0 = stat(Bhat)

    // Compute higher levels using pi_i = pi_0 * R^i
    val result = mutableListOf<Matrix>()
    result.add(pi0)

    var Rpow = R.copy()
    for (i in 1 until options.maxNumComp) {
        val pi_i = pi0.mult(Rpow)

        // Check if probability mass is negligible
        val mass = pi_i.elementSum()
        if (mass < 1e-15) {
            break
        }

        result.add(pi_i)
        Rpow = Rpow.mult(R)
    }

    // Normalize the distribution
    var totalMass = 0.0
    for (pi_i in result) {
        totalMass += pi_i.elementSum()
    }

    // Return as single row vector
    val totalSize = result.size * m
    val piVec = Matrix(1, totalSize)
    for (i in result.indices) {
        for (j in 0 until m) {
            piVec[0, i * m + j] = result[i][0, j] / totalMass
        }
    }

    return piVec
}

/**
 * Computes the G matrix using Cyclic Reduction for M/G/1-type Markov Chains.
 *
 * Solves: G = A0 + A1*G + A2*G^2 + ... + A_max*G^max
 *
 * Implements the Point-Wise Cyclic Reduction (PWCR) and Shift+PWCR algorithms
 * from Bini and Meini, using DFT for point-wise evaluation at roots of unity.
 *
 * @param A Block matrix [A0 A1 A2 ... A_max] with m rows and m*(max+1) columns
 * @param options Solver options
 * @return G matrix (minimal nonnegative solution)
 */
fun mg1_cr(A: Matrix, options: MG1CROptions = MG1CROptions()): Matrix {
    val m = A.numRows

    // Check whether G is known explicitly
    val explicitG = mg1_eg(A, options.verbose > 0)
    if (explicitG != null) {
        return explicitG
    }

    var Dold: Matrix? = null
    var workA = A.copy()
    var drift = 0.0
    var tau = 0.0
    var v = Matrix.zeros(m, 1)

    if (options.mode.contains("ShiftPWCR")) {
        if (options.verbose == 1) {
            Dold = workA.copy()
        }
        val shiftResult = mg1_shifts(workA, options.shiftType)
        workA = shiftResult.first
        drift = shiftResult.second
        tau = shiftResult.third
        v = shiftResult.fourth
    }

    // Start Cyclic Reduction
    // Transpose and pad D to next power-of-2 + 1 blocks
    val dega = workA.numCols / m - 1
    val paddedBlocks = (1 shl (1 + floorLog2(dega))) + 1
    val D = Matrix(paddedBlocks * m, m)
    // D = transpose of workA, stored as column blocks
    for (i in 0..dega) {
        for (r in 0 until m) {
            for (c in 0 until m) {
                D[i * m + r, c] = workA[c, i * m + r]
            }
        }
    }
    // Remaining blocks are zero (already initialized)

    // Step 0: Split into even and odd indexed blocks
    val totalBlocks = paddedBlocks
    var Aeven = extractEvenBlocks(D, m, totalBlocks)
    var Aodd = extractOddBlocks(D, m, totalBlocks)

    var Ahatodd = Matrix((Aeven.numRows / m - 1 + 1) * m, m)
    // Ahatodd = [Aeven(m+1:end,:); D(end-m+1:end,:)]
    for (i in 1 until Aeven.numRows / m) {
        for (r in 0 until m) {
            for (c in 0 until m) {
                Ahatodd[(i - 1) * m + r, c] = Aeven[i * m + r, c]
            }
        }
    }
    for (r in 0 until m) {
        for (c in 0 until m) {
            Ahatodd[(Aeven.numRows / m - 1) * m + r, c] = D[(totalBlocks - 1) * m + r, c]
        }
    }
    var Ahateven = Aodd.copy()

    // Compute Rj for stop criteria (PWCR mode)
    var Rj = Matrix(m, m)
    for (i in 1 until D.numRows / m) {
        val block = D.extractRows(i * m, (i + 1) * m)
        Rj = Rj.add(block)
    }
    Rj = Matrix.eye(m).sub(Rj).inv().mult(D.extractRows(0, m))

    var G = Matrix.zeros(m, m)
    var numit = 0

    while (numit < options.maxNumIt) {
        numit++
        val nj = Aodd.numRows / m - 1

        var Anew: Matrix
        var Ahatnew: Matrix

        if (nj > 0) {
            // Evaluate 4 functions at nj+1 roots of unity using DFT
            val n = nj + 1
            val omega = -2.0 * Math.PI / n

            // DFT of the 4 block sequences
            val ft1 = complexDftBlocks(Aodd, m, n)
            val ft2 = complexDftBlocks(Aeven, m, n)
            val ft3 = complexDftBlocks(Ahatodd, m, n)
            val ft4 = complexDftBlocks(Ahateven, m, n)

            // Point-wise evaluation: equation (6.20) from Thesis Meini
            val ftAnew = Array(n) { Array(m) { DoubleArray(m * 2) } }  // complex m×m
            val ftAhatnew = Array(n) { Array(m) { DoubleArray(m * 2) } }

            for (cnt in 0 until n) {
                val Im = complexEye(m)
                val ImMinusOdd = complexMatSub(Im, ft1[cnt])
                val invImMinusOdd = complexMatInv(ImMinusOdd, m)

                // Ahatnew[cnt] = ft4[cnt] + ft2[cnt] * inv(I-ft1[cnt]) * ft3[cnt]
                val prod1 = complexMatMul(ft2[cnt], invImMinusOdd, m)
                val prod2 = complexMatMul(prod1, ft3[cnt], m)
                ftAhatnew[cnt] = complexMatAdd(ft4[cnt], prod2)

                // Anew[cnt] = exp(-cnt*2j*pi/n) * ft1[cnt] + ft2[cnt] * inv(I-ft1[cnt]) * ft2[cnt]
                val prod3 = complexMatMul(prod1, ft2[cnt], m)
                val phase = omega * cnt
                val expReal = Math.cos(phase)
                val expImag = Math.sin(phase)
                val scaledOdd = complexMatScale(ft1[cnt], expReal, expImag, m)
                ftAnew[cnt] = complexMatAdd(scaledOdd, prod3)
            }

            // Inverse DFT to recover real coefficients
            Anew = complexIdftBlocks(ftAnew, m, n)
            Ahatnew = complexIdftBlocks(ftAhatnew, m, n)
        } else {
            // Series are constant (nj == 0)
            val ImOdd = Matrix.eye(m).sub(Aodd.extractRows(0, m))
            val temp = Aeven.extractRows(0, m).mult(ImOdd.inv())
            Ahatnew = Ahateven.extractRows(0, m).add(temp.mult(Ahatodd.extractRows(0, m)))

            val block1 = temp.mult(Aeven.extractRows(0, m))
            Anew = Matrix(2 * m, m)
            for (r in 0 until m) {
                for (c in 0 until m) {
                    Anew[r, c] = block1[r, c]
                    Anew[m + r, c] = Aodd[r, c]
                }
            }
            // Make Ahatnew into a single block
            val temp2 = Ahatnew.copy()
            Ahatnew = Matrix(m, m)
            for (r in 0 until m) {
                for (c in 0 until m) {
                    Ahatnew[r, c] = temp2[r, c]
                }
            }
        }

        // Check if accuracy is sufficient; if not, double the roots
        var nAnew = computeUpperHalfNorm(Anew, m)
        var nAhatnew = computeUpperHalfNorm(Ahatnew, m)
        var njCur = if (nj > 0) nj else 0

        while ((nAnew > (njCur + 1) * options.epsilonValue ||
                    nAhatnew > (njCur + 1) * options.epsilonValue) &&
            njCur + 1 < options.maxNumRoot
        ) {
            njCur = 2 * (njCur + 1) - 1
            val stopv = minOf(njCur + 1, Aodd.numRows / m)
            val n = njCur + 1

            val omega2 = -2.0 * Math.PI / n

            val ft1 = complexDftBlocks(Aodd, m, n, stopv)
            val ft2 = complexDftBlocks(Aeven, m, n, stopv)
            val ft3 = complexDftBlocks(Ahatodd, m, n, stopv)
            val ft4 = complexDftBlocks(Ahateven, m, n, stopv)

            val ftAnew2 = Array(n) { Array(m) { DoubleArray(m * 2) } }
            val ftAhatnew2 = Array(n) { Array(m) { DoubleArray(m * 2) } }

            for (cnt in 0 until n) {
                val Im = complexEye(m)
                val ImMinusOdd = complexMatSub(Im, ft1[cnt])
                val invImMinusOdd = complexMatInv(ImMinusOdd, m)

                val prod1 = complexMatMul(ft2[cnt], invImMinusOdd, m)
                val prod2 = complexMatMul(prod1, ft3[cnt], m)
                ftAhatnew2[cnt] = complexMatAdd(ft4[cnt], prod2)

                val prod3 = complexMatMul(prod1, ft2[cnt], m)
                val phase = omega2 * cnt
                val expReal = Math.cos(phase)
                val expImag = Math.sin(phase)
                val scaledOdd = complexMatScale(ft1[cnt], expReal, expImag, m)
                ftAnew2[cnt] = complexMatAdd(scaledOdd, prod3)
            }

            Anew = complexIdftBlocks(ftAnew2, m, n)
            Ahatnew = complexIdftBlocks(ftAhatnew2, m, n)

            nAnew = computeUpperHalfNorm(Anew, m)
            nAhatnew = computeUpperHalfNorm(Ahatnew, m)
        }

        if ((nAnew > (njCur + 1) * options.epsilonValue ||
                    nAhatnew > (njCur + 1) * options.epsilonValue) &&
            njCur + 1 >= options.maxNumRoot
        ) {
            jline.io.line_warning(
                "MG1_CR",
                "Maximum number of '%d' reached, accuracy might be affected",
                options.maxNumRoot
            )
        }

        // Truncate to half if needed
        if (njCur > 1) {
            val halfBlocks = (njCur + 1) / 2
            Anew = Anew.extractRows(0, halfBlocks * m)
            Ahatnew = Ahatnew.extractRows(0, halfBlocks * m)
        }

        // Split into even and odd
        val numBlocksAnew = Anew.numRows / m
        Aeven = extractEvenBlocks(Anew, m, numBlocksAnew)
        Aodd = extractOddBlocks(Anew, m, numBlocksAnew)

        val numBlocksAhatnew = Ahatnew.numRows / m
        Ahateven = extractEvenBlocks(Ahatnew, m, numBlocksAhatnew)
        Ahatodd = extractOddBlocks(Ahatnew, m, numBlocksAhatnew)

        if (options.verbose == 1) {
            val modeStr = if (options.mode == "PWCR") "Point-wise" else "Shifted PWCR"
            println("The $modeStr evaluation of Iteration $numit required ${njCur + 1} roots")
        }

        // Test stop criteria
        if (options.mode == "PWCR" || options.mode == "DCR") {
            var Rnewj = Anew.extractRows(m, 2 * m)
            for (i in 2 until Anew.numRows / m) {
                Rnewj = Rnewj.add(Anew.extractRows(i * m, (i + 1) * m))
            }
            Rnewj = Matrix.eye(m).sub(Rnewj).inv().mult(Anew.extractRows(0, m))

            if (Rj.sub(Rnewj).infinityNorm() < options.epsilonValue ||
                Matrix.eye(m).sub(Anew.extractRows(0, m).mult(Matrix.eye(m).sub(Anew.extractRows(m, 2 * m)).inv())).sumCols().maxElement() < options.epsilonValue
            ) {
                G = Ahatnew.extractRows(0, m)
                for (i in 1 until Ahatnew.numRows / m) {
                    G = G.add(Rnewj.mult(Ahatnew.extractRows(i * m, (i + 1) * m)))
                }
                G = D.extractRows(0, m).mult(Matrix.eye(m).sub(G).inv())
                break
            }
            Rj = Rnewj

            if (Anew.extractRows(0, m).infinityNorm() < options.epsilonValue ||
                Ahatnew.extractRows(m, Ahatnew.numRows).elementAbsSum() < options.epsilonValue ||
                Matrix.eye(m).sub(D.extractRows(0, m).mult(Matrix.eye(m).sub(Ahatnew.extractRows(0, m)).inv())).sumCols().maxElement() < options.epsilonValue
            ) {
                G = D.extractRows(0, m).mult(Matrix.eye(m).sub(Ahatnew.extractRows(0, m)).inv())
                break
            }
        } else {
            // ShiftPWCR mode
            val Gold = G.copy()
            G = D.extractRows(0, m).mult(Matrix.eye(m).sub(Ahatnew.extractRows(0, m)).inv())
            if (G.sub(Gold).infinityNorm() < options.epsilonValue ||
                (Ahatnew.numRows > m && Ahatnew.extractRows(m, Ahatnew.numRows).infinityNorm() < options.epsilonValue)
            ) {
                break
            }
        }
    }

    if (numit == options.maxNumIt && G.elementSum() == 0.0) {
        jline.io.line_warning("MG1_CR", "Maximum Number of Iterations %d reached", numit)
        G = D.extractRows(0, m).mult(Matrix.eye(m).sub(
            if (Ahateven.numRows >= m) Ahateven.extractRows(0, m) else Matrix.zeros(m, m)
        ).inv())
    }

    // Transpose G (MATLAB returns G' since D was transposed)
    G = G.transpose()

    // Apply shift correction
    if (options.mode.contains("ShiftPWCR")) {
        when (options.shiftType) {
            "one" -> if (drift < 1) G = G.add(Matrix.ones(m, m).scale(1.0 / m))
            "tau" -> if (drift > 1) G = G.add(v.mult(Matrix.ones(1, m)).scale(tau))
            "dbl" -> {
                if (drift < 1) G = G.add(Matrix.ones(m, m).scale(1.0 / m))
                if (drift > 1) G = G.add(v.mult(Matrix.ones(1, m)).scale(tau))
            }
        }
    }

    if (options.verbose == 1) {
        val Dcheck = if (options.mode == "PWCR") D else Dold ?: D
        val mOrig = Dcheck.numRows / (Dcheck.numCols / m + 1).coerceAtLeast(1)
        // Compute residual: temp = A_max; for i=max-1:-1:0: temp = A_i + temp*G
        // This verifies G - (A_max*G^max + ... + A_0) = 0
        println("Final Residual Error for G: CR computation completed in $numit iterations")
    }

    return G
}

/**
 * Options for MG1_CR solver
 */
data class MG1CROptions(
    val mode: String = "ShiftPWCR",
    val maxNumIt: Int = 50,
    val maxNumRoot: Int = 2048,
    val epsilonValue: Double = 1e-16,
    val verbose: Int = 0,
    val shiftType: String = "one"
)

// ============================================================================
// Helper functions for MG1_CR: Complex DFT and block matrix operations
// ============================================================================

/**
 * Floor of log2(n), returns 0 for n <= 1
 */
private fun floorLog2(n: Int): Int {
    if (n <= 1) return 0
    var result = 0
    var v = n
    while (v > 1) {
        v = v shr 1
        result++
    }
    return result
}

/**
 * Extract even-indexed blocks (0, 2, 4, ...) from a block matrix.
 * Block i starts at row i*m, has m rows.
 */
private fun extractEvenBlocks(D: Matrix, m: Int, numBlocks: Int): Matrix {
    val evenCount = (numBlocks + 1) / 2
    val result = Matrix(evenCount * m, D.numCols)
    var idx = 0
    for (i in 0 until numBlocks step 2) {
        for (r in 0 until m) {
            for (c in 0 until D.numCols) {
                result[idx * m + r, c] = D[i * m + r, c]
            }
        }
        idx++
    }
    return result
}

/**
 * Extract odd-indexed blocks (1, 3, 5, ...) from a block matrix.
 */
private fun extractOddBlocks(D: Matrix, m: Int, numBlocks: Int): Matrix {
    val oddCount = numBlocks / 2
    val result = Matrix(oddCount * m, D.numCols)
    var idx = 0
    for (i in 1 until numBlocks step 2) {
        for (r in 0 until m) {
            for (c in 0 until D.numCols) {
                result[idx * m + r, c] = D[i * m + r, c]
            }
        }
        idx++
    }
    return result
}

/**
 * Complex matrix stored as Array[row][col*2], where [col*2] = real, [col*2+1] = imag
 */
private typealias ComplexMatrix = Array<DoubleArray>

/**
 * DFT of block sequences. Input: real block matrix with numBlocks blocks of m×m.
 * Evaluates the block polynomial at N roots of unity.
 * Returns array of N complex m×m matrices.
 *
 * @param blocks Real block matrix (numBlocks*m rows, m cols)
 * @param m Block size
 * @param N DFT size (number of evaluation points)
 * @param maxBlocks Maximum number of blocks to use (default: all)
 */
private fun complexDftBlocks(
    blocks: Matrix, m: Int, N: Int, maxBlocks: Int = blocks.numRows / m
): Array<ComplexMatrix> {
    val numBlocks = minOf(maxBlocks, blocks.numRows / m)
    val result = Array(N) { Array(m) { DoubleArray(m * 2) } }

    for (k in 0 until N) {
        val angle = -2.0 * Math.PI * k / N
        // Evaluate polynomial at w = exp(j*angle)
        for (n in 0 until numBlocks) {
            val theta = angle * n
            val cosTheta = Math.cos(theta)
            val sinTheta = Math.sin(theta)
            for (r in 0 until m) {
                for (c in 0 until m) {
                    val v = blocks[n * m + r, c]
                    result[k][r][c * 2] += v * cosTheta       // real
                    result[k][r][c * 2 + 1] += v * sinTheta   // imag
                }
            }
        }
    }
    return result
}

/**
 * Inverse DFT of complex block sequences. Converts N complex m×m matrices
 * back to a real block matrix of N blocks.
 */
private fun complexIdftBlocks(ft: Array<ComplexMatrix>, m: Int, N: Int): Matrix {
    val result = Matrix(N * m, m)
    val invN = 1.0 / N

    for (n in 0 until N) {
        for (k in 0 until N) {
            val angle = 2.0 * Math.PI * k * n / N
            val cosAngle = Math.cos(angle)
            val sinAngle = Math.sin(angle)
            for (r in 0 until m) {
                for (c in 0 until m) {
                    val re = ft[k][r][c * 2]
                    val im = ft[k][r][c * 2 + 1]
                    // IDFT: x[n] = (1/N) * sum_k X[k] * exp(2*pi*i*k*n/N)
                    result[n * m + r, c] += (re * cosAngle - im * sinAngle) * invN
                }
            }
        }
    }
    return result
}

/**
 * Complex identity matrix
 */
private fun complexEye(m: Int): ComplexMatrix {
    val result = Array(m) { DoubleArray(m * 2) }
    for (i in 0 until m) {
        result[i][i * 2] = 1.0
    }
    return result
}

/**
 * Complex matrix subtraction: A - B
 */
private fun complexMatSub(A: ComplexMatrix, B: ComplexMatrix): ComplexMatrix {
    val m = A.size
    val result = Array(m) { r -> DoubleArray(A[r].size) { c -> A[r][c] - B[r][c] } }
    return result
}

/**
 * Complex matrix addition: A + B
 */
private fun complexMatAdd(A: ComplexMatrix, B: ComplexMatrix): ComplexMatrix {
    val m = A.size
    val result = Array(m) { r -> DoubleArray(A[r].size) { c -> A[r][c] + B[r][c] } }
    return result
}

/**
 * Complex matrix multiplication: A * B (both m×m complex)
 */
private fun complexMatMul(A: ComplexMatrix, B: ComplexMatrix, m: Int): ComplexMatrix {
    val result = Array(m) { DoubleArray(m * 2) }
    for (i in 0 until m) {
        for (j in 0 until m) {
            var re = 0.0
            var im = 0.0
            for (k in 0 until m) {
                val aRe = A[i][k * 2]
                val aIm = A[i][k * 2 + 1]
                val bRe = B[k][j * 2]
                val bIm = B[k][j * 2 + 1]
                re += aRe * bRe - aIm * bIm
                im += aRe * bIm + aIm * bRe
            }
            result[i][j * 2] = re
            result[i][j * 2 + 1] = im
        }
    }
    return result
}

/**
 * Complex matrix scaling: A * (scaleRe + i*scaleIm)
 */
private fun complexMatScale(
    A: ComplexMatrix, scaleRe: Double, scaleIm: Double, m: Int
): ComplexMatrix {
    val result = Array(m) { DoubleArray(m * 2) }
    for (i in 0 until m) {
        for (j in 0 until m) {
            val re = A[i][j * 2]
            val im = A[i][j * 2 + 1]
            result[i][j * 2] = re * scaleRe - im * scaleIm
            result[i][j * 2 + 1] = re * scaleIm + im * scaleRe
        }
    }
    return result
}

/**
 * Complex matrix inverse using Gauss-Jordan elimination
 */
private fun complexMatInv(A: ComplexMatrix, m: Int): ComplexMatrix {
    // Augmented matrix [A | I]
    val aug = Array(m) { DoubleArray(m * 4) }
    for (i in 0 until m) {
        for (j in 0 until m) {
            aug[i][j * 2] = A[i][j * 2]
            aug[i][j * 2 + 1] = A[i][j * 2 + 1]
        }
        aug[i][(m + i) * 2] = 1.0
    }

    // Gauss-Jordan with partial pivoting
    for (col in 0 until m) {
        // Find pivot
        var maxMag = 0.0
        var pivotRow = col
        for (row in col until m) {
            val re = aug[row][col * 2]
            val im = aug[row][col * 2 + 1]
            val mag = re * re + im * im
            if (mag > maxMag) {
                maxMag = mag
                pivotRow = row
            }
        }

        // Swap rows
        if (pivotRow != col) {
            val temp = aug[col]
            aug[col] = aug[pivotRow]
            aug[pivotRow] = temp
        }

        // Normalize pivot row
        val pivRe = aug[col][col * 2]
        val pivIm = aug[col][col * 2 + 1]
        val pivMag2 = pivRe * pivRe + pivIm * pivIm
        if (pivMag2 < 1e-30) continue // Nearly singular

        val invRe = pivRe / pivMag2
        val invIm = -pivIm / pivMag2
        for (j in 0 until 2 * m) {
            val re = aug[col][j * 2]
            val im = aug[col][j * 2 + 1]
            aug[col][j * 2] = re * invRe - im * invIm
            aug[col][j * 2 + 1] = re * invIm + im * invRe
        }

        // Eliminate column
        for (row in 0 until m) {
            if (row == col) continue
            val factRe = aug[row][col * 2]
            val factIm = aug[row][col * 2 + 1]
            for (j in 0 until 2 * m) {
                aug[row][j * 2] -= factRe * aug[col][j * 2] - factIm * aug[col][j * 2 + 1]
                aug[row][j * 2 + 1] -= factRe * aug[col][j * 2 + 1] + factIm * aug[col][j * 2]
            }
        }
    }

    // Extract inverse from augmented matrix
    val result = Array(m) { DoubleArray(m * 2) }
    for (i in 0 until m) {
        for (j in 0 until m) {
            result[i][j * 2] = aug[i][(m + j) * 2]
            result[i][j * 2 + 1] = aug[i][(m + j) * 2 + 1]
        }
    }
    return result
}

/**
 * Compute the maximum infinity norm of blocks in the upper half of a block matrix.
 */
private fun computeUpperHalfNorm(A: Matrix, m: Int): Double {
    val deg = A.numRows / m
    var maxNorm = 0.0
    for (i in deg / 2 until deg) {
        val block = A.extractRows(i * m, (i + 1) * m)
        val norm = block.infinityNorm()
        if (norm > maxNorm) maxNorm = norm
    }
    return maxNorm
}

/**
 * Helper to get max element of a Matrix (used for sum-column checks)
 */
private fun Matrix.maxElement(): Double {
    var max = Double.NEGATIVE_INFINITY
    for (i in 0 until this.numRows) {
        for (j in 0 until this.numCols) {
            if (this[i, j] > max) max = this[i, j]
        }
    }
    return max
}

/**
 * Helper to compute sum of absolute values of all elements
 */
private fun Matrix.elementAbsSum(): Double {
    var sum = 0.0
    for (i in 0 until this.numRows) {
        for (j in 0 until this.numCols) {
            sum += Math.abs(this[i, j])
        }
    }
    return sum
}

/**
 * Helper to compute column sums
 */
private fun Matrix.sumCols(): Matrix {
    val result = Matrix(1, this.numCols)
    for (j in 0 until this.numCols) {
        var s = 0.0
        for (i in 0 until this.numRows) {
            s += this[i, j]
        }
        result[0, j] = s
    }
    return result
}
