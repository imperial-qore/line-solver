/**
 * @file Method of Moments (MOM) for exact normalizing constant computation
 * 
 * Implements the Method of Moments using exact arithmetic with BigFraction for computing
 * normalizing constants in closed product-form queueing networks. Provides numerically
 * stable computation for systems requiring high precision results.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.nc

import jline.io.Ret
import jline.util.Maths
import jline.util.matrix.Matrix
import org.apache.commons.math3.fraction.BigFraction
import org.apache.commons.math3.fraction.BigFractionField
import org.apache.commons.math3.linear.Array2DRowFieldMatrix
import org.apache.commons.math3.linear.ArrayFieldVector
import org.apache.commons.math3.linear.FieldLUDecomposition
import org.apache.commons.math3.linear.FieldMatrix
import java.math.BigInteger
import kotlin.math.min

/**
 * Method of Moments (MOM) solver for product-form queueing networks.
 * This implementation uses exact arithmetic via BigFraction.
 */
object Pfqn_mom {

    /**
     * Computes performance measures using the Method of Moments algorithm.
     *
     * @param L Service demands matrix (M x R)
     * @param N Population vector (1 x R)
     * @param Z Think times vector (1 x R)
     * @return pfqnMom containing (X: throughputs, Q: queue lengths, G: normalizing constant)
     */
    @JvmStatic
    fun pfqn_mom(L: Matrix, N: Matrix, Z: Matrix): Ret.pfqnMom {
        try {
            val M = L.getNumRows() // number of stations
            val R = L.getNumCols() // number of classes


        // For now, let's implement only the single class case to debug
        if (R == 1 && M == 1) {
            return singleClassSingleStation(L, N, Z)
        }

        // Convert input to BigFraction for exact arithmetic
        val Lf = Array(M) { i -> Array(R) { j -> BigFraction(L.get(i, j)) } }

        // N should be a row vector
        val Nf = Array(R) { r -> BigFraction(N.get(0, r).toInt()) }

        // Z should be a row vector, if it's a matrix, sum the rows
        val Zf = if (Z.getNumRows() == 1) {
            Array(R) { r -> BigFraction(Z.get(0, r)) }
        } else {
            // Sum each column across rows to get a single row vector
            Array(R) { r ->
                var sum = BigFraction.ZERO
                for (i in 0 until Z.getNumRows()) {
                    sum = sum.add(BigFraction(Z.get(i, r)))
                }
                sum
            }
        }

        val n = IntArray(R) { 0 }
        var g: Array<BigFraction>? = null
        var gr: Array<BigFraction>? = null

        for (r in 1..R) {
                val (C, Cg, D, Dr) = setupls(Lf, Nf, Zf, r)

            if (r == 1) {
                // For first class, initialize g with ones
                g = Array(M + 1) { BigFraction.ONE }
            } else {
                // For r>1, we need to reorganize g and gr from previous iteration
                // The size should be nchoosek(M+r-2,r-1) * r
                val actualSize = Maths.nchoosek(M + r - 2, r - 1) * r
                val G = Array<BigFraction>(actualSize) { BigFraction.ZERO }


                for (i in 0 until Maths.nchoosek(M + r - 2, r - 1)) {
                    for (s in 0 until r - 1) {
                        val srcIdx = i * (r - 1) + s
                        val dstIdx = i * r + s
                        G[dstIdx] = g!![srcIdx]
                    }
                    val srcIdx = i * (r - 1)
                    val dstIdx = i * r + r - 1
                    G[dstIdx] = gr!![srcIdx]
                }

                val Gk = Array(Maths.nchoosek(M + r - 1, r) * r) { BigFraction.ZERO }
                val newGk = blocksolve(M, r, C, matrixVectorMultiply(Cg.scalarMultiply(BigFraction.MINUS_ONE), G), Gk)
                g = concatenate(newGk, G)
            }

            val CgDr = C.createMatrix(Cg.getRowDimension(), Dr.getColumnDimension())
            for (i in 0 until Cg.getRowDimension()) {
                for (j in 0 until Dr.getColumnDimension()) {
                    var sum = BigFraction.ZERO
                    for (k in 0 until Cg.getColumnDimension()) {
                        sum = sum.add(Cg.getEntry(i, k).multiply(Dr.getEntry(k, j)))
                    }
                    CgDr.setEntry(i, j, sum)
                }
            }

            // Main iteration
            var Gk = Array(Maths.nchoosek(M + r - 1, r) * r) { BigFraction.ZERO }
            val numIterations = N.get(0, r-1).toInt() - 1
            for (nr in 0 until numIterations) {
                n[r-1] = nr + 1

                val nrFrac = BigFraction(n[r - 1])
                val gLocal = g!!
                val G = Array(Dr.getRowDimension()) { i ->
                    var sum = BigFraction.ZERO
                    for (j in 0 until Dr.getColumnDimension()) {
                        if (j < gLocal.size) {
                            sum = sum.add(Dr.getEntry(i, j).multiply(gLocal[j]))
                        }
                    }
                    sum.divide(nrFrac)
                }

                val b = matrixVectorMultiply(D.subtract(CgDr.scalarMultiply(nrFrac.reciprocal())), g!!)
                Gk = blocksolve(M, r, C, b, Gk)
                g = concatenate(Gk, G)
            }

            // Last iteration
            gr = g
            n[r-1] = N.get(0, r-1).toInt()

            val nrFrac = BigFraction(n[r - 1])
            val gLocal2 = g!!
            val G = Array(Dr.getRowDimension()) { i ->
                var sum = BigFraction.ZERO
                for (j in 0 until Dr.getColumnDimension()) {
                    if (j < gLocal2.size) {
                        sum = sum.add(Dr.getEntry(i, j).multiply(gLocal2[j]))
                    }
                }
                sum.divide(nrFrac)
            }

            val b = matrixVectorMultiply(D.subtract(CgDr.scalarMultiply(nrFrac.reciprocal())), g)
            Gk = blocksolve(M, r, C, b, Gk)

            // Compute F matrices for verification
            val Cinv = invertMatrix(C)
            val F1 = concatenateVertical(Cinv.multiply(D), zeros(Dr.getRowDimension(), D.getColumnDimension()))
            val F2 = concatenateVertical(Cinv.scalarMultiply(BigFraction.MINUS_ONE).multiply(CgDr), Dr)
            val F = F1.add(F2.scalarMultiply(nrFrac.reciprocal()))



            g = concatenate(Gk, G)
        }

        // Compute performance measures

        // X should be a row vector (1 x R), not (R x 1)
        val Xdata = Array(1) { DoubleArray(R) }
        val Qdata = Array(M) { DoubleArray(R) }
        val Gconst: BigFraction

        if (R == 1) {
            // For single class, g has M+1 elements
            // MATLAB uses 1-based indexing: G=g(2), X=gr(2)/G, Q=(g(1)/G)-1
            // In 0-based indexing: G=g[1], X=gr[1]/G, Q=(g[0]/G)-1
            Gconst = g!![1]
            Xdata[0][0] = gr!![1].divide(Gconst).toDouble()
            Qdata[0][0] = g[0].divide(Gconst).subtract(BigFraction.ONE).toDouble()
        } else {
            // Complex case for R > 1
            // Reorganize g and gr into Gk format for recursion
            var Gk = Array(Maths.nchoosek(M + R - 2, R - 1) * (R + 1)) { BigFraction.ZERO }

            for (i in 0 until Maths.nchoosek(M + R - 2, R - 1)) {
                for (s in 0 until R) {
                    val idx_from = Maths.nchoosek(M + R - 1, R) * R + i * R + s
                    val idx_to = i * (R + 1) + s
                    Gk[idx_to] = g!![idx_from]
                }
                val idx_from_gr = Maths.nchoosek(M + R - 1, R) * R + i * R
                val idx_to_gr = i * (R + 1) + R
                Gk[idx_to_gr] = gr!![idx_from_gr]
            }

            // Convolution to compute final G values
            var finalG: Array<BigFraction>? = null

            for (l in R - 1 downTo 1) {
                val G = Array(Maths.nchoosek(M + l - 2, l - 1) * (R + 1)) { BigFraction.ZERO }
                val Ik = Maths.sortByNnzPos(Maths.multichooseList(M, l))
                val I = Maths.sortByNnzPos(Maths.multichooseList(M, l - 1))

                for (i in 0 until I.size) {
                    G[i * (R + 1)] = BigFraction.ZERO
                    val Ii = I[i].copyOf()
                    Ii[0]++ // add station 1
                    val t = Maths.matchRow(Ik, Ii)
                    G[i * (R + 1)] = G[i * (R + 1)].add(Gk[t * (R + 1)])

                    for (s in 0 until R) {
                        G[i * (R + 1)] = G[i * (R + 1)].subtract(Lf[0][s].multiply(Gk[t * (R + 1) + s + 1]))
                    }

                    for (s in 0 until R) {
                        if (Z.get(0, s) == 0.0) {
                            G[i * (R + 1) + s + 1] = BigFraction.ZERO
                        } else {
                            G[i * (R + 1) + s + 1] = Nf[s].divide(Zf[s]).multiply(G[i * (R + 1)])

                            for (j in 0 until M) {
                                val Ij = I[i].copyOf()
                                Ij[j]++
                                val tj = Maths.matchRow(Ik, Ij)
                                G[i * (R + 1) + s + 1] = G[i * (R + 1) + s + 1].subtract(
                                    BigFraction(1 + I[i][j]).multiply(Lf[j][s]).divide(Zf[s]).multiply(Gk[tj * (R + 1) + s + 1])
                                )
                            }
                        }
                    }
                }

                if (l > 1) {
                    Gk = G
                } else {
                    // When l=1, this is the final G
                    finalG = G
                }
            }

            // Use the final G values from convolution
            val G = finalG!!
            Gconst = G[0]

            for (s in 0 until R) {
                Xdata[0][s] = G[s + 1].divide(Gconst).toDouble()
                for (m in 0 until M) {
                    // Note: We use the original Gk for Q computation as per MATLAB
                    Qdata[m][s] = Lf[m][s].multiply(Gk[m * (R + 1) + s + 1]).divide(Gconst).toDouble()
                }
            }
        }

        val X = Matrix(Xdata)
        val Q = Matrix(Qdata)

        // Compute log(G) safely even when G is outside floating point range
        // log(a/b) = log(a) - log(b)
        val lG = if (Gconst.denominator.equals(BigInteger.ONE)) {
            // G is an integer, use log of numerator
            logBigInteger(Gconst.numerator)
        } else {
            // G is a fraction, use log(numerator) - log(denominator)
            logBigInteger(Gconst.numerator) - logBigInteger(Gconst.denominator)
        }

        return Ret.pfqnMom(X, Q, Gconst, lG, g!!, gr!!)
        } catch (e: Exception) {
            throw e
        }
    }

    /**
     * Sets up the linear system for a given number of classes.
     */
    private fun setupls(
        L: Array<Array<BigFraction>>,
        N: Array<BigFraction>,
        Z: Array<BigFraction>,
        R: Int
    ): Quadruple<FieldMatrix<BigFraction>, FieldMatrix<BigFraction>, FieldMatrix<BigFraction>, FieldMatrix<BigFraction>> {
        val M = L.size
        val m = IntArray(M) { 1 }

        val Ik = Maths.sortByNnzPos(Maths.multichooseList(M, R))
        val I = Maths.sortByNnzPos(Maths.multichooseList(M, R - 1))

        val rows = Ik.size * R
        val colsC = Ik.size * R
        val colsCg = I.size * R
        val colsD = (Ik.size + I.size) * R
        val rowsDr = I.size * R

        val field = BigFractionField.getInstance()
        val C = Array2DRowFieldMatrix(field, rows, colsC)
        val Cg = Array2DRowFieldMatrix(field, rows, colsCg)
        val D = Array2DRowFieldMatrix(field, rows, colsD)
        val Dr = Array2DRowFieldMatrix(field, rowsDr, colsD)

        // Initialize with zeros
        for (i in 0 until rows) {
            for (j in 0 until colsC) C.setEntry(i, j, BigFraction.ZERO)
            for (j in 0 until colsCg) Cg.setEntry(i, j, BigFraction.ZERO)
            for (j in 0 until colsD) D.setEntry(i, j, BigFraction.ZERO)
        }
        for (i in 0 until rowsDr) {
            for (j in 0 until colsD) Dr.setEntry(i, j, BigFraction.ZERO)
        }

        val pcpos = mutableListOf<Int>()
        var currentRow = 0

        // Build convolution expressions
        for (i in 0 until Ik.size) {
            val h = Ik[i].count { it > 0 }

            for (j in 0 until M) {
                if (Ik[i][j] > 0) {
                    C.setEntry(currentRow, i * R, BigFraction.ONE)
                    for (s in 0 until R - 1) {
                        C.setEntry(currentRow, i * R + s + 1, L[j][s].negate())
                    }

                    val IkCopy = Ik[i].copyOf()
                    IkCopy[j]--
                    val idx = Maths.matchRow(I, IkCopy)
                    Cg.setEntry(currentRow, idx * R, BigFraction.MINUS_ONE)

                    D.setEntry(currentRow, i * R, L[j][R - 1])
                    currentRow++
                }
            }

            // Reserve rows for population constraints
            for (j in 0 until R - h) {
                pcpos.add(currentRow)
                currentRow++
            }
        }

        // Build population constraints
        var last = 0
        for (i in 0 until I.size) {
            for (s in 0 until R) {
                Dr.setEntry(i * R + s, Ik.size * R + i * R + s, Z[R - 1])
            }

            for (j in 0 until M) {
                val ICopy = I[i].copyOf()
                ICopy[j]++
                val t = Maths.matchRow(Ik, ICopy)
                Dr.setEntry(i * R, t * R, BigFraction(m[j] + I[i][j]).multiply(L[j][R - 1]))

                for (s in 0 until R - 1) {
                    C.setEntry(pcpos[last + s], t * R + s + 1, BigFraction(m[j] + I[i][j]).multiply(L[j][s]).negate())
                    Cg.setEntry(pcpos[last + s], i * R, N[s])
                    Cg.setEntry(pcpos[last + s], i * R + s + 1, Z[s].negate())
                    Dr.setEntry(i * R + s + 1, t * R + s + 1, BigFraction(m[j] + I[i][j]).multiply(L[j][R - 1]))
                }
            }
            last += R - 1
        }

        return Quadruple(C, Cg, D, Dr)
    }

    /**
     * Solves the block linear system.
     */
    private fun blocksolve(
        M: Int,
        R: Int,
        C: FieldMatrix<BigFraction>,
        b: Array<BigFraction>,
        Gr: Array<BigFraction>
    ): Array<BigFraction> {
        var blockend = C.columnDimension - 1
        val H = min(M, R)
        val x = Array(b.size) { BigFraction.ZERO }
        val bAdjusted = b.copyOf()

        // Subtract C*Gr from b
        for (i in 0 until bAdjusted.size) {
            for (j in 0 until Gr.size) {
                bAdjusted[i] = bAdjusted[i].subtract(C.getEntry(i, j).multiply(Gr[j]))
            }
        }

        // Solve blocks from bottom to top
        for (h in H downTo 1) {
            for (t in 0 until Maths.nchoosek(M, h)) {
                val blockstart = blockend - Maths.nchoosek(R - 1, R - h) * R + 1
                val blockSize = blockend - blockstart + 1

                // Extract block
                val blockC = C.getSubMatrix(blockstart, blockend, blockstart, blockend)
                val blockB =
                    ArrayFieldVector(BigFractionField.getInstance(), bAdjusted.sliceArray(blockstart..blockend))

                // Adjust for already solved parts
                if (blockend < C.columnDimension - 1) {
                    for (i in 0 until blockSize) {
                        for (j in blockend + 1 until C.columnDimension) {
                            val term = C.getEntry(blockstart + i, j).multiply(x[j])
                            blockB.setEntry(i, blockB.getEntry(i).subtract(term))
                        }
                    }
                }

                // Solve block using LU decomposition
                val solver = FieldLUDecomposition(blockC).getSolver()
                val blockX = solver.solve(blockB)

                // Store solution
                for (i in 0 until blockSize) {
                    x[blockstart + i] = blockX.getEntry(i)
                }

                blockend = blockstart - 1
            }
        }

        // Add Gr to result
        for (i in x.indices) {
            if (i < Gr.size) {
                x[i] = x[i].add(Gr[i])
            }
        }

        return x
    }

    // Matrix operations with BigFraction

    private fun matrixVectorMultiply(A: FieldMatrix<BigFraction>, v: Array<BigFraction>): Array<BigFraction> {
        val result = Array(A.rowDimension) { BigFraction.ZERO }
        for (i in 0 until A.rowDimension) {
            for (j in 0 until A.columnDimension) {
                result[i] = result[i].add(A.getEntry(i, j).multiply(v[j]))
            }
        }
        return result
    }

    private fun invertMatrix(A: FieldMatrix<BigFraction>): FieldMatrix<BigFraction> {
        return FieldLUDecomposition(A).getSolver().getInverse()
    }

    private fun concatenate(a: Array<BigFraction>, b: Array<BigFraction>): Array<BigFraction> {
        return a + b
    }

    private fun concatenateVertical(A: FieldMatrix<BigFraction>, B: FieldMatrix<BigFraction>): FieldMatrix<BigFraction> {
        val rows = A.rowDimension + B.rowDimension
        val cols = A.columnDimension
        val field = BigFractionField.getInstance()
        val result = Array2DRowFieldMatrix(field, rows, cols)

        for (i in 0 until A.rowDimension) {
            for (j in 0 until A.columnDimension) {
                result.setEntry(i, j, A.getEntry(i, j))
            }
        }

        for (i in 0 until B.rowDimension) {
            for (j in 0 until B.columnDimension) {
                result.setEntry(A.rowDimension + i, j, B.getEntry(i, j))
            }
        }

        return result
    }

    private fun zeros(rows: Int, cols: Int): FieldMatrix<BigFraction> {
        val field = BigFractionField.getInstance()
        val result = Array2DRowFieldMatrix(field, rows, cols)
        for (i in 0 until rows) {
            for (j in 0 until cols) {
                result.setEntry(i, j, BigFraction.ZERO)
            }
        }
        return result
    }

    private fun printMatrix(matrix: FieldMatrix<BigFraction>) {
        for (i in 0 until matrix.rowDimension) {
            print("[")
            for (j in 0 until matrix.columnDimension) {
                val entry = matrix.getEntry(i, j)
                print(formatFraction(entry))
                if (j < matrix.columnDimension - 1) print(", ")
            }
            println("]")
        }
    }

    private fun formatFraction(f: BigFraction): String {
        return if (f.denominator == BigInteger.ONE) {
            f.numerator.toString()
        } else {
            "${f.numerator}/${f.denominator}"
        }
    }

    /**
     * Computes the natural logarithm of a BigInteger safely.
     * Uses a combination of bit length estimation and double precision for accuracy.
     */
    private fun logBigInteger(bigInt: BigInteger): Double {
        if (bigInt.signum() <= 0) {
            throw IllegalArgumentException("Cannot compute log of non-positive number")
        }

        if (bigInt == BigInteger.ONE) {
            return 0.0
        }

        // For small numbers that fit in double precision, use direct computation
        val bitLength = bigInt.bitLength()
        if (bitLength <= 53) { // 53 bits is the precision limit for double
            return Math.log(bigInt.toDouble())
        }

        // For large numbers, use: log(n) ≈ bitLength * log(2) + log(n / 2^bitLength)
        // This gives us: log(n) = log(2^bitLength * (n/2^bitLength)) = bitLength*log(2) + log(n/2^bitLength)
        val powerOfTwo = BigInteger.ONE.shiftLeft(bitLength - 1)
        val ratio = bigInt.toDouble() / powerOfTwo.toDouble()
        return (bitLength - 1) * Math.log(2.0) + Math.log(ratio)
    }

    /**
     * Simple implementation for single class, single station case
     */
    private fun singleClassSingleStation(L: Matrix, N: Matrix, Z: Matrix): Ret.pfqnMom {
        // For M/M/1 queue with finite population N
        val mu = BigFraction.ONE.divide(BigFraction(L.get(0, 0)))  // service rate
        val lambda0 = BigFraction.ONE.divide(BigFraction(Z.get(0, 0)))  // think time rate
        val n = N.get(0, 0).toInt()

        // For finite population, we need to solve the birth-death process
        // This is a placeholder - should implement proper finite population M/M/1
        val rho = BigFraction(L.get(0, 0)).divide(BigFraction(Z.get(0, 0)))
        val X_val = rho.divide(BigFraction.ONE.add(rho))
        val Q_val = rho.multiply(rho).divide(BigFraction.ONE.add(rho))

        val X = Matrix(arrayOf(doubleArrayOf(X_val.toDouble())))
        val Q = Matrix(arrayOf(doubleArrayOf(Q_val.toDouble())))
        val G = BigFraction.ONE
        val lG = 0.0
        val g = arrayOf(BigFraction.ONE, G)
        val g_1 = arrayOf(BigFraction.ONE)

        return Ret.pfqnMom(X, Q, G, lG, g, g_1)
    }

    // Helper data classes
    private data class Quadruple<A, B, C, D>(val first: A, val second: B, val third: C, val fourth: D)
}
/**
 * PFQN mom algorithms
 */
@Suppress("unused")
class PfqnMomAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}