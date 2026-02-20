package jline.lib.smc

import jline.io.line_warning
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath
import kotlin.math.abs
import kotlin.math.min
import kotlin.math.sign

/**
 * Invariant Subspace for Quasi-Birth-Death Markov Chains [Akar, Sohraby]
 *
 * DISCRETE TIME CASE:
 * Computes the minimal nonnegative solution to the matrix equation G = A0 + A1 G + A2 G^2,
 * where A0, A1 and A2 are square nonnegative matrices, with (A0+A1+A2) irreducible and stochastic
 *
 * CONTINUOUS TIME CASE:
 * Computes the minimal nonnegative solution to the matrix equation 0 = A0 + A1 G + A2 G^2,
 * where A0, A1 and A2 are square nonnegative matrices, with (A0+A1+A2) having row sums equal to zero
 *
 * @param A0 transition matrix A0
 * @param A1 transition matrix A1
 * @param A2 transition matrix A2
 * @param MaxNumIt_ maximum number of iterations (default: 50)
 * @param Verbose_ verbose output flag (default: 0)
 * @param Mode_ solution mode: 'MSignStandard', 'MSignBalzer', 'Schur' (default: 'Schur')
 * @param RAPComp_ RAP computation flag (default: 0)
 * @return map containing matrices G, R, and U
 */
fun QBD_IS(A0: Matrix,
           A1: Matrix,
           A2: Matrix,
           MaxNumIt_: Int?,
           Verbose_: Int?,
           Mode_: String?,
           RAPComp_: Int?): Map<String?, Matrix?> {
    var A1 = A1
    var Mode = "Schur"
    var MaxNumIt = 50
    var Verbose = false
    var RAPComp = false
    val m = A1.numRows

    if (MaxNumIt_ != null) {
        MaxNumIt = MaxNumIt_
    }

    if (Mode_ != null) {
        if (Mode_ == "MSignStandard" || Mode_ == "MSignBalzer" || Mode_ == "Schur") {
            Mode = Mode_
        } else {
            throw RuntimeException("QBD_LR mode not recognized")
        }
    }

    if (Verbose_ != null) {
        if (Verbose_ == 1) {
            Verbose = true
        }
    }

    if (RAPComp_ != null) {
        if (RAPComp_ == 1) {
            RAPComp = true
        }
    }
    val A1_diag = Matrix(0, 0, 0)
    Matrix.extractDiag(A1, A1_diag)
    var continues = true
    val lamb = Matrix.negative(A1_diag).elementMax()

    if (!RAPComp) {
        continues = false
        if (A1_diag.elementSum() < 0) {
            continues = true
            A0.scaleEq(1 / lamb)
            A1.scaleEq(1 / lamb)
            A1 = A1.add(1.0, Matrix.eye(m))
            A2.scaleEq(1 / lamb)
        }
        QBD_ParsePara(A0, A1, A2)
    } else {
        QBD_ParsePara(A0, A1, A2)
        A0.scaleEq(1 / lamb)
        A1.scaleEq(1 / lamb)
        A1 = A1.add(1.0, Matrix.eye(m))
        A2.scaleEq(1 / lamb)
    }
    var result = QBD_EG(A0, A1, A2, Verbose)
    val GCheck = result["G"]
    if (GCheck != null && GCheck.length() > 0) {
        return result
    }

    val epsilon = FastMath.pow(10.0, -12)
    val f = 2

    val theta = stat(A0.add(1.0, A1).add(1.0, A2))
    val drift = theta.mult(A0.sumRows())[0] - theta.mult(A2.sumRows())[0]
    val F: MutableMap<Int, Matrix> = HashMap()
    F[1] = Matrix.scaleMult(A0, -1.0)
    F[2] = Matrix.eye(m).add(-1.0, A1)
    F[3] = Matrix.scaleMult(A2, -1.0)
    val H: MutableMap<Int, Matrix> = HashMap()
    for (i in 0..f) {
        H[i + 1] = Matrix(m, m)
    }
    for (i in 0..f) {
        var con1 = doubleArrayOf(1.0)
        var con2 = doubleArrayOf(1.0)
        val temp1 = doubleArrayOf(1.0, -1.0)
        val temp2 = doubleArrayOf(1.0, 1.0)

        // Convolution for (1-s)^(f-i)
        for (j in 1..(f - i)) {
            con1 = convolution(con1, temp1)
        }

        // Convolution for (1+s)^i
        for (j in 1..i) {
            con2 = convolution(con2, temp2)
        }

        val contrib = convolution(con1, con2)
        for (j in 0..f) {
            if (j < contrib.size) {
                H[j + 1] = H[j + 1]!!.add(contrib[j], F[i + 1]!!)
            }
        }
    }

    // Step 3: \hat{H}_i = H_f^-1 * H_i
    val HfInv = H[f + 1]!!.inv()
    val hatH: MutableMap<Int, Matrix> = HashMap()
    for (i in 0..(f - 1)) {
        hatH[i + 1] = HfInv.mult(H[i + 1]!!)
    }

    // Step 4: y, xT
    val y = Matrix(m * f, 1)
    for (i in 0..<m) {
        y[i, 0] = 1.0
    }

    val tempMatrix = Matrix(m, m + 1)
    for (i in 0..<m) {
        for (j in 0..<m) {
            tempMatrix[i, j] = hatH[1]!![i, j]
        }
        tempMatrix[i, m] = 1.0
    }
    val x0T = Matrix(1, m + 1)
    for (i in 0..<m) {
        x0T[0, i] = 0.0
    }
    x0T[0, m] = 1.0
    val x0TSolved = x0T.mult(tempMatrix.inv())

    val xT = Matrix(1, m * f)
    for (i in 1..(f - 1)) {
        val hatHi = hatH[i]!!
        for (j in 0..<m) {
            var sum = 0.0
            for (k in 0..<m) {
                sum += x0TSolved[0, k] * hatHi[k, j]
            }
            xT[0, (i - 1) * m + j] = sum
        }
    }
    for (i in 0..<m) {
        xT[0, (f - 1) * m + i] = x0TSolved[0, i]
    }

    // Step 5: E_m in Zold
    val Zold = Matrix(m * f, m * f)
    for (i in 1..(f - 1)) {
        for (j in 0..<m) {
            Zold[(i - 1) * m + j, i * m + j] = 1.0
        }
    }
    for (i in 0..(f - 1)) {
        val hatHi = hatH[i + 1]!!
        for (row in 0..<m) {
            for (col in 0..<m) {
                Zold[m * (f - 1) + row, i * m + col] = -hatHi[row, col]
            }
        }
    }

    val yNorm = y.copy()
    val xTy = xT.mult(y)[0, 0]
    yNorm.scaleEq(1.0 / xTy)

    val yxT = yNorm.mult(xT)
    val signDrift = sign(drift)
    Zold.addEq(-signDrift, yxT)

    // Step 6: Matrix sign function algorithm or Schur decomposition
    var T: Matrix
    if (Mode != "Schur") {
        var numit = 0
        var check = 1.0
        var Znew = Zold.copy()

        while (check > epsilon && numit < MaxNumIt) {
            numit++
            val determ = if (Mode == "MSignStandard") {
                0.5
            } else { // MSignBalzer
                val detVal = abs(Znew.det())
                val detTerm = FastMath.pow(detVal, 1.0 / (m * f))
                min(1.0 / (1.0 + detTerm), 1.0 - 1e-3)
            }

            val ZnewNext = Matrix.scaleMult(Znew, determ).add(1.0 - determ, Znew.inv())
            check = Matrix.firstNorm(ZnewNext.add(-1.0, Znew)) / Matrix.firstNorm(Znew)

            Znew = ZnewNext
        }

        if (numit == MaxNumIt && check > epsilon) {
            line_warning("QBD_IS", "Maximum Number of Iterations %d reached: T may not have m columns", numit)
        }

        // Step 7: Orthogonal basis
        T = orthogonalBasis(Znew.add(-1.0, Matrix.eye(m * f)))
    } else {
        // Schur decomposition approach
        val (TSchur, _) = schurDecomposition(Zold)
        T = Matrix.extractColumns(TSchur, 0, m)
    }

    // Step 8: Compute G
    val T1 = Matrix.extractRows(T, 0, m)
    val T2 = Matrix.extractRows(T, m, 2 * m)
    val GMatrix = (T1.add(1.0, T2)).mult((T1.add(-1.0, T2)).inv())
    // Compute R
    val R = A2.mult(Matrix.eye(m).add(-1.0, A1.add(1.0, A2.mult(GMatrix))).inv())

    // Compute U
    var U = A1.add(1.0, R.mult(A0))
    if (continues) {
        U = Matrix.scaleMult(U.add(-1.0, Matrix.eye(m)), lamb)
    }

    result = HashMap()
    result["G"] = GMatrix
    result["R"] = R
    result["U"] = U
    return result
}

/**
 * Performs convolution of two arrays.
 */
private fun convolution(a: DoubleArray, b: DoubleArray): DoubleArray {
    val result = DoubleArray(a.size + b.size - 1)
    for (i in a.indices) {
        for (j in b.indices) {
            result[i + j] += a[i] * b[j]
        }
    }
    return result
}

/**
 * Computes an orthogonal basis for the column space of a matrix.
 * This is a simplified implementation - in practice, you might want to use QR decomposition.
 */
private fun orthogonalBasis(matrix: Matrix): Matrix {
    // Simplified orthogonalization using Gram-Schmidt process
    val m = matrix.numRows
    val n = matrix.numCols
    val result = matrix.copy()

    for (j in 0..<n) {
        // Normalize column j
        var norm = 0.0
        for (i in 0..<m) {
            norm += result[i, j] * result[i, j]
        }
        norm = FastMath.sqrt(norm)

        if (norm > 1e-12) {
            for (i in 0..<m) {
                result[i, j] /= norm
            }
        }

        // Orthogonalize remaining columns
        for (k in (j + 1)..<n) {
            var dot = 0.0
            for (i in 0..<m) {
                dot += result[i, j] * result[i, k]
            }
            for (i in 0..<m) {
                result[i, k] -= dot * result[i, j]
            }
        }
    }

    return result
}

/**
 * Simplified Schur decomposition.
 * In a full implementation, this would use eigenvalue computation.
 */
private fun schurDecomposition(matrix: Matrix): Pair<Matrix, Matrix> {
    // This is a placeholder - actual Schur decomposition would require eigenvalue computation
    val n = matrix.numRows
    val T = Matrix.eye(n)
    val D = matrix.copy()
    return Pair(T, D)
}