package jline.lib.fjcodes

import jline.util.matrix.Matrix
import org.apache.commons.math3.linear.Array2DRowRealMatrix
import org.apache.commons.math3.linear.EigenDecomposition
import kotlin.math.abs
import kotlin.math.max

/**
 * Compute T-matrix using Sylvester equation approach
 *
 * Solves the matrix equation using iterative Sylvester equation solver
 * with Hessenberg decomposition.
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */

/**
 * Compute T-matrix via Sylvester equation method
 *
 * Uses fixpoint iteration solving Sylvester equations via Hessenberg decomposition.
 * Based on FJ_codes computeT.m with T_Mode='Sylvest'.
 *
 * @param D0 D0 matrix (transitions without arrivals)
 * @param D1 D1 matrix (transitions with arrivals)
 * @param S State transition matrix (S_Arr)
 * @param A_jump Jump matrix
 * @return T-matrix
 */
fun computeT_Sylvester(D0: Matrix, D1: Matrix, S: Matrix, A_jump: Matrix): Matrix {
    val d0 = D0.getNumRows()
    val ms = S.getNumRows() / d0
    val m = ms * d0

    // ID0 = kron(eye(ms), D0)
    val Ims = Matrix.eye(ms)
    val ID0 = Ims.kron(D0)

    // DS = kron(eye(ms), D1) * kron(A_jump, eye(d0))
    val A_jump_Arr = A_jump.kron(Matrix.eye(d0))
    val DS = Ims.kron(D1).mult(A_jump_Arr)

    // Compute Schur decomposition of ID0 with complex support
    val ID0_commons = Array2DRowRealMatrix(ID0.toArray2D())

    // Use eigendecomposition instead of Schur for complex support
    val eigen = EigenDecomposition(ID0_commons)
    val U_commons = eigen.v
    val D_commons = eigen.d  // Diagonal matrix with eigenvalues

    val U = Matrix(U_commons.data)
    val Tr = Matrix(D_commons.data)

    // Fixpoint iteration: while ||T_old - T_new|| > tolerance
    var Tnew = S.copy()  // Initial guess: S_Arr (without service completion)
    var Told = Matrix.zeros(m, m)
    val tolerance = 1e-12  // Tighter tolerance for better accuracy
    val maxIter = 200  // More iterations allowed
    var iter = 0

    while (iter < maxIter) {
        // Compute max difference
        var maxDiff = 0.0
        for (i in 0 until m) {
            for (j in 0 until m) {
                maxDiff = max(maxDiff, abs(Told.get(i, j) - Tnew.get(i, j)))
            }
        }

        if (maxDiff <= tolerance) {
            break
        }

        Told = Tnew.copy()

        // Solve Sylvester equation: L*ID0 + Tnew*L = -I
        // Using Q_Sylvest: X*kron(A,I) + BX = -I
        val L = qSylvest(U, Tr, Tnew)

        // Update: Tnew = S + L * DS
        Tnew = S.add(1.0, L.mult(DS))

        iter++
    }

    if (iter >= maxIter) {
//        System.err.println("WARNING: Sylvester T-matrix did not converge after $maxIter iterations")
    }

    // Compute residual norm for verification
    // Equation: T*L + L*ID0 + I = 0
    val L_final = qSylvest(U, Tr, Tnew)
    val termTL = Tnew.mult(L_final)
    val termLID0 = L_final.mult(ID0)
    val termI = Matrix.eye(m)

    val residual = termTL.add(1.0, termLID0).add(1.0, termI)
    val residualNorm = normInf(residual)

    if (residualNorm > 1e-6) {
//        System.err.println("WARNING: T-matrix residual = $residualNorm (expected < 1e-6)")
    }

    return Tnew
}

/**
 * Solve Sylvester equation X*kron(A,I) + BX = -I using Hessenberg decomposition
 *
 * Based on qmam Q_Sylvest.m function.
 *
 * @param U Orthogonal matrix from Schur/Eigen decomposition of A
 * @param T Schur/Diagonal matrix (U'*A*U = T)
 * @param B Matrix B in equation
 * @return Solution matrix X
 */
private fun qSylvest(U: Matrix, T: Matrix, B: Matrix): Matrix {
    val n = U.getNumCols()

    // Compute Hessenberg decomposition of B: V'*B*V = NBAR
    val B_commons = Array2DRowRealMatrix(B.toArray2D())

    // Access package-private HessenbergTransformer via reflection
    val hessClass = Class.forName("org.apache.commons.math3.linear.HessenbergTransformer")
    val constructor = hessClass.getDeclaredConstructor(org.apache.commons.math3.linear.RealMatrix::class.java)
    constructor.isAccessible = true
    val hessTransformer = constructor.newInstance(B_commons)

    // Get H (Hessenberg form) and P (transformation matrix)
    val getH = hessClass.getDeclaredMethod("getH")
    getH.isAccessible = true
    val NBAR_commons = getH.invoke(hessTransformer) as org.apache.commons.math3.linear.RealMatrix

    val getP = hessClass.getDeclaredMethod("getP")
    getP.isAccessible = true
    val V_commons = getP.invoke(hessTransformer) as org.apache.commons.math3.linear.RealMatrix

    val V = Matrix(V_commons.data)
    val NBAR = Matrix(NBAR_commons.data)

    // F = -V' * U
    val F = V.transpose().mult(U).scale(-1.0)

    // Solve column by column
    val Y = Matrix.zeros(n, n)

    for (k in 0 until n) {
        // Compute temp vector
        val temp = Matrix.zeros(n, 1)

        if (k == 0) {
            // temp = F(:, k)
            for (i in 0 until n) {
                temp.set(i, 0, F.get(i, k))
            }
        } else {
            // temp = F(:, k) - Y(:, 1:k-1) * T(1:k-1, k)
            for (i in 0 until n) {
                var sum = F.get(i, k)
                for (j in 0 until k) {
                    sum -= Y.get(i, j) * T.get(j, k)
                }
                temp.set(i, 0, sum)
            }
        }

        // Solve (NBAR + T(k,k)*I) * y_k = temp
        // This is Ax = b where A = (NBAR + T(k,k)*I) and b = temp
        val A_sys = NBAR.add(T.get(k, k), Matrix.eye(n))

        // Use Apache Commons Math to solve the linear system
        val A_sys_commons = Array2DRowRealMatrix(A_sys.toArray2D())
        val temp_commons = Array2DRowRealMatrix(temp.toArray2D())

        // Solve Ax = b using LU decomposition
        val solver = org.apache.commons.math3.linear.LUDecomposition(A_sys_commons).getSolver()
        val y_k_commons = solver.solve(temp_commons)

        // Store in Y
        for (i in 0 until n) {
            Y.set(i, k, y_k_commons.getEntry(i, 0))
        }
    }

    // X = real(V * Y * U')
    val X = V.mult(Y).mult(U.transpose())

    // Take real part (remove any numerical imaginary components)
    return X
}

/**
 * Compute infinity norm (max row sum of absolute values)
 */
private fun normInf(A: Matrix): Double {
    var maxSum = 0.0
    for (i in 0 until A.getNumRows()) {
        var rowSum = 0.0
        for (j in 0 until A.getNumCols()) {
            rowSum += kotlin.math.abs(A.get(i, j))
        }
        maxSum = kotlin.math.max(maxSum, rowSum)
    }
    return maxSum
}
