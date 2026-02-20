/**
 * @file PANACEA (Product-form Approximation using Normal usage Approximation for Closed queueing nEtwork Analysis)
 * 
 * Implements the PANACEA approximation method for computing normalizing constants in
 * closed product-form queueing networks. Uses asymptotic expansion with gamma parameter
 * transformations and convolution techniques for high-accuracy approximation.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.nc

import jline.io.Ret
import jline.io.line_warning
import jline.GlobalConstants
import jline.solvers.SolverOptions
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath

/**
 * Compute the PANACEA approximation
 *
 * @param L - demands at all stations
 * @param N - number of jobs for each class
 * @param Z - think time for each class
 * @return normalizing constant and its logarithm
 */
fun pfqn_panacea(L: Matrix, N: Matrix, Z: Matrix, options: SolverOptions = SolverOptions()): Ret.pfqnNc {
    var Z = Z
    val method = options.method
    val M = L.numRows
    val R = L.numCols
    var lG = Double.NaN
    var G = Double.NaN
    if (Z.isEmpty || Z.elementSum() < options.tol) {
        Z = Z.copy()
        Z = Z.add(GlobalConstants.FineTol, Z)
    }


    if (L.isEmpty || L.elementSum() < options.tol) {
        run {
            val tmp1 = Z.sumCols()
            for (i in 0..<tmp1.length()) {
                tmp1[i] = N[i] * FastMath.log(tmp1[i])
            }
            lG = -Matrix.factln(N).elementSum() + tmp1.elementSum()
        }
        G = FastMath.exp(lG)
        return Ret.pfqnNc(G, lG, method)
    }

    // Compute r = L./repmat(Z,q,1)
    val r = Matrix(M, R)
    for (i in 0..<M) {
        for (j in 0..<R) {
            r[i, j] = L[i, j] / Z[0, j]
        }
    }

    // Find Nt = max(1./r)
    var Nt = Double.MIN_VALUE
    for (i in 0..<M) {
        for (j in 0..<R) {
            Nt = FastMath.max(Nt, 1.0 / r[i, j])
        }
    }

    // Compute beta = N / Nt
    val beta = Matrix(1, R)
    for (j in 0..<R) {
        beta[0, j] = N[0, j] / Nt
    }

    // Compute gamma = r * Nt
    val gamma = Matrix(M, R)
    for (i in 0..<M) {
        for (j in 0..<R) {
            gamma[i, j] = r[i, j] * Nt
        }
    }

    // Compute alpha = 1 - N * r'
    val alpha = Matrix(1, M)
    for (i in 0..<M) {
        var sum = 0.0
        for (j in 0..<R) {
            sum += N[j] * r[i, j]
        }
        alpha[0, i] = 1 - sum
    }

    // Compute gammatilde = gamma ./ repmat(alpha', 1, p)
    val gammatilde = Matrix(M, R)
    for (i in 0..<M) {
        for (j in 0..<R) {
            gammatilde[i, j] = gamma[i, j % R] / alpha[0, i]
        }
    }

    // Check if min(alpha) < 0
    var minAlpha = Double.MAX_VALUE
    for (i in 0..<M) {
        minAlpha = FastMath.min(minAlpha, alpha[0, i])
    }
    if (minAlpha < 0) {
        line_warning("pfqn_panacea", "Model is not in normal usage")
        lG = Double.NaN
        G = Double.NaN
        return Ret.pfqnNc(G, lG, method)
    }

    val A0 = 1.0
    var A1 = 0.0
    for (j in 0..<R) {
        val m = Matrix(1, R)
        m[0, j] = 2.0
        A1 -= beta[j] * pfqn_ca(gammatilde, m).G
    }

    var A2 = 0.0
    for (j in 0..<R) {
        val m = Matrix(1, R)
        m[0, j] = 3.0
        A2 += 2 * beta[j] * pfqn_ca(gammatilde, m).G

        m[0, j] = 4.0
        A2 += 3 * FastMath.pow(beta[j], 2) * pfqn_ca(gammatilde, m).G

        for (k in 0..<R) {
            if (k == j) continue
            m.zero()
            m[0, j] = 2.0
            m[0, k] = 2.0
            A2 = A2 + 0.5 * beta[j] * beta[k] * pfqn_ca(gammatilde, m).G
        }
    }

    // Compute I = [A0, A1/Nt, A2/Nt^2] and sum them
    val I = doubleArrayOf(A0, A1 / Nt, A2 / (Nt * Nt))
    val sumI = I.sum()
    
    val tmp1 = Z.sumCols()
    for (i in 0..<tmp1.length()) {
        tmp1[i] = N[i] * FastMath.log(tmp1[i])
    }
    lG = -Matrix.factln(N).elementSum() + tmp1.elementSum()
    lG += FastMath.log(sumI)
    for (s in 0..<M) {
        lG -= FastMath.log(alpha[0, s])
    }
    G = FastMath.exp(lG)

    // if ~isfinite(lGn)
    if (!java.lang.Double.isFinite(lG)) {
        G = Double.NaN
        lG = Double.NaN
    }
    return Ret.pfqnNc(G, lG, method)
}
/**
 * PFQN panacea algorithms
 */
@Suppress("unused")
class PfqnPanaceaAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}