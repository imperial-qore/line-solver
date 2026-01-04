/**
 * @file COMOM normalizing constant method for load-dependent repairman models
 * 
 * Implements the Convolution Method of Moments (COMOM) for computing normalizing constants in
 * load-dependent repairman queueing models. Handles state-dependent service rates and provides
 * exact normalization for closed networks with load-dependent service mechanisms.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.ld

import jline.io.Ret
import jline.GlobalConstants
import jline.api.pfqn.nc.pfqn_ca
import jline.api.pfqn.nc.pfqn_nc_sanitize
import jline.solvers.SolverOptions
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath

/**
 * Run the COMOM normalizing constant solution method on a repairman model
 *
 * @param L       - demands at all stations
 * @param N       - number of jobs for each class
 * @param Z       - think time for each class
 * @param mu      - load-depedent scalings
 * @param options - solver options
 * @return normalizing constant and its logarithm
 */

fun pfqn_comomrm_ld(L: Matrix, N: Matrix, Z: Matrix, mu: Matrix, options: SolverOptions): Ret.pfqnComomrmLd {
    var L = L.copy()
    var N = N.copy()
    var Z = Z.copy()
    val mu = mu.copy()
    val atol = options.tol
    N.ceilEq()
    var M = L.numRows
    var R = L.numCols
    val Nt = N.elementSum().toInt()
    Z = Z.sumCols()

    if (Z.elementSum() < GlobalConstants.Zero) {
        val zset = BooleanArray(M)
        val ignoreIndex = HashSet<Int>()
        val colset = BooleanArray(R)  // Fix: size should be R (number of classes), not M
        val OneToNt = Matrix.createLike(mu.getRow(0))
        for (nt in 1..OneToNt.length()) {
            OneToNt[nt - 1] = nt.toDouble()
        }

        // Set all columns (classes) to be selected
        for (j in 0..<R) {
            colset[j] = true
        }

        for (i in 0..<M) {
            if (mu.getRow(i).sub(OneToNt).norm() < atol) {
                zset[i] = true
                ignoreIndex.add(i)
            } else {
                zset[i] = false
            }
        }
        Z = L.getSlice(zset, colset)
        mu.removeRows(ignoreIndex)
        L.removeRows(ignoreIndex)
    }

    if (L.elementSum() < GlobalConstants.Zero) {
        val ca = pfqn_ca(L, N, Z)
        val G = ca.G
        val lG = ca.lG
        val prob = Matrix(Nt + 1, 1)
        prob[Nt] = 1.0
        return Ret.pfqnComomrmLd(G, lG, prob)
    }

    val lG0: Double

    val sanitizedResult = pfqn_nc_sanitize(Matrix(1, R), L, N, Z, atol)
    L = sanitizedResult.L
    N = sanitizedResult.N
    Z = sanitizedResult.Z
    lG0 = sanitizedResult.lGremaind
    M = L.numRows
    R = L.numCols

    if (Z.isEmpty) {
        if (L.isEmpty) {
            val G = FastMath.exp(lG0)
            val lG = lG0
            val prob = Matrix(Nt + 1, 1)
            prob[0] = 1.0
            return Ret.pfqnComomrmLd(G, lG, prob)
        }
        Z = Matrix(1, R)
    } else if (L.isEmpty) {
        L = Matrix(1, R)
    }

    if (M == 0) {
        val G = FastMath.exp(lG0)
        val lG = lG0
        val prob = Matrix(Nt + 1, 1)
        prob[Nt] = 1.0
        return Ret.pfqnComomrmLd(G, lG, prob)
    }

    require(M == 1) { "The solver accepts at most a single queueing station." }

    var h = Matrix(Nt + 1, 1)
    h[Nt] = 1.0
    val scale = Matrix(Nt, 1)
    var nt = 0

    for (r in 0..<R) {
        val Tr = Matrix.eye(Nt + 1).scale(Z[r])
        for (i in 0..<Tr.numCols - 1) {
            Tr[i, i + 1] = L[r] * (Nt - i) / mu[Nt - i - 1]
        }
        var hT: Matrix
        var nr = 0
        while (nr < N[r]) {
            hT = Tr.copy().scale(1.0 / (1.0 + nr))
            h = hT.mult(h)
            scale[nt] = FastMath.abs(h.sort().elementSum())
            h.absEq()
            h.scaleEq(1.0 / scale[nt])
            nt++
            nr++
        }
    }

    val lG = lG0 + scale.log().elementSum()
    val G = FastMath.exp(lG)
    val prob = h.reverse().scale(1 / G)
    prob.divideEq(prob.elementSum())

    return Ret.pfqnComomrmLd(G, lG, prob)
}
/**
 * PFQN comomrm ld algorithms
 */
@Suppress("unused")
class PfqnComomrmLdAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}