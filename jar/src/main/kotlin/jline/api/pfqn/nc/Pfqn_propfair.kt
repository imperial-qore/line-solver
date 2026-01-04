/**
 * @file Proportionally fair allocation approximation for normalizing constants
 * 
 * Implements the proportionally fair allocation method using convex optimization
 * for computing normalizing constants in closed queueing networks. Based on Schweitzer's
 * approach and Walton's proportional fairness theory for multi-class networks.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.nc

import de.xypron.jcobyla.Calcfc
import de.xypron.jcobyla.Cobyla
import jline.io.Ret
import jline.GlobalConstants
import jline.util.Maths
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath

/**
 * Compute the proportionally fair allocation approximation
 *
 * Estimate the normalizing constant using a convex optimization program
 * that is asymptotically exact in models with single-server PS queues only.
 * The underlying optimization program is convex.
 * The script implements a heuristic to estimate the solution in the
 * presence of delay stations.
 *
 * References:
 * Schweitzer, P. J. (1979). Approximate analysis of multiclass closed networks of
 * queues. In Proceedings of the International Conference on Stochastic Control
 * and Optimization. Free Univ., Amsterdam.
 *
 * Walton, Proportional fairness and its relationship with multi-class
 * queueing networks, 2009.
 *
 * @param L - demands at all stations
 * @param N - number of jobs for each class
 * @param Z - think time for each class
 * @return normalizing constant, its logarithm, and performance metrics
 */
fun pfqn_propfair(L: Matrix, N: Matrix, Z: Matrix): Ret.pfqnNcXQ {
    val M = L.numRows
    val R = L.numCols

    val Nvec = N.toArray1D()
    val Zvec = Z.toArray1D()
    val Lmat = L.toArray2D()

    /* ---------- objective & constraints for COBYLA ---------------- */
    val objFun = Calcfc { n, m, x, con ->
        var obj = 0.0
        for (r in 0..<R) {
            obj += ((Nvec[r] - x[r] * Zvec[r]) * FastMath.log(FastMath.abs(x[r]) + GlobalConstants.FineTol))
        }
        for (i in 0..<M) {
            con[i] = 1.0
            for (r in 0..<R) {
                con[i] -= Lmat[i][r] * x[r]
            }
        }
        System.arraycopy(x, 0, con, M, R) // x ≥ 0
        -obj // maximise obj
    }

    val Xasy = DoubleArray(R)
    Cobyla.findMinimum(objFun, R, M + R, Xasy, 1.0, 1.0e-8, 2, 10000)

    /* -------------- updated lG -------------- */
    var lG = 0.0
    for (r in 0..<R) {
        val x = Xasy[r]
        lG += ((N[0, r] - x * Z[0, r]) * FastMath.log(1.0 / (x + GlobalConstants.FineTol)))
    }
    for (r in 0..<R) {
        lG -= FastMath.log(Maths.fact(Xasy[r] * Z[0, r]))
    }

    val G = FastMath.exp(lG)

    val Xa = Matrix(Xasy)
    val Qa = Xa.copy()
    Qa.fill(Double.NaN)

    return Ret.pfqnNcXQ(G, lG, Xa, Qa, "propfair")
}
/**
 * PFQN propfair algorithms
 */
@Suppress("unused")
class PfqnPropfairAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}