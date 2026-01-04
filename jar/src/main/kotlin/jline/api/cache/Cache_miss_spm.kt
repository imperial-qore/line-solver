/**
 * @file Cache Miss Analysis via SPM
 * 
 * Estimates cache miss rates using SPM method for partial differential equations.
 * Provides PDE-based approximations for cache performance metrics in large-scale
 * systems where exact analysis becomes computationally intractable.
 * 
 * @since LINE 3.0
 */
package jline.api.cache

import jline.io.Ret
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import org.apache.commons.math3.util.FastMath

/**
 * Estimates the cache miss rate and related metrics using the SPM method for PDEs.
 * This method computes the miss rate (M), user-specific miss rates (MU), item-specific miss rates (MI),
 * initial state probabilities (pi0), and the logarithm of the normalizing constant (lZ).
 *
 * @param gamma  - Matrix representing the cache access factors.
 * @param m      - Matrix representing the cache capacity vector.
 * @param lambda - MatrixCell representing the request rates for different users or items.
 * @return cacheMissSpmReturn - An object containing the miss rate metrics and probabilities.
 *//* Example:
    Matrix gamma = new Matrix("[0.4,0.2; 0.4,0.2; 0.4,0.1; 0.4,0.2; 0.4,0.2]");
    Matrix m = new Matrix("[2,1]");
    MatrixCell lambda = new MatrixCell(2);
    lambda.set(0,new Matrix("[1;1;1;1;1]"));
    lambda.set(1,new Matrix("[0;0;1;1;1]"));
    cacheMissSpmReturn ret = cache_miss_spm(gamma, m, lambda);
    new Matrix(ret.MU).print();
 */

fun cache_miss_spm(gamma: Matrix, m: Matrix, lambda: MatrixCell): Ret.cacheMissSpm {
    val ma = m.copy()
    ma[0, 0] = ma.value() + 1

    val lE = cache_spm(gamma, m).lZ
    val lEa = cache_spm(gamma, ma).lZ

    val M = FastMath.exp(lEa - lE)

    val u = lambda.size()
    val n = lambda[0].numRows
    val pi0 = DoubleArray(n)
    val MU = DoubleArray(u)
    val lE1 = DoubleArray(n)

    for (k in 0..<n) {
        if (gamma.getRow(k).elementSum() > 0) {
            val subGamma = gamma.copy()
            subGamma.removeRows(setOf(k))

            lE1[k] = cache_spm(subGamma, m).lZ
            pi0[k] = FastMath.exp(lE1[k] - lE)

            //                if (pi0[k] > 1 || pi0[k] < 0) {
//                    result = cache_spm(subGamma, m);
//                    lE1[k] = result.lZ;
//                    pi0[k] = FastMath.exp(lE1[k] - lE);
//                }
            for (v in 0..<u) {
                MU[v] += lambda[v][k, 0] * pi0[k]
            }
        }
    }

    val MI = DoubleArray(n)
    for (k in 0..<n) {
        if (gamma.getRow(k).elementSum() > 0) {
            MI[k] = lambda.cellsum(k, 0) * FastMath.exp(lE1[k] - lE)
        } else {
            MI[k] = 0.0
        }
    }

    return Ret.cacheMissSpm(M, MU, MI, pi0, lE)
}
/**
 * Cache miss rayint algorithms
 */
@Suppress("unused")
class CacheMissRayintAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}