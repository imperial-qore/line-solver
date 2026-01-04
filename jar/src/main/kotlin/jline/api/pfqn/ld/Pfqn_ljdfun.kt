/**
 * @file Limited joint-dependent function evaluation
 *
 * Evaluates limited joint-dependent (LJD) scaling functions using a lookup table
 * indexed by per-class population vector (n1, n2, ..., nR).
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.ld

import jline.util.Maths
import jline.util.Utils
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath
import java.lang.Double

/**
 * Evaluate limited joint-dependent (LJD) function
 *
 * @param nvec       Per-class population vector [n1, n2, ..., nK]
 * @param ljdscaling Linearized scaling tables indexed by station index, or null if none
 * @param ljdcutoffs Per-class cutoffs [N1, N2, ..., NK] indexed by station index, or null if none
 * @param M          Number of stations
 * @param nservers   Number of servers per station (optional)
 * @return Scaling factor vector (M x 1)
 */
fun pfqn_ljdfun(
    nvec: Matrix,
    ljdscaling: List<Matrix?>?,
    ljdcutoffs: List<Matrix?>?,
    M: Int,
    nservers: Matrix?
): Matrix {
    val K = nvec.length()
    val r = Matrix(M, 1)
    r.fill(1.0)
    val alpha = 20.0

    if (ljdscaling == null || ljdscaling.isEmpty()) {
        return r
    }

    for (i in 0 until M) {
        val scaling = ljdscaling.getOrNull(i)
        val cutoffs = ljdcutoffs?.getOrNull(i)

        if (scaling != null && !scaling.isEmpty && cutoffs != null) {
            // Clamp population to cutoffs
            val nClamped = Matrix(1, K)
            for (k in 0 until K) {
                nClamped.set(0, k, FastMath.min(nvec.get(k), cutoffs.get(k)))
            }

            // Compute linearized index (0-indexed for Java)
            val idx = ljd_linearize(nClamped, cutoffs)

            // Look up scaling value
            if (idx < scaling.length()) {
                r.set(i, 0, scaling.get(idx))
            }
        }

        // Handle servers
        if (nservers != null && !nservers.isEmpty) {
            if (!Utils.isInf(nservers.get(i))) {
                val totalN = nvec.elementSum()
                val softminValue = r.get(i, 0) / Maths.softmin(totalN, nservers.get(i), alpha)
                if (Double.isNaN(softminValue)) {
                    r.set(i, 0, 1.0 / FastMath.min(totalN, nservers.get(i)))
                } else {
                    r.set(i, 0, softminValue)
                }
            }
        }
    }
    return r
}

/**
 * Convert per-class population vector to linearized index (0-indexed)
 *
 * @param nvec Per-class population vector [n1, n2, ..., nK]
 * @param cutoffs Per-class cutoffs [N1, N2, ..., NK]
 * @return 0-indexed linearized index
 */
fun ljd_linearize(nvec: Matrix, cutoffs: Matrix): Int {
    val K = nvec.length()
    var idx = 0
    var multiplier = 1

    for (k in 0 until K) {
        val nk = FastMath.min(nvec.get(k), cutoffs.get(k)).toInt()
        idx += nk * multiplier
        multiplier *= (cutoffs.get(k).toInt() + 1)
    }
    return idx
}

/**
 * Convert linearized index back to per-class population vector
 *
 * @param idx 0-indexed linearized index
 * @param cutoffs Per-class cutoffs [N1, N2, ..., NK]
 * @return Per-class population vector [n1, n2, ..., nK]
 */
fun ljd_delinearize(idx: Int, cutoffs: Matrix): Matrix {
    val K = cutoffs.length()
    val nvec = Matrix(1, K)
    var remaining = idx

    for (k in 0 until K) {
        val base = cutoffs.get(k).toInt() + 1
        nvec.set(0, k, (remaining % base).toDouble())
        remaining /= base
    }
    return nvec
}

/**
 * PFQN ljdfun algorithms
 */
@Suppress("unused")
class PfqnLjdfunAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}
