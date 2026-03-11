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
 * Multi-linear interpolation for LJCD scaling tables.
 *
 * Performs multi-linear interpolation of a throughput value from an LJCD
 * scaling table for non-integer population vectors. For K classes, interpolates
 * between 2^K corner points of the hypercube containing the population vector.
 *
 * The interpolation uses floor/ceil corners in each dimension with fractional
 * weights, producing smooth transitions between discrete table entries.
 *
 * @param nvec Continuous population vector [n1, n2, ..., nK] (clamped to cutoffs)
 * @param cutoffs Per-class cutoffs [N1, N2, ..., NK]
 * @param table Linearized throughput table (1-D matrix indexed by ljd_linearize)
 * @param K Number of classes
 * @return Interpolated throughput value
 */
fun ljcd_interpolate(nvec: Matrix, cutoffs: Matrix, table: Matrix, K: Int): kotlin.Double {
    // Get floor and ceiling for each dimension
    val nFloor = IntArray(K)
    val nCeil = IntArray(K)
    val frac = DoubleArray(K)

    for (i in 0 until K) {
        val v = nvec.get(i)
        val c = cutoffs.get(i).toInt()
        nFloor[i] = FastMath.max(0, FastMath.min(Math.floor(v).toInt(), c))
        nCeil[i] = FastMath.max(0, FastMath.min(Math.ceil(v).toInt(), c))
        frac[i] = v - Math.floor(v)
    }

    // Handle edge case: all integer values
    var allInteger = true
    for (i in 0 until K) {
        if (FastMath.abs(frac[i]) >= 1e-10) {
            allInteger = false
            break
        }
    }

    if (allInteger) {
        val cornerVec = Matrix(1, K)
        for (i in 0 until K) {
            cornerVec.set(0, i, nFloor[i].toDouble())
        }
        val idx = ljd_linearize(cornerVec, cutoffs)
        return if (idx < table.length()) {
            table.get(idx)
        } else {
            0.0
        }
    }

    // Multi-linear interpolation over 2^K corners
    var result = 0.0
    val numCorners = 1 shl K // 2^K

    for (corner in 0 until numCorners) {
        // Build corner point: bit i determines floor (0) or ceil (1) for class i
        val cornerVec = Matrix(1, K)
        var weight = 1.0

        for (i in 0 until K) {
            if ((corner shr i) and 1 == 1) {
                // Use ceiling for this dimension
                cornerVec.set(0, i, nCeil[i].toDouble())
                weight *= frac[i]
            } else {
                // Use floor for this dimension
                cornerVec.set(0, i, nFloor[i].toDouble())
                weight *= (1.0 - frac[i])
            }
        }

        // Look up table value at corner point
        val idx = ljd_linearize(cornerVec, cutoffs)
        if (idx < table.length()) {
            result += weight * table.get(idx)
        }
    }

    return result
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
