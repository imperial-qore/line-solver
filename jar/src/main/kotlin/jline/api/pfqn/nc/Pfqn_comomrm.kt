/**
 * @file Convolution Method of Moments for Repairman models (COMOM-RM)
 * 
 * Implements the Convolution Method of Moments specialized for repairman queueing models
 * with a single queueing station and delay server. Provides exact normalizing constant
 * computation for multi-class closed networks with think times.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.nc

import jline.io.Ret
import jline.GlobalConstants
import jline.util.Maths
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath

/**
 * Compute the normalizing constant of a repairmen model using COMOM
 *
 * @param L    - demands at all stations
 * @param N    - number of jobs for each class
 * @param Z    - think time for each class
 * @param m    - multiplicy of queueing station
 * @param atol - absolute numerical tolerance
 * @return sanitized parameters
 */

fun pfqn_comomrm(L: Matrix, N: Matrix, Z: Matrix, m: Int?, atol: Double): Ret.pfqnComomrm {
    var m = m
    val M = L.numRows
    val R = L.numCols
    if (M != 1) {
        throw RuntimeException("pfqn_comomrm: The solver accepts at most a single queueing station.")
    }
    if (m == null) {
        m = 1
    }
    val lambda = N.copy()
    lambda.fill(0.0)
    val ret = pfqn_nc_sanitize(lambda, L, N, Z, atol)
    val L_new = ret.L
    val N_new = ret.N
    val Z_new = ret.Z
    val lG0 = ret.lGremaind
    val zerothinktimes: MutableList<Int> = ArrayList()

    for (i in 0..<Z_new.length()) {
        if (Z_new[i] < GlobalConstants.FineTol) {
            zerothinktimes.add(i)
        }
    }
    val nvec = Matrix(1, R)
    nvec.fill(0.0)

    var lh: Matrix

    if (!zerothinktimes.isEmpty()) {
        for (i in zerothinktimes.indices) {
            nvec[zerothinktimes[i]] = N_new[zerothinktimes[i]]
        }
        lh = Matrix(0, 1)
        val tmp = Matrix(1, 1)
        tmp[0] = Maths.factln(nvec.elementSum() + m) - Matrix.factln(nvec).elementSum()
        lh = Matrix.concatRows(lh, tmp, null)

        for (i in zerothinktimes.indices) {
            val nvec_s = nvec.copy()
            nvec_s[i] = nvec_s[i] - 1
            tmp[0] = Maths.factln(nvec_s.elementSum() + m) - Matrix.factln(nvec_s).elementSum()
            lh = Matrix.concatRows(lh, tmp, null)
        }
        tmp[0] = Maths.factln(nvec.elementSum() + m - 1) - Matrix.factln(nvec).elementSum()
        lh = Matrix.concatRows(lh, tmp, null)
        for (i in zerothinktimes.indices) {
            val nvec_s = nvec.copy()
            nvec_s[i] = nvec_s[i] - 1
            tmp[0] = Maths.factln(nvec_s.elementSum() + m - 1) - Matrix.factln(nvec_s).elementSum()
            lh = Matrix.concatRows(lh, tmp, null)
        }
    } else {
        lh = Matrix(2, 1)
        lh.fill(0.0)
    }
    var h = lh.copy()
    for (i in 0..<h.length()) {
        h[i] = FastMath.exp(h[i])
    }

    val lG: Double
    val lGbasis: Matrix

    if (zerothinktimes.size == R) {
        lGbasis = h.copy()
        for (i in 0..<lGbasis.length()) {
            lGbasis[i] = FastMath.log(lGbasis[i])
        }
        lG = lG0 + FastMath.log(h[h.length() - 1 - R])
    } else {
        val scale = Matrix(1, N_new.elementSum().toInt())
        scale.fill(1.0)
        var nt = nvec.elementSum()
        var h_1 = h.copy()
        for (r in zerothinktimes.size + 1..R) {
            var F1r: Matrix? = null
            var F2r: Matrix? = null
            var Nr = 1
            while (Nr <= N_new[r - 1]) {
                nvec[r - 1] = nvec[r - 1] + 1
                if (Nr == 1) {
                    if (r > zerothinktimes.size + 1) {
                        val hr = Matrix(2 * r, 1)
                        hr.fill(0.0)
                        for (i in 0..<r - 1) {
                            hr[i] = h[i]
                        }
                        for (i in r..<2 * r - 1) {
                            hr[i] = h[i - 1]
                        }
                        h = hr
                        if (nt > 0) {
                            h[r - 1] = h_1[0] / scale[nt.toInt() - 1]
                            h[h.length() - 1] = h_1[r - 1] / scale[nt.toInt() - 1]
                        }
                    }

                    val A12 = Matrix(r, r)
                    A12.fill(0.0)
                    A12[0, 0] = -1
                    for (s in 1..<r) {
                        A12[s, 0] = N_new[s - 1]
                        A12[s, s] = -Z_new[s - 1]
                    }

                    var B2r = Matrix.eye(r)
                    val B2r_tmp = Matrix.eye(r)
                    for (i in 0..<r) {
                        B2r[i, i] = m * L_new[0, r - 1]
                        B2r_tmp[i, i] = Z_new[r - 1]
                    }
                    B2r = Matrix.concatColumns(B2r, B2r_tmp, null)

                    val iC = Matrix.eye(r)
                    for (i in 0..<r) {
                        iC[i, i] = 1.0 / m
                        iC[0, i] = 1.0 / m
                    }
                    iC[0, 0] = -1

                    F1r = Matrix(2 * r, 2 * r)
                    F1r.fill(0.0)
                    F1r[0, 0] = 1

                    F2r = iC.mult(A12).mult(B2r)
                    F2r = Matrix.concatRows(F2r, B2r, null)
                }
                h_1 = h
                val tmp_mat = F1r!!.copy()
                for (i in 0..<tmp_mat.numRows) {
                    for (j in 0..<tmp_mat.numCols) {
                        tmp_mat[i, j] = F1r[i, j] + F2r!![i, j] / nvec[r - 1]
                    }
                }
                h = tmp_mat.mult(h_1)
                nt = nvec.elementSum()
                scale[nt.toInt() - 1] = FastMath.abs(h.elementSum())
                for (i in 0..<h.length()) {
                    h[i] = FastMath.abs(h[i]) / scale[nt.toInt() - 1]
                }
                Nr++
            }
        }
        val log_scale = scale.copy()
        for (i in 0..<log_scale.length()) {
            log_scale[i] = FastMath.log(log_scale[i])
        }
        lG = lG0 + FastMath.log(h[h.length() - 1 - (R - 1)]) + log_scale.elementSum()
        lGbasis = h.copy()
        for (i in 0..<lGbasis.length()) {
            lGbasis[i] = FastMath.log(h[i]) + log_scale.elementSum()
        }
    }
    return Ret.pfqnComomrm(lG, lGbasis)
}
/**
 * PFQN comomrm algorithms
 */
@Suppress("unused")
class PfqnComomrmAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}