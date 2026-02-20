/**
 * @file Main load-dependent normalizing constant computation with method selection
 * 
 * Provides the main entry point for computing normalizing constants in load-dependent
 * queueing networks with automatic method selection and preprocessing. Handles scaling,
 * class separation, and method dispatch for various load-dependent solution algorithms.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.ld

import jline.io.Ret
import jline.io.line_warning
import jline.io.mfilename
import jline.GlobalConstants
import jline.api.pfqn.nc.pfqn_nrl
import jline.api.pfqn.nc.pfqn_nrp
import jline.api.pfqn.nc.pfqn_rd
import jline.solvers.SolverOptions
import jline.util.Maths
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath

/**
 * Main method to compute the normalizing constant of a load-dependent model
 *
 * @param L       - demands at all stations
 * @param N       - number of jobs for each class
 * @param Z       - think time for each class
 * @param mu      - load-dependent scalings
 * @param options - solver options
 * @return normalizing constant and its logarithm
 */
fun pfqn_ncld(L: Matrix, N: Matrix, Z: Matrix, mu: Matrix, options: SolverOptions): Ret.pfqnNc {
    var lG = Double.NaN
    var G = Double.NaN
    var method = options.method

    var mu_new: Matrix
    if (N.elementSum().toInt() >= mu.numCols) {
        mu_new = mu.copy()
    } else {
        mu_new = Matrix(mu.numRows, 0)
        var i = 0
        while (i < N.elementSum()) {
            val mu_col_i = Matrix(mu.numRows, 1)
            Matrix.extract(mu, 0, mu.numRows, i, i + 1, mu_col_i, 0, 0)
            mu_new = Matrix.concatColumns(mu_new, mu_col_i, null)
            i++
        }
    }

    var L_new = Matrix(L.numRows, 0)
    var N_new = Matrix(N.numRows, 0)
    var Z_new = Matrix(Z.numRows, 0)
    for (i in 0..<N.length()) {
        if (FastMath.abs(N[i]) >= GlobalConstants.FineTol) {
            val L_col_i = Matrix.extractColumn(L, i, null)
            val N_col_i = Matrix.extractColumn(N, i, null)
            val Z_col_i = Matrix.extractColumn(Z, i, null)
            L_new = Matrix.concatColumns(L_new, L_col_i, null)
            N_new = Matrix.concatColumns(N_new, N_col_i, null)
            Z_new = Matrix.concatColumns(Z_new, Z_col_i, null)
        }
    }

    var R = N_new.numCols
    val scalevec = Matrix(1, R)
    scalevec.fill(1.0)
    for (r in 0..<R) {
        val L_col_r = Matrix.extractColumn(L_new, r, null)
        val Z_col_r = Matrix.extractColumn(Z_new, r, null)
        scalevec[r] = FastMath.max(L_col_r.elementMax(), Z_col_r.elementMax())
    }

    for (i in 0..<L_new.numRows) {
        for (j in 0..<L_new.numCols) {
            L_new[i, j] = L_new[i, j] / scalevec[j]
        }
    }

    for (j in 0..<Z_new.numCols) {
        Z_new[j] = Z_new[j] / scalevec[j]
    }

    val Lsum = Matrix(L_new.numRows, 1)
    val Lmax = Matrix(L_new.numRows, 1)

    for (i in 0..<L_new.numRows) {
        val L_row_i = Matrix(1, L_new.numCols)
        Matrix.extract(L_new, i, i + 1, 0, L_new.numCols, L_row_i, 0, 0)
        Lsum[i] = L_new.sumRows(i)
        Lmax[i] = L_row_i.elementMax()
    }

    val demStations: MutableList<Int> = ArrayList()
    var L_tmp = Matrix(0, L_new.numCols)
    var mu_tmp = Matrix(0, mu_new.numCols)

    for (i in 0..<L_new.numRows) {
        if (!java.lang.Double.isNaN(Lmax[i] / Lsum[i]) && Lmax[i] / Lsum[i] > GlobalConstants.FineTol) {
            demStations.add(i)
            val L_row_i = Matrix(1, L_new.numCols)
            Matrix.extract(L_new, i, i + 1, 0, L_new.numCols, L_row_i, 0, 0)
            L_tmp = Matrix.concatRows(L_tmp, L_row_i, null)
            val mu_row_i = Matrix.extractRows(mu_new, i, i + 1, null)
            mu_tmp = Matrix.concatRows(mu_tmp, mu_row_i, null)
        }
    }
    L_new = L_tmp.copy()
    mu_new = mu_tmp.copy()

    var flag = false
    for (i in 0..<N_new.numCols) {
        if (FastMath.abs(L_new.sumCols(i) + Z_new.sumCols(i)) < GlobalConstants.FineTol && N_new[i] > GlobalConstants.FineTol) {
            flag = true
            break
        }
    }

    if (flag) {
        println("pfqn_ncld warning: The model has no positive demands in any class.")
        if (Z_new.isEmpty || Z_new.elementSum() < options.tol) {
            lG = 0.0
        } else {
            val tmp1 = Z_new.sumCols()
            val tmp2 = scalevec.copy()
            for (i in 0..<tmp1.length()) {
                tmp1[i] = FastMath.log(tmp1[i])
                tmp2[i] = FastMath.log(tmp2[i])
            }
            lG = -Matrix.factln(N_new).elementSum() + N_new.elementMult(tmp1, null)
                .elementSum() + N_new.mult(tmp2.transpose())[0]
        }
        G = Double.NaN
        return Ret.pfqnNc(G, lG, method)
    }

    val M = L_new.numRows
    R = L_new.numCols

    if (L_new.isEmpty || L_new.elementSum() < options.tol) {
        if (Z_new.isEmpty || Z_new.elementSum() < options.tol) {
            lG = 0.0
        } else {
            val tmp1 = Z_new.sumCols()
            val tmp2 = scalevec.copy()
            for (i in 0..<tmp1.length()) {
                tmp1[i] = FastMath.log(tmp1[i])
                tmp2[i] = FastMath.log(tmp2[i])
            }
            lG = (-Matrix.factln(N_new).elementSum() + N_new.elementMult(tmp1, null)
                .elementSum() + N_new.mult(tmp2.transpose())[0])
        }
        return Ret.pfqnNc(G, lG, method)
    } else if (M == 1 && (Z_new.isEmpty || Z_new.elementSum() < options.tol)) {
        val tmp1 = L_new.sumCols()
        val tmp2 = scalevec.copy()
        for (i in 0..<tmp1.length()) {
            tmp1[i] = FastMath.log(tmp1[i])
            tmp2[i] = FastMath.log(tmp2[i])
        }
        val tmp3 = mu_new.copy()
        for (i in 0..<tmp3.length()) {
            tmp3[i] = FastMath.log(tmp3[i])
        }
        lG = (Maths.factln(N_new.elementSum()) - Matrix.factln(N_new).elementSum() + N_new.elementMult(tmp1, null)
            .elementSum() + N_new.mult(tmp2.transpose())[0]) - tmp3.elementSum()
        return Ret.pfqnNc(G, lG, method)
    }

    val zeroDemandClasses: MutableList<Int> = ArrayList()
    val nonzeroDemandClasses: MutableList<Int> = ArrayList()

    for (i in 0..<R) {
        if (L_new.sumCols(i) < options.tol) {
            zeroDemandClasses.add(i)
        } else {
            nonzeroDemandClasses.add(i)
        }
    }

    val lGzdem: Double
    var Nz: Matrix
    var Zz = Matrix(Z_new.numRows, 0)
    for (i in zeroDemandClasses) {
        val Z_col_i = Matrix(Z_new.numRows, 1)
        Matrix.extract(Z_new, 0, Z_new.numRows, i, i + 1, Z_col_i, 0, 0)
        Zz = Matrix.concatColumns(Zz, Z_col_i, null)
    }

    flag = true
    for (i in 0..<Zz.numCols) {
        if (Zz.sumCols(i) >= options.tol) {
            flag = false
            break
        }
    }
    if (Z_new.isEmpty || flag) {
        lGzdem = 0.0
        Nz = Matrix(1, 1)
        Nz.fill(0.0)
    } else {
        if (zeroDemandClasses.isEmpty()) {
            lGzdem = 0.0
            Nz = Matrix(1, 1)
            Nz.fill(0.0)
        } else {
            Nz = Matrix(1, 0)
            for (i in zeroDemandClasses) {
                val N_col_i = Matrix(1, 1)
                Matrix.extract(N_new, 0, 1, i, i + 1, N_col_i, 0, 0)
                Nz = Matrix.concatColumns(Nz, N_col_i, null)
            }

            val tmp1 = Zz.sumCols()
            var tmp2 = Matrix(1, 0)
            for (i in zeroDemandClasses) {
                val scalevec_col_i = Matrix(1, 1)
                Matrix.extract(scalevec, 0, 1, i, i + 1, scalevec_col_i, 0, 0)
                tmp2 = Matrix.concatColumns(tmp2, scalevec_col_i, null)
            }
            for (i in 0..<tmp1.length()) {
                tmp1[i] = FastMath.log(tmp1[i])
                tmp2[i] = FastMath.log(tmp2[i])
            }

            lGzdem = (-Matrix.factln(Nz).elementSum() + Nz.elementMult(tmp1, null)
                .elementSum() + Nz.mult(tmp2.transpose())[0])
        }
    }

    L_tmp = Matrix(L_new.numRows, 0)
    var N_tmp = Matrix(1, 0)
    var Z_tmp = Matrix(Z_new.numRows, 0)
    var scalevecz = Matrix(1, 0)

    for (i in nonzeroDemandClasses) {
        val L_col_i = Matrix(L_new.numRows, 1)
        val N_col_i = Matrix(1, 1)
        val Z_col_i = Matrix(Z_new.numRows, 1)
        val scalevec_col_i = Matrix(1, 1)
        Matrix.extract(L_new, 0, L_new.numRows, i, i + 1, L_col_i, 0, 0)
        Matrix.extract(N_new, 0, 1, i, i + 1, N_col_i, 0, 0)
        Matrix.extract(Z_new, 0, Z_new.numRows, i, i + 1, Z_col_i, 0, 0)
        Matrix.extract(scalevec, 0, 1, i, i + 1, scalevec_col_i, 0, 0)
        L_tmp = Matrix.concatColumns(L_tmp, L_col_i, null)
        N_tmp = Matrix.concatColumns(N_tmp, N_col_i, null)
        Z_tmp = Matrix.concatColumns(Z_tmp, Z_col_i, null)
        scalevecz = Matrix.concatColumns(scalevecz, scalevec_col_i, null)
    }
    L_new = L_tmp
    N_new = N_tmp
    Z_new = Z_tmp

    val lGnnzdem: Double
    if (N_new.elementMin() < 0.0) {
        lGnnzdem = 0.0
    } else {
        val ret = compute_norm_const_ld(L_new, N_new, Z_new, mu_new, options)
        lGnnzdem = ret.lG
        method = ret.method
    }

    val tmp = scalevecz.copy()
    for (i in 0..<tmp.length()) {
        tmp[i] = FastMath.log(tmp[i])
    }
    lG = lGnnzdem + lGzdem + N_new.mult(tmp.transpose())[0]
    G = FastMath.exp(lG)
    return Ret.pfqnNc(G, lG, method)
}


/**
 * Run a normalizing constant solution method in a load-dependent model
 *
 * @param L       - demands at all stations
 * @param N       - number of jobs for each class
 * @param Z       - think time for each class
 * @param mu      - load-depedent scalings
 * @param options - solver options
 * @return normalizing constant and its logarithm
 */

fun compute_norm_const_ld(L: Matrix, N: Matrix, Z: Matrix, mu: Matrix, options: SolverOptions): Ret.pfqnNc {
    val M = L.numRows
    val R = L.numCols
    var method = options.method
    var lG: Double? = null

    when (options.method) {
        "default", "exact" -> {
            val Lz: Matrix
            val muz: Matrix
            if (Z.elementSum() < GlobalConstants.FineTol) {
                Lz = L
                muz = mu
            } else {
                val D = Z.numRows
                Lz = Matrix.concatRows(L, Z, null)
                val tmp = Matrix(1, mu.numCols)
                var i = 0
                while (i < tmp.length()) {
                    tmp[i] = (i + 1).toDouble()
                    i++
                }
                muz = Matrix.concatRows(mu, tmp.repmat(D, 1), null)
            }
            if (R == 1) {
                lG = pfqn_gldsingle(Lz, N, muz, options).lG
                method = "exact/gld"
            } else if (M == 1  && Z.elementMax() > 0) {
                val ret = pfqn_comomrm_ld(L, N, Z, muz, options)
                lG = ret.lG
                method = "exact/comomld"
            } else if (M == 2 && Z.elementMax() < GlobalConstants.FineTol) {
                val zeroZ = Matrix(N.numRows, N.numCols)
                zeroZ.fill(0.0)
                val ret = pfqn_comomrm_ld(L, N, zeroZ, mu, options)
                lG = ret.lG
                method = "exact/comomld"
            } else {
                val ret = pfqn_gld(Lz, N, muz, options)
                lG = ret.lG
                method = "exact/gld"
            }
        }

        "rd" -> {
            lG = pfqn_rd(L, N, Z, mu, options).lG
        }

        "nrp" -> {
            lG = pfqn_nrp(L, N, Z, mu, options)
        }

        "nrl" -> {
            lG = pfqn_nrl(L, N, Z, mu, options)
        }

        "comomld" -> {
            if (M <= 1 || Z.elementSum() <= GlobalConstants.Zero) {
                lG = pfqn_comomrm_ld(L, N, Z, mu, options).lG
            } else {
                line_warning(mfilename(object : Any() {}),
                    "Load-dependent CoMoM is available only in models with a delay and m identical stations, running the \"rd\" algorithm instead.\n")
                lG = pfqn_rd(L, N, Z, mu, options).lG
                method = "rd"
            }
        }

        else -> throw RuntimeException("Unrecognized method: " + options.method)
    }
    val G = FastMath.exp(lG!!)
    return Ret.pfqnNc(G, lG, method)
}
/**
 * PFQN ncld algorithms
 */
@Suppress("unused")
class PfqnNcldAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}