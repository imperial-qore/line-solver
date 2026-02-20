/**
 * Normalizing Constant Methods for Product-Form Networks
 * 
 * Provides a comprehensive suite of normalizing constant algorithms including convolution,
 * Laguerre expansion, cubature, and sampling methods. Automatically selects the most
 * appropriate algorithm based on network characteristics and computational constraints.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc

import jline.io.Ret
import jline.io.line_warning
import jline.io.mfilename
import jline.GlobalConstants
import jline.GlobalConstants.NegInf
import jline.VerboseLevel
import jline.solvers.SolverOptions
import jline.util.Maths
import jline.util.Maths.binomialCoeff
import jline.util.Utils
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath
import kotlin.math.ceil
import kotlin.math.min
import kotlin.math.pow

/**
 * Main method to compute the normalizing constant of a load-independent model
 *
 * @param lambda  - arrival rate of open classes
 * @param L       - demands at all stations
 * @param N       - number of jobs for each class
 * @param Z       - think time for each class
 * @param options - solver options
 * @return normalizing constant, its logarithm, and mean performance metrics computed as a by-product
 */

fun pfqn_nc(lambda: Matrix, L: Matrix, N: Matrix, Z: Matrix, options: SolverOptions): Ret.pfqnNcXQ {/*
		SolverOptions options = Solver.parseOptions();
     */
    //L.print(); N.print(); Z.print();
    var method = "exact" //if early return is triggered
    N.length()

    var lG: Double? = null
    var X = Matrix(0, 0)
    var Q = Matrix(0, 0)
    var lambda_new = lambda.copy()
    var L_new = L.copy()
    var N_new = N.copy()
    var Z_new = Z.copy()

    // Check if any element in N is negative
    var hasNegative = false
    for (i in 0 until N_new.length()) {
        if (N_new.get(i) < 0) {
            hasNegative = true
            break
        }
    }

    if (hasNegative || N_new.isEmpty) {
        lG = NegInf
        return Ret.pfqnNcXQ(lG, X, Q, method)
    }

    if (N_new.elementSum() < GlobalConstants.FineTol) {
        lG = 0.0
        return Ret.pfqnNcXQ(lG, X, Q, method)
    }

    if (lambda_new.isEmpty) {
        lambda_new = N_new.copy()
        lambda_new.fill(0.0)
    }

    var Qopen = Matrix(0, lambda_new.length())
    val Ut = Matrix(1, L_new.numRows)
    val lGopen = 0.0
    for (i in 0..<L_new.numRows) {
        val L_row_i = Matrix(1, L_new.numCols)
        Matrix.extract(L_new, i, i + 1, 0, L_new.numCols, L_row_i, 0, 0)
        Ut[i] = 1 - lambda_new.mult(L_row_i.transpose())[0]
        if (java.lang.Double.isNaN(Ut[i])) {
            Ut[i] = 0.0
        }
        for (j in 0..<L_row_i.length()) {
            L_new[i, j] = L_new[i, j] / Ut[i]
        }
        Matrix.extract(L_new, i, i + 1, 0, L_new.numCols, L_row_i, 0, 0)
        for (j in 0..<L_row_i.length()) {
            L_row_i[j] = L_row_i[j] / Ut[i]
        }
        val tmp = lambda_new.elementMult(L_row_i, null)
        Qopen = Matrix.concatRows(Qopen, tmp, null)
    }
    Qopen.removeNaN()

    val ocl: MutableList<Int> = ArrayList()
    for (i in 0..<N_new.length()) {
        if (Utils.isInf(N_new[i])) {
            ocl.add(i)
        }
    }

    for (i in 0..<N_new.length()) {
        if (Utils.isInf(N_new[i])) {
            N_new[i] = 0.0
        }
    }

    var L_tmp = Matrix(L_new.numRows, 0)
    var N_tmp = Matrix(N_new.numRows, 0)
    var Z_tmp = Matrix(Z_new.numRows, 0)
    var lambda_tmp: Matrix? = Matrix(lambda_new.numRows, 0)
    for (i in 0..<N_new.length()) {
        if (FastMath.abs(N_new[i]) >= GlobalConstants.FineTol) {
            val L_col = Matrix.extractColumn(L_new, i, null)
            val N_col = Matrix.extractColumn(N_new, i, null)
            val Z_col = Matrix.extractColumn(Z_new, i, null)
            val lambda_col = Matrix.extractColumn(lambda_new, i, null)
            L_tmp = Matrix.concatColumns(L_tmp, L_col, null)
            N_tmp = Matrix.concatColumns(N_tmp, N_col, null)
            Z_tmp = Matrix.concatColumns(Z_tmp, Z_col, null)
            lambda_tmp = Matrix.concatColumns(lambda_tmp, lambda_col, null)
        }
    }
    L_new = L_tmp
    N_new = N_tmp
    Z_new = Z_tmp

    var R = N_new.length()
    val scalevec = Matrix(1, R)
    scalevec.fill(1.0)
    for (r in 0..<R) {
        val L_col_r = Matrix(L_new.numRows, 1)
        if (L_new.numCols > 0) {
            Matrix.extractColumn(L_new, r, L_col_r)
        } else {
            L_col_r.fill(0.0)
        }
        val Z_col_r = Matrix(Z_new.numRows, 1)
        if (Z_new.numCols > 0) {
            Matrix.extractColumn(Z_new, r, Z_col_r)
        } else {
            Z_col_r.fill(0.0)
        }
        scalevec[r] = FastMath.max(L_col_r.elementMax(), Z_col_r.elementMax())
    }

    for (i in 0..<L_new.numRows) {
        for (j in 0..<L_new.numCols) {
            L_new[i, j] = L_new[i, j] / scalevec[j]
        }
    }

    for (i in 0..<Z_new.numCols) {
        Z_new[i] = Z_new[i] / scalevec[i]
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
    val noDemStations: MutableList<Int> = ArrayList()
    L_tmp = Matrix(0, L_new.numCols)

    for (i in 0..<L_new.numRows) {
        if (!java.lang.Double.isNaN(Lmax[i] / Lsum[i]) && Lmax[i] / Lsum[i] > GlobalConstants.FineTol) {
            demStations.add(i)
            val L_row_i = Matrix(1, L_new.numCols)
            Matrix.extract(L_new, i, i + 1, 0, L_new.numCols, L_row_i, 0, 0)
            L_tmp = Matrix.concatRows(L_tmp, L_row_i, null)
        } else {
            noDemStations.add(i)
        }
    }
    L_new = L_tmp

    var flag = false
    for (i in 0..<N_new.numCols) {
        if (FastMath.abs(L_new.sumCols(i) + Z_new.sumCols(i)) < GlobalConstants.FineTol && N_new[i] > GlobalConstants.FineTol) {
            flag = true
            break
        }
    }

    if (flag) {
        if (options.verbose != VerboseLevel.SILENT) {
            line_warning(mfilename(object : Any() {}),
                "pfqn_nc warning: The model has no positive demands in any class.")
        }
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
        return Ret.pfqnNcXQ(lG, X, Q, method)
    }

    val M = L_new.numRows
    R = L_new.numCols

    if (L_new.isEmpty || L_new.elementSum() < options.tol) {
        if (Z_new.isEmpty || Z_new.elementSum() < options.tol) {
            lG = lGopen
        } else {
            val tmp1 = Z_new.sumCols()
            val tmp2 = scalevec.copy()
            for (i in 0..<tmp1.length()) {
                tmp1[i] = FastMath.log(tmp1[i])
                tmp2[i] = FastMath.log(tmp2[i])
            }
            lG = (lGopen - Matrix.factln(N_new).elementSum() + N_new.elementMult(tmp1, null).elementSum() + N_new.mult(
                tmp2.transpose())[0])
        }
        return Ret.pfqnNcXQ(lG, X, Q, method)
    } else if (M == 1 && (Z_new.isEmpty || Z_new.elementSum() < options.tol)) {
        val tmp1 = L_new.sumCols()
        val tmp2 = scalevec.copy()
        for (i in 0..<tmp1.length()) {
            tmp1[i] = FastMath.log(tmp1[i])
            tmp2[i] = FastMath.log(tmp2[i])
        }
        lG = (Maths.factln(N_new.elementSum()) - Matrix.factln(N_new).elementSum() + N_new.elementMult(tmp1, null)
            .elementSum() + N_new.mult(tmp2.transpose())[0])
        return Ret.pfqnNcXQ(lG, X, Q, method)
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
    N_tmp = Matrix(1, 0)
    Z_tmp = Matrix(Z_new.numRows, 0)
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

    val ret = compute_norm_const(L_new, N_new, Z_new, options)
    val lGnzdem = ret.lG
    val Xnnzdem = ret.X
    ret.Q
    method = ret.method

    if (Xnnzdem.isEmpty) {
        X = Matrix(0, 0)
        Q = Matrix(0, 0)
    }

    val tmp = scalevecz.copy()
    for (i in 0..<tmp.length()) {
        tmp[i] = FastMath.log(tmp[i])
    }
    lG = lGopen + lGnzdem + lGzdem + N_new.mult(tmp.transpose())[0]

    //System.out.println(lG);
    return Ret.pfqnNcXQ(lG, X, Q, method)
}


/**
 * Run a normalizing constant solution method in a load-independent model
 *
 * @param L       - demands at all stations
 * @param N       - number of jobs for each class
 * @param Z       - think time for each class
 * @param options - solver options
 * @return normalizing constant, its logarithm, and mean performance metrics computed as a by-product
 */

fun compute_norm_const(L: Matrix, N: Matrix, Z: Matrix, options: SolverOptions): Ret.pfqnNcXQ {
    val M = L.numRows
    val R = L.numCols
    val X = Matrix(0, 0)
    val Q = Matrix(0, 0)
    var method = options.method
    var lG: Double? = null

    when (options.method) {
        "ca" -> {
            val ret = pfqn_ca(L, N, Z.sumCols())
            lG = ret.lG
        }

        "default", "adaptive" -> {
            val Z_colSum = Z.sumCols()

            if (M > 1) {/* ------------------------------------------
                     *  multi-station model
                     *  – small population  ⇒ Grundmann-Mueller
                     *  – otherwise         ⇒ Laguerre expansion
                     * ------------------------------------------ */
                if (N.elementSum() < 1e3) {
                    val Cmax = M * R * 50.0.pow(3) // M stations, R classes, (50)^3 budget
                    // maximum possible cubature order
                    val maxOrder = min(ceil((N.elementSum() - 1) / 2.0).toInt(), 16)

                    var totCost = 0.0
                    var order = 0

                    // raise order as far as the cost‐budget allows
                    while (order < maxOrder) {
                        val nextCost = R * binomialCoeff(M + 2 * (order + 1), M - 1)
                        if (totCost + nextCost <= Cmax) {
                            order += 1
                            totCost += nextCost
                        } else {
                            break
                        }
                    }

                    val ret = pfqn_cub(L, N, Z_colSum, order, GlobalConstants.FineTol)
                    lG = ret.lG
                    method = "cub"
                } else {
                    val ret = pfqn_le(L, N, Z_colSum)
                    lG = ret.lG
                    method = "le"
                }
            } else {   /* ---------- single-station model ---------- */
                if (Z_colSum.numCols == 1 && FastMath.abs(Z_colSum[0]) < GlobalConstants.FineTol) {/* M = 1, no delay → exact */

                    val tmp = L.copy()
                    var i = 0
                    while (i < tmp.length()) {
                        tmp[i] = FastMath.log(tmp[i])
                        i++
                    }
                    lG = -N.mult(tmp.transpose())[0]
                    method = "exact"
                } else {        /* repairman */
                    if (N.elementSum() < 10000) {
                        val ret = pfqn_comomrm(L, N, Z_colSum, 1, GlobalConstants.Zero)
                        lG = ret.lG
                        method = "comom"
                    } else {
                        val ret = pfqn_le(L, N, Z_colSum)
                        lG = ret.lG
                        method = "le"
                    }
                }
            }
        }

        "sampling" -> {
            if (M == 1) {
                val ret = pfqn_mmsample2(L, N, Z.sumCols(), options.samples)
                lG = ret.lG
                method = "sampling"
            } else if (M > R) {
                val ret = pfqn_mci(L, N, Z.sumCols(), options.samples, "imci")
                lG = ret.lG
                method = "imci"
            } else {
                val ret = pfqn_ls(L, N, Z.sumCols(), options.samples.toLong(), options.seed.toLong())
                lG = ret.lG
                method = "ls"
            }
        }

        "cub", "gm" -> {
            val order = FastMath.ceil((N.elementSum() - 1) / 2).toInt() // exact
            val ret = pfqn_cub(L, N, Z.sumCols(), order, GlobalConstants.FineTol)
            lG = ret.lG
        }

        "kt" -> {
            val ret = pfqn_kt(L, N, Z.sumCols())
            lG = ret.lG
            method = "kt"
        }

        "mmint2", "gleint" -> {
            if (L.numRows > 1) {
                throw RuntimeException("The " + options.method + " method requires a model with a delay and a single queueing station.")
            } else {
                val ret = pfqn_mmint2_gausslegendre(L, N, Z.sumCols(), null)
                lG = ret.lG
            }
        }

        "le" -> {
            val ret = pfqn_le(L, N, Z.sumCols())
            lG = ret.lG
        }

        "ls" -> {
            val ret = pfqn_ls(L, N, Z.sumCols(), options.samples.toLong(), options.seed.toLong())
            lG = ret.lG
        }

        "mci", "imci" -> {
            val ret = pfqn_mci(L, N, Z.sumCols(), options.samples, options.method)
            lG = ret.lG
        }

        "exact" -> {
            if (M >= R || N.elementSum() > 10 || Z.elementSum() > 0) {
                val ret = pfqn_ca(L, N, Z.sumCols())
                lG = ret.lG
                method = "exact/ca"
            } else {
                val ret = pfqn_recal(L, N, Z.sumCols())
                lG = ret.lG
                method = "exact/recal"
            }
        }

        "comom" -> {
            if (R > 1) {
                try {
                    if (M <= 1) {
                        val ret = pfqn_comomrm(L, N, Z, 1, options.tol)
                        lG = ret.lG
                    }
                } catch (e: Exception) {
                    e.printStackTrace()
                    lG = Double.NaN
                }
            } else {
                val ret = pfqn_ca(L, N, Z.sumCols())
                lG = ret.lG
                method = "ca"
            }
        }

        "comomld" -> {
            // For comomld method, fall back to comom if applicable, otherwise use ca
            if (R > 1 && M <= 1) {
                try {
                    val ret = pfqn_comomrm(L, N, Z, 1, options.tol)
                    lG = ret.lG
                    method = "comom"
                } catch (e: Exception) {
                    e.printStackTrace()
                    val ret = pfqn_ca(L, N, Z.sumCols())
                    lG = ret.lG
                    method = "ca"
                }
            } else {
                val ret = pfqn_ca(L, N, Z.sumCols())
                lG = ret.lG
                method = "ca"
            }
        }

        else -> {
            line_warning("pfqn_nc", "unrecognized method \"%s\"", options.method)
            lG = Double.NaN // indicate failure but keep running
            method = options.method // echo back the unsupported name
        }
    }
    return Ret.pfqnNcXQ(lG, X, Q, method)
}
/**
 * PFQN nc algorithms
 */
@Suppress("unused")
class PfqnNcAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}