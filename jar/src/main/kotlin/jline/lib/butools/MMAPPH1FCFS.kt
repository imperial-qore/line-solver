package jline.lib.butools

import jline.api.mc.ctmc_solve
import jline.util.Maths
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import org.apache.commons.math3.util.CombinatoricsUtils
import org.apache.commons.math3.util.FastMath
import kotlin.math.pow

fun MMAPPH1FCFS(D: MatrixCell,
                sigma: Map<Int?, Matrix>,
                S: Map<Int?, Matrix>,
                numOfQLMoms: Int?,
                numOfQLProbs: Int?,
                numOfSTMoms: Int?,
                stDistr: Matrix?,
                stDistrME: Boolean,
                stDistrPH: Boolean,
                prec: Double?,
                classes_: Matrix?): Map<String, MutableMap<Int, Matrix>> {
    val K = D.size() - 1
    var precision = 1e-14
    var classes = Matrix(1, K, K)
    for (i in 0..<K) {
        classes[i] = i.toDouble()
    }
    if (prec != null) {
        precision = prec
    }

    if (classes_ != null) {
        classes = classes_
    }

    val D0 = D[0]
    val N = D0.numRows
    val Ia = Matrix.eye(N)
    var Da = Matrix(N, N, N * N)
    for (q in 0..<K) {
        Da = Da.add(1.0, D[q + 1])
    }
    val beta: MutableMap<Int, Matrix> = HashMap()
    val theta = ctmc_solve(D0.add(1.0, Da))
    val lambda = Matrix(K, 1, K)
    val mu = Matrix(K, 1, K)
    val Nsk = Matrix(1, K, K)
    var ro = 0.0
    for (k in 0..<K) {
        lambda[k] = theta.mult(D[k + 1]).elementSum()
        beta[k] = ctmc_solve(S[k]!!.add(-1.0, S[k]!!.sumRows().mult(sigma[k]!!)))
        val neg_sk = S[k]!!.copy()
        neg_sk.scaleEq(-1.0)
        mu[k] = beta[k]!!.mult(neg_sk).elementSum()
        Nsk[k] = S[k]!!.numRows.toDouble()
        ro = ro + lambda[k] / mu[k]
    }

    val alpha = theta.mult(Da)
    alpha.scaleEq(1 / lambda.elementSum())
    val D0i = D0.inv().copy()
    D0i.scaleEq(-1.0)

    var Sa = S[0]
    val sa: MutableMap<Int, Matrix?> = HashMap()
    val ba: MutableMap<Int, Matrix?> = HashMap()
    val sv: MutableMap<Int, Matrix> = HashMap()
    sa[0] = sigma[0]
    ba[0] = beta[0]
    val sv0 = S[0]!!.sumRows()
    sv0.scaleEq(-1.0)
    sv[0] = sv0

    val Pk: MutableMap<Int, Matrix> = HashMap()
    Pk[0] = D0i.mult(D[1])

    for (q in 1..<K) {
        sa[q] = Matrix(1, sigma[0]!!.length(), sigma[0]!!.length())
        ba[q] = Matrix(1, beta[0]!!.length(), beta[0]!!.length())
        sv[q] = Matrix(sigma[0]!!.length(), 1, sigma[0]!!.length())
        Pk[q] = D0i.mult(D[q + 1])
    }

    for (k in 1..<K) {
        Sa = Sa!!.createBlockDiagonal(S[k])
        for (q in 0..<K) {
            if (q == k) {
                sa[q] = Matrix.concatColumns(sa[q]!!, sigma[k]!!, null)
                ba[q] = Matrix.concatColumns(ba[q]!!, beta[k]!!, null)
                val sk_neg_sum = S[k]!!.sumRows()
                sk_neg_sum.scaleEq(-1.0)
                sv[q] = Matrix.concatRows(sv[q]!!, sk_neg_sum, null)
            } else {
                sa[q] = Matrix.concatColumns(sa[q]!!, Matrix(sigma[k]!!.numRows, sigma[k]!!.numCols, 0), null)
                ba[q] = Matrix.concatColumns(ba[q]!!, Matrix(beta[k]!!.numRows, beta[k]!!.numCols, 0), null)
                sv[q] = Matrix.concatRows(sv[q]!!, Matrix(sigma[k]!!.length(), 1, 0), null)
            }
        }
    }
    D0i.mult(Da)
    var iVec = D[1].kron(sa[0]!!)

    for (k in 1..<K) {
        iVec = iVec.add(1.0, D[k + 1].kron(sa[k]!!))
    }

    val Ns = Sa!!.numRows
    val Is = Matrix.eye(Ns)
    val neg_Sa_row_sum = Sa.sumRows()
    neg_Sa_row_sum.scaleEq(-1.0)
    val Y0 = FluidFundamentalMatrices(Ia.kron(Sa),
        Ia.kron(neg_Sa_row_sum),
        iVec,
        D0,
        precision,
        null,
        null)["P"]
    val T = Ia.kron(Sa).add(1.0, Y0!!.mult(iVec))
    var pi0 = Matrix(1, T.numRows, T.numRows)
    for (k in 0..<K) {
        val ba_mu = ba[k]!!.copy()
        ba_mu.scaleEq(1 / mu[k])
        pi0 = pi0.add(1.0, theta.mult(D[k + 1]).kron(ba_mu))
    }
    pi0 = pi0.mult(T)
    pi0.scaleEq(-1.0)

    val iT = T.inv()
    iT.scaleEq(-1.0)
    val oa = Matrix.ones(N, 1)

    val result: MutableMap<String, MutableMap<Int, Matrix>> = HashMap()
    if (numOfSTMoms != null) {
        result["stNoms"] = HashMap()
    }
    if (stDistr != null) {
        result["stDistr"] = HashMap()
    }
    if (stDistrME) {
        result["stDistrME_alpha"] = HashMap()
        result["stDistrME_A"] = HashMap()
    }
    if (stDistrPH) {
        result["stDistrPH_alpha"] = HashMap()
        result["stDistrPH_A"] = HashMap()
    }

    if (numOfQLProbs != null) {
        result["ncDistr"] = HashMap()
    }

    if (numOfQLMoms != null) {
        result["ncMoms"] = HashMap()
    }


    for (i in 0..<classes.length()) {
        val k = classes[i].toInt()
        val clo = iT.mult(oa.kron(sv[k]!!))
        if (numOfSTMoms != null) {
            val rtMoms = Matrix(1, numOfSTMoms)
            for (m in 0..<numOfSTMoms) {
                val rtMoms_m = pi0.mult(Matrix.pow(iT, m)).mult(clo)[0] / pi0.mult(clo)[0]
                rtMoms[m] = CombinatoricsUtils.factorial(m + 1) * rtMoms_m
            }
            result["stNoms"]!![k] = rtMoms
        }

        if (stDistr != null) {
            var cdf = Matrix(0, 0, 0)
            for (p in 0..<stDistr.length()) {
                val Tt = T.copy()
                Tt.scaleEq(stDistr[p])
                val pr = 1 - pi0.mult(Maths.matrixExp(Tt)).mult(clo)[0] / pi0.mult(clo)[0]
                val cdf_ = Matrix(1, 1, 1)
                cdf_[0, 0] = pr
                cdf = if (p == 0) {
                    cdf_
                } else {
                    Matrix.concatColumns(cdf, cdf_, null)
                }
            }
            result["stDistr"]!![k] = cdf
        }

        if (stDistrME) {
            val Bm = SimilarityMatrixForVectors(clo.mult(pi0.mult(clo).inv()), Matrix.ones(N * Ns, 1))
            val Bmi = Bm.inv()
            val A = Bm.mult(T).mult(Bmi)
            pi0.mult(Bmi)
            result["stDistrME_alpha"]!![k] = alpha
            result["stDistrME_A"]!![k] = A
        }


        if (stDistrPH) {
            val vv = pi0.mult(iT)

            val nz: MutableList<Double> = ArrayList()
            val nz_index: MutableList<Int> = ArrayList()
            for (n in 0..<vv.length()) {
                if (vv[n] > precision) {
                    nz.add(vv[n])
                    nz_index.add(n)
                }
            }
            val delta = Matrix.diag(*nz.stream().mapToDouble { obj: Double -> obj.toDouble() }.toArray())
            val neg_T = T.copy()
            neg_T.scaleEq(-1.0)
            val cl = neg_T.mult(clo)
            cl.scaleEq(pi0.mult(clo).value())
            var alpha_stDistrPH = Matrix(0, 0, 0)
            for (n in nz_index.indices) {
                if (alpha_stDistrPH.length() == 0) {
                    alpha_stDistrPH = Matrix.extractRows(cl, nz_index[n], nz_index[n] + 1, null)
                } else {
                    Matrix.concatRows(alpha_stDistrPH, Matrix.extractRows(cl, nz_index[n], nz_index[n] + 1, null), null)
                }
            }
            var A = Matrix(nz_index.size, nz_index.size, nz_index.size.toDouble().pow(2.0).toInt())
            for (n in nz_index.indices) {
                for (m in nz_index.indices) {
                    A[nz_index[n], nz_index[m]] = T[nz_index[n], nz_index[m]]
                }
            }
            A = delta.inv().mult(A.transpose()).mult(delta)
            result["stDistrPH_alpha"]!![k] = alpha
            result["stDistrPH_A"]!![k] = A
        }

        if (numOfQLProbs != null) {
            val value = Matrix(1, numOfQLProbs, numOfQLProbs)
            val jm = Matrix(Ns, 1, Ns)
            run {
                var n = Nsk.sumRows(0, k).elementSum().toInt()
                while (n < Nsk.sumRows(0, k + 1).elementSum()) {
                    jm[n, 0] = 1.0
                    n++
                }
            }
            var jmc = Matrix.ones(Ns, 1)
            jmc = jmc.add(-1.0, jm)
            var LmCurr = Matrix.lyap(T, D0.add(1.0, Da).add(-1.0, D[k + 1]).kron(Is), Matrix.eye(N * Ns), null)
            value[0] = 1 - ro + pi0.mult(LmCurr).mult(oa.kron(jmc))[0]
            for (n in 0..<numOfQLProbs - 1) {
                val LmPrev = LmCurr.copy()
                LmCurr =
                    Matrix.lyap(T, D0.add(1.0, Da).add(-1.0, D[k + 1]).kron(Is), LmPrev.mult(D[k + 1].kron(Is)), null)
                value[n + 1] = pi0.mult(LmCurr).mult(oa.kron(jmc))[0] + pi0.mult(LmPrev).mult(oa.kron(jm))[0]
            }
            result["ncDistr"]!![k] = value
        }

        if (numOfQLMoms != null) {
            val jm = Matrix(Ns, 1, Ns)
            run {
                var n = Nsk.sumRows(0, k).elementSum().toInt()
                while (n < Nsk.sumRows(0, k + 1).elementSum()) {
                    jm[n, 0] = 1.0
                    n++
                }
            }
            val ELn: MutableMap<Int, Matrix> = HashMap()
            ELn[0] = Matrix.lyap(T, D0.add(1.0, Da).kron(Is), Matrix.eye(N * Ns), null)
            val qlMoms = Matrix(1, numOfQLMoms, numOfQLMoms)
            for (n in 0..<numOfQLMoms) {
                var bino = 1
                var Btag = Matrix(N * Ns, N * Ns, FastMath.pow((N * Ns).toDouble(), 2).toInt())
                for (m in -1..n - 1) {
                    Btag = Btag.add(bino.toDouble(), ELn[m + 1]!!)
                    bino = bino * (n - m) / (m + 2)
                }
                ELn[n + 1] = Matrix.lyap(T, D0.add(1.0, Da).kron(Is), Btag.mult(D[k + 1].kron(Is)), null)
                qlMoms[n] = pi0.mult(ELn[n + 1]!!).elementSum() + pi0.mult(Btag).mult(oa.kron(jm))[0]
            }
            result["ncMoms"]!![k] = qlMoms
        }
    }


    return result
}