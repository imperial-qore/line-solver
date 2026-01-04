package jline.solvers.nc.handlers

import jline.api.npfqn.npfqn_nonexp_approx
import jline.api.pfqn.nc.pfqn_nc
import jline.api.sn.snDeaggregateChainResults
import jline.api.sn.snGetDemandsChain
import jline.lang.NetworkStruct
import jline.GlobalConstants
import jline.lang.constant.SchedStrategy
import jline.VerboseLevel
import jline.solvers.SolverOptions
import jline.solvers.nc.SolverNC
import jline.util.Utils
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath
import kotlin.math.exp

fun solver_nc(sn: NetworkStruct, options: SolverOptions): SolverNC.SolverNCReturn {
    val M = sn.nstations
    val nservers = sn.nservers
    val NK = sn.njobs.transpose()
    val sched = sn.sched
    var C = sn.nchains
    val SCV = sn.scv
    val K = sn.nclasses
    val startTime = System.nanoTime()

    // Check for LCFS scheduling - not supported in this version
    for (i in 0..<M) {
        when (sched[sn.stations[i]]) {
            SchedStrategy.LCFS, SchedStrategy.LCFSPR ->
                throw RuntimeException("LCFS queueing networks are not supported in this version.")
            else -> {}
        }
    }

    var V = Matrix(sn.nstateful, K)
    for (i in 0..<sn.visits.size) {
        V = V.add(1.0, sn.visits.get(i))
    }

    //V = cellsum(sn.visits);
    val rates = sn.rates
    var ST = rates.copy()
    for (i in 0..<ST.numRows) {
        for (j in 0..<ST.numCols) {
            ST.set(i, j, 1.0 / ST.get(i, j))
        }
    }
    ST.removeNaN()
    val ST0 = ST.copy()
    var Nchain = Matrix(1, C)
    Nchain.fill(0.0)

    for (c in 0..<C) {
        val inchain: Matrix = sn.inchain.get(c)!!
        var sum = 0.0
        for (col in 0..<sn.inchain.get(c)!!.numCols) {
            sum += NK.get(inchain.get(0, col).toInt())
        }
        Nchain.set(c, sum)
    }

    val openChains: MutableList<Int?> = ArrayList<Int?>()
    val closedChains: MutableList<Int?> = ArrayList<Int?>()

    for (c in 0..<C) {
        if (Utils.isInf(Nchain.get(c))) {
            openChains.add(c)
        } else {
            closedChains.add(c)
        }
    }

    var gamma = Matrix(1, M)
    gamma.fill(0.0)
    var eta_1 = Matrix(1, M)
    eta_1.fill(0.0)
    var eta = Matrix(1, M)
    eta.fill(1.0)
    C = sn.nchains

    if (!sched.containsValue(SchedStrategy.FCFS)) {
        options.iter_max = 1
    }

    var it = 0
    val tmp_eta = Matrix(1, M)
    for (i in 0..<M) {
        tmp_eta.set(i, FastMath.abs(1 - eta.get(i) / eta_1.get(i)))
    }

    var lambda: Matrix? = null
    var Lchain: Matrix? = null
    var STchain: Matrix? = null
    var Vchain: Matrix? = null
    var alpha: Matrix? = null
    var Q: Matrix? = null
    var U: Matrix? = null
    var R: Matrix? = null
    var T: Matrix? = null
    var X: Matrix? = null
    var STeff: Matrix? = null
    var lG: Double? = null
    var method: String? = null

    while (tmp_eta.elementMax() > options.iter_tol && it < options.iter_max) {
        it++
        eta_1 = eta

        if (it == 1) {
            lambda = Matrix(1, C)
            lambda.fill(0.0)
            val ret = snGetDemandsChain(sn)
            Lchain = ret.Dchain
            STchain = ret.STchain
            Vchain = ret.Vchain
            alpha = ret.alpha
            Nchain = ret.Nchain

            for (c in 0..<C) {
                val inchain: Matrix = sn.inchain.get(c)!!
                var isOpenChain = false
                for (col in 0..<inchain.numCols) {
                    if (Utils.isInf(NK.get(inchain.get(0, col).toInt()))) {
                        isOpenChain = true
                        break
                    }
                }
                for (i in 0..<M) {
                    if (isOpenChain && FastMath.abs(i - sn.refstat.get(inchain.get(0).toInt())) < 1e-6) {
                        lambda.set(c, 1.0 / STchain.get(i, c))
                    }
                }
            }
        } else {
            for (c in 0..<C) {
                val inchain: Matrix = sn.inchain.get(c)!!
                for (i in 0..<M) {
                    val ST_tmp = Matrix(1, inchain.numCols)
                    val alpha_tmp = Matrix(1, inchain.numCols)
                    for (col in 0..<inchain.numCols) {
                        ST_tmp.set(col, ST.get(i, inchain.get(0, col).toInt()))
                        alpha_tmp.set(col, alpha!!.get(i, inchain.get(0, col).toInt()))
                    }
                    STchain!!.set(i, c, ST_tmp.mult(alpha_tmp.transpose()).get(0))
                    Lchain!!.set(i, c, Vchain!!.get(i, c) * STchain.get(i, c))
                }
            }
        }

        STchain!!.removeInfinity()
        Lchain!!.removeInfinity()
        val Lms = Matrix(M, C)
        val Z = Matrix(M, C)
        val Zms = Matrix(M, C)
        Lms.fill(0.0)
        Z.fill(0.0)
        Zms.fill(0.0)

        val infServers: MutableList<Int?> = ArrayList<Int?>()
        for (i in 0..<M) {
            if (Utils.isInf(nservers.get(i))) {
                infServers.add(i)
                for (j in 0..<C) {
                    Z.set(i, j, Lchain.get(i, j))
                }
            } else {
                if (options.method != null && options.method.equals("exact",
                        ignoreCase = true) && nservers.get(i) > 1) {
                    if (options.verbose != VerboseLevel.SILENT) println("SolverNC does not support exact multiserver yet. Switching to approximate method.")
                }
                for (j in 0..<C) {
                    Lms.set(i, j, Lchain.get(i, j) / nservers.get(i))
                    Zms.set(i, j, Lchain.get(i, j) * (nservers.get(i) - 1) / nservers.get(i))
                }
            }
        }

        val Z_new = Matrix(1, C)
        for (i in 0..<C) {
            Z_new.set(i, Z.sumCols(i) + Zms.sumCols(i))
        }
        val Z_tmp_append_0 = Matrix(1, Z_new.numCols + 1)
        Z_tmp_append_0.set(Z_new.numCols, 0.0)
        Matrix.extract(Z_new, 0, 1, 0, Z_new.length(), Z_tmp_append_0, 0, 0)

        val ret = pfqn_nc(lambda!!, Lms, Nchain, Z_new, options)
        lG = ret.lG
        var Xchain = ret.X
        var Qchain = ret.Q
        method = ret.method

        if (Zms.sumCols().elementMin() > GlobalConstants.FineTol) {
            Xchain = Matrix(0, 0)
            Qchain = Matrix(0, 0)
        }

        if (lG == null) {
            val runtime = (System.nanoTime() - startTime) / 1000000000.0
            return SolverNC.SolverNCReturn(Q, U, R, T, C, X, java.lang.Double.NaN, STeff, it, runtime, method)
        }

        if (Xchain.isEmpty) {
            Xchain = lambda.copy()
            Qchain = Matrix(M, C)
            Qchain.fill(0.0)
            for (r in closedChains) {
                val Nchain_tmp = Nchain.copy()
                Nchain_tmp.set(r!!, Nchain_tmp.get(r) - 1)
                val Nchain_tmp_append_1 = Matrix(1, Nchain.numCols + 1)
                Nchain_tmp_append_1.set(Nchain.numCols, 1.0)
                Matrix.extract(Nchain_tmp, 0, 1, 0, Nchain_tmp.length(), Nchain_tmp_append_1, 0, 0)
                Xchain.set(r, FastMath.exp(pfqn_nc(lambda, Lms, Nchain_tmp, Z_new, options).lG - lG))
                for (i in 0..<M) {
                    if (Lchain.get(i, r) > 1e-6) {
                        if (Utils.isInf(nservers.get(i))) {
                            Qchain.set(i, r, Lchain.get(i, r) * Xchain.get(r))
                        } else {
                            val lambda_tmp = Matrix(1, lambda.length() + 1)
                            lambda_tmp.fill(0.0)
                            Matrix.extract(lambda, 0, 1, 0, lambda.length(), lambda_tmp, 0, 0)

                            var L_tmp = Matrix(0, Lms.numCols + 1)
                            for (row in 0..<Lms.numRows) {
                                if (row != i) {
                                    val L_tmp_row = Matrix(1, Lms.numCols + 1)
                                    L_tmp_row.set(Lms.numCols, 0.0)
                                    Matrix.extract(Lms, row, row + 1, 0, Lms.numCols, L_tmp_row, 0, 0)
                                    L_tmp = Matrix.concatRows(L_tmp, L_tmp_row, null)
                                }
                            }
                            val L_tmp_row = Matrix(1, Lms.numCols + 1)
                            L_tmp_row.set(Lms.numCols, 1.0)
                            Matrix.extract(Lms, i, i + 1, 0, Lms.numCols, L_tmp_row, 0, 0)
                            L_tmp = Matrix.concatRows(L_tmp, L_tmp_row, null)

                            val ret_tmp = pfqn_nc(lambda_tmp, L_tmp, Nchain_tmp_append_1, Z_tmp_append_0, options)
                            val res = ret_tmp.lG
                            method = ret_tmp.method
                            Qchain.set(i, r, Zms.get(i, r) * Xchain.get(r) + Lms.get(i, r) * exp(res - lG))
                        }
                    }
                }
                Qchain.removeNaN()
            }

            for (r in openChains) {
                for (i in 0..<M) {
                    val lambda_open = Matrix(1, openChains.size)
                    val Lchain_i_open = Matrix(1, openChains.size)
                    val Qchain_i_closed = Matrix(1, closedChains.size)
                    for (j in openChains.indices) {
                        lambda_open.set(j, lambda.get(openChains.get(j)!!))
                        Lchain_i_open.set(j, Lchain.get(i, openChains.get(j)!!))
                    }
                    for (j in closedChains.indices) {
                        Qchain_i_closed.set(j, Qchain.get(i, closedChains.get(j)!!))
                    }
                    Qchain.set(i,
                        r!!,
                        lambda.get(r) * Lchain.get(i, r) / (1 - lambda_open.mult(Lchain_i_open.transpose())
                            .get(0) / nservers.get(i)) * (1 + Qchain_i_closed.elementSum()))
                }
            }
        } else {
            for (r in 0..<C) {
                for (i in 0..<M) {
                    if (Lchain.get(i, r) > 1e-6) {
                        if (Utils.isInf(nservers.get(i))) {
                            Qchain.set(i, r, Lchain.get(i, r) * Xchain.get(r))
                        }
                    }
                }
            }
        }

        var allNaN = true
        for (i in 0..<Xchain.length()) {
            if (!java.lang.Double.isNaN(Xchain.get(i))) {
                allNaN = false
                break
            }
        }

        if (allNaN) {
            if (options.verbose != VerboseLevel.SILENT) println("Normalizing constant computations produced a floating-point range exception. Model is likely too large.")
        }

        //Z = sum(Z(1:M,:),1);
        val Rchain = Matrix(Qchain.numRows, Qchain.numCols)
        for (i in 0..<Rchain.numRows) {
            for (j in 0..<Rchain.numCols) {
                if (Qchain.get(i, j) < GlobalConstants.Zero)  {
                    Rchain.set(i, j, 0.0)
                } else {
                    Rchain.set(i, j, Qchain.get(i, j) / Xchain.get(j) / Vchain!!.get(i, j))
                }
            }
        }

        for (i in infServers.indices) {
            val row: Int = infServers.get(i)!!
            for (j in 0..<Rchain.numCols) {
                Rchain.set(row, j, Lchain.get(row, j) / Vchain!!.get(row, j))
            }
        }

        val Tchain = Matrix(Vchain!!.numRows, Vchain.numCols)
        for (i in 0..<Tchain.numRows) {
            for (j in 0..<Tchain.numCols) {
                Tchain.set(i, j, Xchain.get(j) * Vchain.get(i, j))
            }
        }

        val ret1 = snDeaggregateChainResults(sn,
            Lchain,
            ST,
            STchain,
            Vchain,
            alpha!!,
            null,
            null,
            Rchain,
            Tchain,
            null,
            Xchain)
        Q = ret1.Q
        U = ret1.U
        R = ret1.R
        T = ret1.T
        X = ret1.X
        STeff = ST.copy()

        val NPFQNret = npfqn_nonexp_approx(if (options.config.highvar == null) "interp" else options.config.highvar,
            sn,
            ST0,
            V,
            SCV,
            T,
            U,
            gamma,
            nservers)
        ST = NPFQNret.ST
        gamma = NPFQNret.gamma
        eta = NPFQNret.eta

        for (i in 0..<M) {
            tmp_eta.set(i, FastMath.abs(1 - eta.get(i) / eta_1.get(i)))
        }
    }
    Q!!.absEq()
    R!!.absEq()
    X!!.absEq()
    U!!.absEq()
    X.removeInfinity()
    U.removeInfinity()
    Q.removeInfinity()
    R.removeInfinity()
    // Chain population corrections (from MATLAB commit 8103b59)
    for (c in 0..<C) {
        val inchain: Matrix = sn.inchain.get(c)!!
        val Nchain_c = Nchain.get(c)
        if (java.lang.Double.isFinite(Nchain_c)) {
            var sumQ = 0.0
            for (col in 0..<inchain.numCols) {
                val classIdx = inchain.get(0, col).toInt()
                sumQ += Q!!.sumCols(classIdx)
            }
            val ratio: Double
            if (sumQ > 0) {
                ratio = Nchain_c / sumQ
            } else {
                ratio = 0.0
            }
            if (sumQ > 0) {
                for (col in 0..<inchain.numCols) {
                    val classIdx = inchain.get(0, col).toInt()
                    for (i in 0..<M) {
                        Q.set(i, classIdx, ratio * Q.get(i, classIdx))
                    }
                    X!!.set(classIdx, ratio * X.get(classIdx))
                    for (i in 0..<M) {
                        T!!.set(i, classIdx, ratio * T.get(i, classIdx))
                    }
                }
            }
        }
    }

    // Handle self-looping classes - jobs stay at reference station
    for (k in 0..<K) {
        if (sn.isslc.get(k) == 1.0) {
            for (m in 0..<M) {
                Q!!.set(m, k, 0.0)
            }
            val ist = sn.refstat.get(k).toInt()
            Q!!.set(ist, k, sn.njobs.get(k))
            T!!.set(ist, k, sn.njobs.get(k) * sn.rates.get(ist, k))
            R!!.set(ist, k, Q.get(ist, k) / T.get(ist, k))
            U!!.set(ist, k, ST.get(ist, k) * T.get(ist, k))
        }
    }

    val runtime = (System.nanoTime() - startTime) / 1000000000.0

    return SolverNC.SolverNCReturn(Q, U, R, T, C, X, lG, STeff, it, runtime, method)
}