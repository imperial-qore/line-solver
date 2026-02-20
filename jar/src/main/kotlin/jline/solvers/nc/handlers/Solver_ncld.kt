package jline.solvers.nc.handlers

import jline.api.npfqn.npfqn_nonexp_approx
import jline.api.pfqn.ld.pfqn_fnc
import jline.api.pfqn.ld.pfqn_mushift
import jline.api.pfqn.ld.pfqn_ncld
import jline.api.sn.snDeaggregateChainResults
import jline.api.sn.snGetDemandsChain
import jline.api.sn.snGetProductFormParams
import jline.api.sn.snHasJointDependence
import jline.io.Ret
import jline.io.line_warning
import jline.lang.NetworkStruct
import jline.lang.constant.SchedStrategy
import jline.solvers.SolverOptions
import jline.solvers.nc.SolverNC
import jline.util.Utils
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath
import java.util.*
import java.util.function.DoublePredicate
import kotlin.math.ceil
import kotlin.math.exp

fun solver_ncld(sn: NetworkStruct, options: SolverOptions): SolverNC.SolverNCLDReturn {
    var M = sn.nstations
    var K = sn.nclasses
    val nservers = sn.nservers
    val nserversFinite = nservers.copy()
    nserversFinite.removeInfinity()
    var minFiniteServer = Double.Companion.MAX_VALUE
    for (i in 0..<nservers.numElements) {
        if (java.lang.Double.isFinite(nservers.get(i)) && nservers.get(i) < minFiniteServer) {
            minFiniteServer = nservers.get(i)
        }
    }
    // Only check for multi-server if finite servers were found (not all delay nodes)
    if (minFiniteServer < Double.MAX_VALUE && minFiniteServer > 1) {
        if (sn.lldscaling.isEmpty && M == 2 && java.lang.Double.isFinite(sn.njobs.elementMaxAbs())) {
            val Nt = sn.njobs.elementSum()
            sn.lldscaling = sn.lldscaling.concatCols(Matrix(M, Nt.toInt()))
            for (i in 0..<M) {
                var j = 0
                while (j < Nt) {
                    sn.lldscaling.set(i, j, FastMath.min((j + 1).toDouble(), sn.nservers.get(i)))
                    j++
                }
            }
        } else {
            throw RuntimeException("The load-dependent solver does not support multi-server stations yet. Specify multi-server stations via limited load-dependence.")
        }
    }
    if (sn.cdscaling != null && !sn.cdscaling.isEmpty() && options.method.equals("exact", ignoreCase = true)) {
        throw RuntimeException("Exact class-dependent solver not yet available in NC.")
    }
    if (snHasJointDependence(sn) && options.method.equals("exact", ignoreCase = true)) {
        throw RuntimeException("Exact joint-dependent solver not yet available in NC.")
    }
    val NK = sn.njobs.transpose()
    if (Utils.isInf(NK.elementMax())) {
        throw RuntimeException("The load-dependent solver does not support open classes yet.")
    }
    var C = sn.nchains
    val SCV = sn.scv
    var gamma = Matrix(M, 1)
    gamma.zero()
    var V = Matrix(sn.nstateful, K)
    for (i in 0..<sn.visits.size) {
        V = V.add(1.0, sn.visits.get(i))
    }
    var ST = Matrix.ones(sn.rates.numRows, sn.rates.numCols).elementDiv(sn.rates)
    for (i in 0..<ST.numRows) {
        for (j in 0..<ST.numCols) {
            if (java.lang.Double.isNaN(ST.get(i, j))) {
                ST.set(i, j, 0)
            }
        }
    }
    val ST0 = ST.copy()
    var lldscaling = sn.lldscaling
    val NKfinite = NK.copy()
    NKfinite.removeInfinity()
    var Nt = NKfinite.elementSum()
    if (lldscaling.isEmpty) {
        lldscaling = Matrix.ones(M, FastMath.ceil(Nt).toInt())
    }
    val demandsChainReturn = snGetDemandsChain(sn)
    val Vchain = demandsChainReturn.Vchain
    val alpha = demandsChainReturn.alpha
    var eta_1 = Matrix(1, M)
    eta_1.zero()
    var eta = Matrix.ones(1, M)
    if (!sn.sched.containsValue(SchedStrategy.FCFS)) {
        options.iter_max = 1
    }
    var iter = 0
    var Tstart = System.nanoTime()
    var snDeaggragatedChains: Ret.snDeaggregateChainResults? = null
    val lambda: Matrix? = null
    var Lchain: Matrix? = null
    var STchain: Matrix? = null
    var Q: Matrix? = null
    var Cmat: Matrix? = null
    var U: Matrix? = null
    var R: Matrix? = null
    var T: Matrix? = null
    var X: Matrix? = null
    var lG: Double? = null
    var method: String? = null
    while (Matrix.ones(eta.numRows, eta.numCols).sub(eta.elementDiv(eta_1))
            .elementMaxAbs() > options.iter_tol && iter < options.iter_max) {
        iter += 1
        eta_1 = eta.copy()
        M = sn.nstations
        K = sn.nclasses
        C = sn.nchains
        Lchain = Matrix(M, C)
        Lchain.zero()
        STchain = Matrix(M, C)
        STchain.zero()
        val SCVchain = Matrix(M, C)
        SCVchain.zero()
        val Nchain = Matrix(1, C)
        Nchain.zero()
        val refstatchain = Matrix(C, 1)
        refstatchain.zero()
        for (c in 0..<C) {
            val inchain: Matrix = sn.inchain.get(c)!!
            var isOpenChain = false
            for (i in 0..<inchain.numElements) {
                if (Utils.isInf(sn.njobs.get(inchain.get(i).toInt()))) {
                    isOpenChain = true
                }
            }
            for (i in 0..<M) {
                val STinchain = Matrix(1, inchain.numElements)
                val alphainchain = Matrix(1, inchain.numElements)
                val SCVinchain = Matrix(1, inchain.numElements)
                for (j in 0..<inchain.numElements) {
                    STinchain.set(j, ST.get(i, inchain.get(j).toInt()))
                    alphainchain.set(j, alpha.get(i, inchain.get(j).toInt()))
                    SCVinchain.set(j, SCV.get(i, inchain.get(j).toInt()))
                }
                val lchainVal = Vchain.get(i, c) * STinchain.mult(alphainchain.transpose()).toDouble()
                Lchain.set(i, c, lchainVal)
                STchain.set(i, c, STinchain.mult(alphainchain.transpose()).toDouble())
                if (isOpenChain && i.toDouble() == sn.refstat.get(inchain.get(0).toInt())) {
                    val STinchainFinite = STinchain.copy()
                    STinchainFinite.removeInfinity()
                    STchain.set(i, c, STinchainFinite.elementSum())
                } else {
                    STchain.set(i, c, STinchain.mult(alphainchain.transpose()).toDouble())
                }
                SCVchain.set(i, c, SCVinchain.mult(alphainchain.transpose()).toDouble())
            }
            val NKinchain = Matrix(1, inchain.numElements)
            for (i in 0..<inchain.numElements) {
                NKinchain.set(i, NK.get(inchain.get(i).toInt()))
            }
            Nchain.set(c, NKinchain.elementSum())
            refstatchain.set(c, sn.refstat.get(inchain.get(0).toInt()))
            if ((sn.refstat.get(inchain.get(0).toInt()) - refstatchain.get(c)) != 0.0) {
                throw RuntimeException(String.format("Classes in chain %d have different reference station.", c))
            }
        }
        for (i in 0..<STchain.numRows) {
            for (j in 0..<STchain.numCols) {
                if (!java.lang.Double.isFinite(STchain.get(i, j))) {
                    STchain.set(i, j, 0)
                }
                if (!java.lang.Double.isFinite(Lchain.get(i, j))) {
                    Lchain.set(i, j, 0)
                }
            }
        }
        Tstart = System.nanoTime()
        val Nchainfinite = Nchain.copy()
        Nchainfinite.removeInfinity()
        Nt = Nchainfinite.elementSum()
        val L = Matrix(M, C)
        L.zero()
        val mu = Matrix(M, FastMath.ceil(Nt).toInt())
        val infServers: MutableList<Int?> = ArrayList<Int?>()
        var Z = L.copy()
         for (i in 0..<M) {
            if (Utils.isInf(nservers.get(i))) {
                infServers.add(i)
                for (j in 0..<C) {
                    L.set(i, j, Lchain.get(i, j))
                    Z.set(i, j, Lchain.get(i, j))
                }
                var j = 0
                while (j < ceil(Nt)) {
                    mu.set(i, j, j + 1)
                    j++
                }
            } else {
                if (options.method.equals("exact", ignoreCase = true) && nservers.get(i) > 1) {
                    line_warning("solver_ncld", "SolverNC does not support exact multiserver yet. Switching to approximate method.")
                }
                for (j in 0..<C) {
                    L.set(i, j, Lchain.get(i, j))
                }
                var j = 0
                while (j < ceil(Nt)) {
                    mu.set(i, j, lldscaling.get(i, j))
                    j++
                }
            }
        }
        val Qchain = Matrix(M, C)
        Qchain.zero()
        val Nchain0 = Nchain.copy()
        Nchain0.zero()
        val ret = pfqn_ncld(L, Nchain, Nchain0, mu, options)
        lG = ret.lG
        method = ret.method
        var Xchain = Matrix(1, 0)
        if (Xchain.isEmpty) {
            val lGr = Matrix(1, C)
            val lGhat_fnci = Matrix(1, C)
            val lGhatir = Matrix(1, C)
            val lGr_i = Matrix(1, C)
            val ldDemand = Matrix(M, C)
            for (r in 0..<C) {
                val Nchain_r = Matrix.oner(Nchain, mutableListOf<Int?>(r))
                lGr.set(r, pfqn_ncld(L, Nchain_r, Nchain0, mu, options).lG)
                val xchainVal = exp(lGr.get(r) - lG)
                Xchain = Xchain.concatCols(Matrix.singleton(xchainVal))
                for (i in 0..<M) {
                    Qchain.set(i, r, 0)
                }
                val CQchain_r = Matrix(M, 1)
                CQchain_r.zero()
                if (M == 2 && Utils.isInf(sn.nservers.elementMaxAbs())) {
                    var firstDelay = -1
                    for (i in 0..<sn.nservers.numElements) {
                        if (Utils.isInf(sn.nservers.get(i))) {
                            firstDelay = i
                            break
                        }
                    }
                    Qchain.set(firstDelay, r, Lchain.get(firstDelay, r) * Xchain.get(r))
                    for (i in 0..<M) {
                        if (i != firstDelay) {
                            Qchain.set(i, r, Nchain.get(r) - Lchain.get(firstDelay, r) * Xchain.get(r))
                        }
                    }
                } else {
                    for (i in 0..<M) {
                        val Lms_i = Matrix.extractRows(L, i, i + 1, null)
                        Matrix.extractRows(mu, i, i + 1, null)
                        val muhati = pfqn_mushift(mu, i)
                        val fncret = pfqn_fnc(Matrix.extractRows(muhati, i, i + 1, null))
                        val muhati_f = fncret.mu
                        val c = fncret.c.toDouble()
                        if (Lchain.get(i, r) > 0) {
                            if (Utils.isInf(nservers.get(i))) {
                                Qchain.set(i, r, Lchain.get(i, r) * Xchain.get(r))
                            } else {
                                if (i == (M - 1) && nserversFinite.elementSum() == 1.0) {
                                    var Lchainsum = 0.0
                                    var Qchainsum = 0.0
                                    for (j in 0..<nservers.numElements) {
                                        if (Utils.isInf(nservers.get(j))) {
                                            Lchainsum += Lchain.get(j, r)
                                        } else {
                                            Qchainsum += Qchain.get(j, r)
                                        }
                                    }
                                    Qchain.set(i,
                                        r,
                                        FastMath.max(0.0, Nchain.get(r) - Lchainsum * Xchain.get(r) - Qchainsum))
                                } else {
                                    lGhat_fnci.set(r,
                                        pfqn_ncld(Matrix.concatRows(L, Matrix.extractRows(L, i, i + 1, null), null),
                                            Nchain_r,
                                            Nchain0,
                                            Matrix.concatRows(muhati, muhati_f, null),
                                            options).lG)
                                    lGhatir.set(r, pfqn_ncld(L, Nchain_r, Nchain0, muhati, options).lG)
                                    lGr_i.set(r, pfqn_ncld(Lms_i, Nchain_r, Nchain0, muhati, options).lG)
                                    val dlGa = lGhat_fnci.get(r) - lGhatir.get(r)
                                    val dlG_i = lGr_i.get(r) - lGhatir.get(r)
                                    CQchain_r.set(i, (exp(dlGa) - 1) + c * (exp(dlG_i) - 1))
                                    ldDemand.set(i,
                                        r,
                                        FastMath.log(L.get(i, r)) + lGhatir.get(r) - FastMath.log(mu.get(i,
                                            0)) - lGr.get(r))
                                    Qchain.set(i,
                                        r,
                                        FastMath.exp(ldDemand.get(i, r)) * Xchain.get(r) * (1 + CQchain_r.get(i)))
                                }
                            }
                        }
                    }
                }
            }
        } else {
            for (r in 0..<C) {
                for (i in 1..<M) {
                    if (Lchain.get(i, r) > 0) {
                        if (Utils.isInf(nservers.get(i))) {
                            Qchain.set(i, r, Lchain.get(i, r) * Xchain.get(r))
                        }
                    }
                }
            }
        }
        if (Arrays.stream(Xchain.toArray1D()).allMatch(DoublePredicate { v: Double -> java.lang.Double.isNaN(v) })) {
            line_warning("solver_ncld", "Normalizing constant computations produced a floating-point range exception. Model is likely too large.")
        }
        Z = Z.sumCols()
        val Rchain = Qchain.elementDivide(Xchain.repmat(M, 1)).elementDivide(Vchain)
        for (i in infServers) {
            for (j in 0..<Rchain.numCols) {
                Rchain.set(i!!, j, Lchain.get(i, j) / Vchain.get(i, j))
            }
        }
        val Tchain = Xchain.repmat(M, 1).elementMult(Vchain, null)
        Tchain.elementMult(Lchain, null)
        Nchain.elementDiv(Xchain).sub(Z)
        snDeaggragatedChains =
            snDeaggregateChainResults(sn, Lchain, ST, STchain, Vchain, alpha, null, null, Rchain, Tchain, null, Xchain)
        Q = snDeaggragatedChains.Q
        U = snDeaggragatedChains.U
        R = snDeaggragatedChains.R
        T = snDeaggragatedChains.T
        Cmat = snDeaggragatedChains.C
        X = snDeaggragatedChains.X
        val NPFQNret = npfqn_nonexp_approx(if (options.config.highvar == null) "default" else options.config.highvar,
            sn,
            ST0.copy(),
            V,
            SCV,
            T,
            U,
            gamma,
            nservers)
        ST = NPFQNret.ST
        gamma = NPFQNret.gamma
        eta = NPFQNret.eta.transpose()
    }
    snGetProductFormParams(sn)
    val runtime = (System.nanoTime() - Tstart) / 1000000000.0
    Q = snDeaggragatedChains!!.Q
    Q.absEq()
    R = snDeaggragatedChains.R
    R.absEq()
    U = snDeaggragatedChains.U
    U.absEq()
    X = snDeaggragatedChains.X
    X.absEq()
    for (i in 0..<M) {
        val openClasses: MutableList<Int?> = ArrayList<Int?>()
        for (j in 0..<NK.numElements) {
            if (Utils.isInf(NK.get(j))) {
                openClasses.add(j)
            }
        }
        if (sn.nservers.get(i) > 1 && !Utils.isInf(sn.nservers.get(i))) {
            for (r in 0..<K) {
                val c = Matrix.extractColumn(sn.chains, r, null).find().value().toInt()
                if ((!openClasses.contains(r) && snDeaggragatedChains.X.get(r) > 0)) {
                    U.set(i,
                        r,
                        X.get(r) * sn.visits.get(c)!!.get(i, r) / sn.visits.get(c)!!
                            .get(sn.refstat.get(r).toInt(), r) * ST.get(i, r) / sn.nservers.get(i))
                } else if (openClasses.contains(r) && lambda!!.get(r) > 0) {
                    U.set(i,
                        r,
                        lambda.get(r) * sn.visits.get(c)!!.get(i, r) / sn.visits.get(c)!!
                            .get(sn.refstat.get(r).toInt(), r) * ST.get(i, r) / sn.nservers.get(i))
                }
            }
        } else if (Utils.isInf(sn.nservers.get(i))) {
            for (r in 0..<K) {
                val c = Matrix.extractColumn(sn.chains, r, null).find().get(0).toInt()
                if ((!openClasses.contains(r) && snDeaggragatedChains.X.get(r) > 0)) {
                    U.set(i,
                        r,
                        X.get(r) * sn.visits.get(c)!!.get(i, r) / sn.visits.get(c)!!
                            .get(sn.refstat.get(r).toInt(), r) * ST.get(i, r))
                } else if (openClasses.contains(r) && lambda!!.get(r) > 0) {
                    U.set(i,
                        r,
                        lambda.get(r) * sn.visits.get(c)!!.get(i, r) / sn.visits.get(c)!!
                            .get(sn.refstat.get(r).toInt(), r) * ST.get(i, r))
                }
            }
        } else {
            for (j in 0..<U.numCols) {
                U.set(i, j, U.get(i, j) / Matrix.extractRows(lldscaling, i, i + 1, null).elementMax())
            }
            if (Matrix.extractRows(U, i, i + 1, null).elementSum() > 1) {
                val Uinotnan = Matrix.extractRows(U, i, i + 1, null)
                Uinotnan.removeNaN()
                for (j in 0..<U.numCols) {
                    U.set(i, j, U.get(i, j) / Uinotnan.elementSum())
                }
            }
        } // Continue from line 213
    }
    X = snDeaggragatedChains.X
    X.apply(Double.Companion.POSITIVE_INFINITY, 0.0, "equal")
    X.apply(Double.Companion.NaN, 0.0, "equal")
    U.apply(Double.Companion.POSITIVE_INFINITY, 0.0, "equal")
    U.apply(Double.Companion.NaN, 0.0, "equal")
    Q.apply(Double.Companion.POSITIVE_INFINITY, 0.0, "equal")
    Q.apply(Double.Companion.NaN, 0.0, "equal")
    R.apply(Double.Companion.POSITIVE_INFINITY, 0.0, "equal")
    R.apply(Double.Companion.NaN, 0.0, "equal")

    // Chain population corrections
    val Nchain_final = Matrix(1, C)
    for (c in 0..<C) {
        val inchain: Matrix = sn.inchain.get(c)!!
        var sum = 0.0
        for (col in 0..<inchain.numCols) {
            sum += NK.get(inchain.get(0, col).toInt())
        }
        Nchain_final.set(c, sum)
    }

    for (c in 0..<C) {
        val inchain: Matrix = sn.inchain.get(c)!!
        val Nchain_c = Nchain_final.get(c)
        if (java.lang.Double.isFinite(Nchain_c)) {
            var sumQ = 0.0
            for (col in 0..<inchain.numCols) {
                val classIdx = inchain.get(0, col).toInt()
                sumQ += Q.sumCols(classIdx)
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
                    X.set(classIdx, ratio * X.get(classIdx))
                    for (i in 0..<M) {
                        snDeaggragatedChains.T.set(i, classIdx, ratio * snDeaggragatedChains.T.get(i, classIdx))
                    }
                    for (i in 0..<M) {
                        U.set(i, classIdx, ratio * U.get(i, classIdx))
                    }
                    for (i in 0..<M) {
                        if (snDeaggragatedChains.T.get(i, classIdx) > 0) {
                            R.set(i, classIdx, Q.get(i, classIdx) / snDeaggragatedChains.T.get(i, classIdx))
                        }
                    }
                }
            }
        }
    }

    // Handle self-looping classes - jobs stay at reference station
    for (k in 0..<K) {
        if (sn.isslc.get(k) == 1.0) {
            for (m in 0..<M) {
                Q.set(m, k, 0.0)
            }
            val ist = sn.refstat.get(k).toInt()
            Q.set(ist, k, sn.njobs.get(k))
            snDeaggragatedChains!!.T.set(ist, k, sn.njobs.get(k) * sn.rates.get(ist, k))
            R.set(ist, k, Q.get(ist, k) / snDeaggragatedChains.T.get(ist, k))
            U.set(ist, k, ST.get(ist, k) * snDeaggragatedChains.T.get(ist, k))
        }
    }

    return SolverNC.SolverNCLDReturn(Q, U, R, snDeaggragatedChains!!.T, Cmat, X, lG!!, runtime, iter, method)
}