/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.mva.handlers

import jline.GlobalConstants
import jline.api.pfqn.mva.ab_amva
import jline.api.pfqn.mva.pfqn_aql
import jline.api.pfqn.mva.pfqn_bs
import jline.api.pfqn.mva.pfqn_conwayms
import jline.api.pfqn.mva.pfqn_linearizermx
import jline.api.pfqn.mva.pfqn_schmidt
import jline.api.pfqn.mva.pfqn_schmidt_ext
import jline.api.pfqn.mva.pfqn_sqni
import jline.lang.nodes.Queue
import jline.api.sn.*
import jline.io.line_error
import jline.io.mfilename
import jline.lang.NetworkStruct
import jline.lang.constant.NodeType
import jline.lang.constant.SchedStrategy
import jline.solvers.SolverOptions
import jline.solvers.mva.MVAResult
import jline.util.matrix.Matrix
import org.apache.commons.lang3.NotImplementedException
import java.lang.Double
import kotlin.Any
import kotlin.Array
import kotlin.Int
import kotlin.arrayOfNulls


/**
 * Handler for the solver_amva function.
 */
fun solver_amva(sn: NetworkStruct, options: SolverOptions): MVAResult {
    val chainReturn = snGetDemandsChain(sn)
    val Lchain = chainReturn.Dchain
    val STchain = chainReturn.STchain
    val Vchain = chainReturn.Vchain
    val alpha = chainReturn.alpha
    val Nchain = chainReturn.Nchain
    chainReturn.SCVchain
    chainReturn.refstatchain
    if (options.config.np_priority == null) {
        options.config.np_priority = "default"
    }
    if (options.config.multiserver == null) {
        options.config.multiserver = "default"
    }
    if (options.config.highvar == null) {
        options.config.highvar = "default"
    }
    when (options.method) {
        "amva.qli" -> options.method = "qli"
        "amva.qd", "amva.qdamva", "qdamva" -> options.method = "qd"
        "amva.lin" -> options.method = "lin"
        "amva.gflin" -> options.method = "gflin"
        "amva.egflin" -> options.method = "egflin"
        "amva.qdlin" -> options.method = "qdlin"
        "amva.fli" -> options.method = "fli"
        "amva.bs" -> options.method = "bs"
        "default", "amva" -> {
            var NchainSum = 0.0
            var flag = false
            var i = 0
            while (i < Nchain.numRows) {
                var j = 0
                while (j < Nchain.numCols) {
                    NchainSum += Nchain.get(i, j)
                    if (Nchain.get(i, j) < 1) flag = true
                    j++
                }
                i++
            }
            if (NchainSum <= 2 || flag) {
                options.method = "qd" // changing to bs degrades accuracy
            } else {
                options.method = "egflin"
                var i = 0
                while (i < sn.nstations) {
                    if (sn.nservers.get(i) > 1 && sn.nservers.get(i) < Int.Companion.MAX_VALUE) {
                        // if multi-server
                        options.method = "lin" // lin seems way worse than aql in test_LQN_8.xml
                    }
                    i++
                }
            }
        }
    }

    // trivial models
    if (snHasHomogeneousScheduling(sn, SchedStrategy.INF)) {
        options.config.multiserver = "default"
        return solver_amvald(sn, options)
    }

    val queueIdx = ArrayList<Int?>()
    val delayIdx = ArrayList<Int?>()
    val sourceIdx = ArrayList<Int?>()
    for (i in sn.nodetype.indices) {
        when (sn.nodetype.get(i)) {
            NodeType.Source -> sourceIdx.add(i)
            NodeType.Queue -> queueIdx.add(i)
            NodeType.Delay -> delayIdx.add(i)
            NodeType.Join -> queueIdx.add(i)
            else -> {}
        }
    }

    // Run amva method
    val M = sn.nstations
    val C = sn.nchains
    val V = Matrix(M, C)
    for (s in sourceIdx) {
        val i = sn.nodeToStation.get(s!!).toInt()
        for (c in 0..<sn.nchains) {
            var rateSum = 0.0
            for (sp in sourceIdx) {
                for (r in 0..<sn.chains.numCols) {
                    if (sn.chains.get(c, r) == 0.0) {
                        continue
                    }
                    if (!Double.isNaN(sn.rates.get(sn.nodeToStation.get(sp!!).toInt(), r))) {
                        rateSum += sn.rates.get(sn.nodeToStation.get(sp).toInt(), r)
                    }
                }
            }
            if (rateSum > 0) {
                V.set(i, c, 1)
            }
        }
    }
    var Q = Matrix(M, C)
    var U = Matrix(M, C)
    var X: Matrix? = null
    var totiter = 0
    if (snHasProductFormNotHetFCFS(sn) && !snHasLoadDependence(sn) && (!snHasOpenClasses(sn) || (snHasProductForm(sn) && snHasOpenClasses(sn) && options.method == "lin"))) {
        val ret = snGetProductFormChainParams(sn)
        // Note: ret returns snGetDemands, need to use correct properties
        val L0 = Matrix(ret.D)
        val L = ret.D
        val N = ret.N
        val Z = ret.Z
        val Z0 = Matrix(ret.Z)
        val nservers = ret.S
        val lambda = ret.lambda

        var vIdx = 0
        for (i in sn.nodetype.indices) {
            if (sn.nodetype.get(i) != NodeType.Queue && sn.nodetype.get(i) != NodeType.Delay) {
                continue
            }
            for (j in 0..<V.numCols) {
                V.set(sn.nodeToStation.get(i).toInt(), j, ret.V.get(vIdx, j))
            }
            vIdx++
        }
        when (options.config.multiserver) {
            "default", "seidmann" -> {
                val nserversRep = nservers.columnMajorOrder().repmat(1, C)
                // apply seidmann
                var i = 0
                while (i < L.numRows) {
                    var j = 0
                    while (j < L.numCols) {
                        L.set(i, j, L.get(i, j) / nserversRep.get(i, j))
                        j++
                    }
                    i++
                }
                var j = 0
                while (j < L.numRows) {
                    var k = 0
                    while (k < Z!!.numCols) {
                        Z.set(0, k, Z.get(0, k) + L0.get(j, k) * (nservers.get(j) - 1) / nservers.get(j))
                        k++
                    }
                    j++
                }
            }

            "softmin" -> return solver_amvald(sn, options)
        }
        when (options.method) {
            "sqni" -> {
                if (sn.nstations == 2) {
                    val result = pfqn_sqni(N, L, Z!!) // Kotlin function call

                    Q = result.Q.copy()
                    U = result.U.copy()
                    X = result.X.copy()

                    totiter = 1
                } else {
                    line_error(mfilename(object : Any() {}),
                        "SQNI cannot handle more than a single queue with a delay.")
                }
                val schd = arrayOfNulls<SchedStrategy>(queueIdx.size)
                var i = 0
                while (i < queueIdx.size) {
                    schd[i] = sn.sched.get(sn.stations.get(queueIdx.get(i)!!))
                    i++
                }
                val bsret = pfqn_bs(L, N, Z!!, options.tol, options.iter_max, null, schd)
                X = bsret.X
                var iResult = 0
                for (i in queueIdx) {
                    var j = 0
                    while (j < C) {
                        Q.set(sn.nodeToStation.get(i!!).toInt(), j, bsret.Q.get(iResult, j))
                        U.set(sn.nodeToStation.get(i).toInt(), j, bsret.U.get(iResult, j))
                        j++
                    }
                    iResult++
                }
                totiter = bsret.totiter
            }

            "bs" -> {
                val schd = arrayOfNulls<SchedStrategy>(queueIdx.size)
                var i = 0
                while (i < queueIdx.size) {
                    schd[i] = sn.sched.get(sn.stations.get(queueIdx.get(i)!!))
                    i++
                }
                val bsret = pfqn_bs(L, N, Z!!, options.tol, options.iter_max, null, schd)
                X = bsret.X
                var iResult = 0
                for (i in queueIdx) {
                    var j = 0
                    while (j < C) {
                        Q.set(sn.nodeToStation.get(i!!).toInt(), j, bsret.Q.get(iResult, j))
                        U.set(sn.nodeToStation.get(i).toInt(), j, bsret.U.get(iResult, j))
                        j++
                    }
                    iResult++
                }
                totiter = bsret.totiter
            }

            "aql" -> {
                if (snHasMultiServer(sn)) {
                    line_error(mfilename(object : Any() {}),
                        "AQL cannot handle multi-server stations. Try with the \"default\" or \"lin\" methods.")
                }
                val aqlret = pfqn_aql(L, N, Z, options.tol, options.iter_max)
                X = aqlret.X
                var idxResult = 0
                for (i in queueIdx) {
                    var j = 0
                    while (j < C) {
                        Q.set(sn.nodeToStation.get(i!!).toInt(), j, aqlret.Q.get(idxResult, j))
                        U.set(sn.nodeToStation.get(i).toInt(), j, aqlret.U.get(idxResult, j))
                        j++
                    }
                    idxResult++
                }
                totiter = aqlret.totiter
            }

            "ab" -> {
                val schedStrategies = sn.stations.map { station ->
                    if (station is Queue) station.schedStrategy else SchedStrategy.INF
                }
                val abres = ab_amva(STchain, N, V, sn.nservers, schedStrategies, false, "ab")
                X = abres.X
                var idxResult = 0
                for (i in queueIdx) {
                    var j = 0
                    while (j < C) {
                        Q.set(sn.nodeToStation.get(i!!).toInt(), j, abres.Q.get(idxResult, j))
                        U.set(sn.nodeToStation.get(i).toInt(), j, abres.U.get(idxResult, j))
                        j++
                    }
                    idxResult++
                }
            }

            "schmidt" -> {
                val schedStrategies = sn.stations.map { station ->
                    if (station is Queue) station.schedStrategy else SchedStrategy.INF
                }
                val schmidtres = pfqn_schmidt(sn.rates, sn.njobs, sn.nservers, V, schedStrategies)
                X = schmidtres.X
                var idxResult = 0
                for (i in queueIdx) {
                    var j = 0
                    while (j < C) {
                        Q.set(sn.nodeToStation.get(i!!).toInt(), j, schmidtres.Q.get(idxResult, j))
                        U.set(sn.nodeToStation.get(i).toInt(), j, schmidtres.U.get(idxResult, j))
                        j++
                    }
                    idxResult++
                }
            }

            "schmidt-ext" -> {
                val schedStrategies = sn.stations.map { station ->
                    if (station is Queue) station.schedStrategy else SchedStrategy.INF
                }
                val schmidtExtRes = pfqn_schmidt_ext(sn.rates, sn.njobs, sn.nservers, V, schedStrategies)
                X = schmidtExtRes.X
                var idxResult = 0
                for (i in queueIdx) {
                    var j = 0
                    while (j < C) {
                        Q.set(sn.nodeToStation.get(i!!).toInt(), j, schmidtExtRes.Q.get(idxResult, j))
                        U.set(sn.nodeToStation.get(i).toInt(), j, schmidtExtRes.U.get(idxResult, j))
                        j++
                    }
                    idxResult++
                }
            }

            "lin", "gflin", "egflin" -> if (nservers.elementMax() == 1.0) {
                val schdi = arrayOfNulls<SchedStrategy>(queueIdx.size + delayIdx.size)
                var idx = 0
                var i = 0
                while (i < sn.nodetype.size) {
                    if (sn.nodetype.get(i) != NodeType.Queue && sn.nodetype.get(i) != NodeType.Delay) {
                        i++
                        continue
                    }
                    schdi[idx] = sn.sched.get(sn.stations.get(sn.nodeToStation.get(i).toInt()))
                    idx++
                    i++
                }

                val res = pfqn_linearizermx(lambda,
                    L,
                    N,
                    Z!!,
                    nservers,
                    @Suppress("UNCHECKED_CAST") (schdi as Array<SchedStrategy>),
                    options.tol,
                    options.iter_max,
                    options.method)

                var iRes = 0
                for (i in queueIdx) {
                    var j = 0
                    while (j < C) {
                        Q.set(sn.nodeToStation.get(i!!).toInt(), j, res.Q.get(iRes, j))
                        U.set(sn.nodeToStation.get(i).toInt(), j, res.U.get(iRes, j))
                        j++
                    }
                    iRes++
                }
                X = res.X
                totiter = res.totiter
            } else {
                when (options.config.multiserver) {
                    "conway" -> {
                        var schdi = arrayOfNulls<SchedStrategy>(queueIdx.size + delayIdx.size)
                        var idx = 0
                        var i = 0
                        while (i < sn.nodetype.size) {
                            if (sn.nodetype.get(i) != NodeType.Queue && sn.nodetype.get(i) != NodeType.Delay) {
                                i++
                                continue
                            }
                            schdi[idx] = sn.sched.get(sn.stations.get(sn.nodeToStation.get(i).toInt()))
                            idx++
                            i++
                        }
                        val res = pfqn_conwayms(L,
                            N,
                            Z!!,
                            nservers.toIntArray1D(),
                            schdi as Array<SchedStrategy?>,
                            options.tol,
                            options.iter_max)
                        var iRes = 0
                        for (i in queueIdx) {
                            var j = 0
                            while (j < C) {
                                Q.set(sn.nodeToStation.get(i!!).toInt(), j, res.Q.get(iRes, j))
                                U.set(sn.nodeToStation.get(i).toInt(), j, res.U.get(iRes, j))
                                j++
                            }
                            iRes++
                        }
                        X = res.X
                        totiter = res.totiter
                    }

                    "erlang" -> {
                        // Erlang multiserver method not yet implemented, falling back to default
                        options.config.multiserver = "default"
                        return solver_amvald(sn, options)
                    }

                    "krzesinski" -> {
                        var schdi = arrayOfNulls<SchedStrategy>(queueIdx.size + delayIdx.size)
                        var idx = 0
                        var i = 0
                        while (i < sn.nodetype.size) {
                            if (sn.nodetype.get(i) != NodeType.Queue && sn.nodetype.get(i) != NodeType.Delay) {
                                i++
                                continue
                            }
                            schdi[idx] = sn.sched.get(sn.stations.get(sn.nodeToStation.get(i).toInt()))
                            idx++
                            i++
                        }
                        val res = pfqn_linearizermx(lambda,
                            L,
                            N,
                            Z!!,
                            nservers,
                            @Suppress("UNCHECKED_CAST") (schdi as Array<SchedStrategy>),
                            options.tol,
                            options.iter_max,
                            "default")
                        var iRes = 0
                        for (i in queueIdx) {
                            var j = 0
                            while (j < C) {
                                Q.set(sn.nodeToStation.get(i!!).toInt(), j, res.Q.get(iRes, j))
                                U.set(sn.nodeToStation.get(i).toInt(), j, res.U.get(iRes, j))
                                j++
                            }
                            iRes++
                        }
                        X = res.X
                        totiter = res.totiter
                    }

                    "default", "softmin", "seidmann" -> return solver_amvald(sn, options)
                }
            }

            else -> {
                when (options.config.multiserver) {
                    "conway", "erlang", "krzesinski" -> options.config.multiserver = "default"
                }
                return solver_amvald(sn, options)
            }
        }

        // Compute performance at delay, then unapply seidmann if needed
        for (i in 0..<Z0.numRows) {
            if (delayIdx.isNotEmpty()) {
                val mult = X!!.repmat(delayIdx.size, 1).elementMult(Z, null)
                var zidx = 0
                for (d in delayIdx) {
                    for (j in 0..<Q.numCols) {
                        Q.set(sn.nodeToStation.get(d!!).toInt(), j, mult.get(zidx, j))
                        U.set(sn.nodeToStation.get(d).toInt(), j, mult.get(zidx, j))
                    }
                    zidx++
                }
            }
            when (options.config.multiserver) {
                "default", "seidmann" -> {
                    // Skip Seidmann un-apply for ab and schmidt methods as it removes queue length
                    if (options.method != "ab" && !options.method.startsWith("schmidt")) {
                        var j = 0
                        while (j < L.numRows) {
                            if (i == 0 && nservers.get(j) > 1) {
                                // un-apply seidmann from first delay and move it to the origin queue
                                var k = 0
                                while (k <= j) {
                                    if (k >= queueIdx.size) break
                                    val jq: Int = queueIdx.get(k)!!
                                    var l = 0
                                    while (l < Q.numCols) {
                                        Q.set(jq,
                                            l,
                                            Q.get(jq, l) + (L0.get(j, l) * (nservers.get(j) - 1) / nservers.get(j)) * X!!.get(
                                                l))
                                        l++
                                    }
                                    k++
                                }
                            }
                            j++
                        }
                    }
                }
            }
        }
        val T = V.elementMult(X!!.repmat(M, 1), null)
        val R = Matrix(Q.numRows, Q.numCols)
        for (i in 0..<R.numRows) {
            for (j in 0..<R.numCols) {
                if (V.get(i,j) <= GlobalConstants.Zero) {
                    R.set(i,j,0.0);
                } else {
                    R.set(i, j, Q.get(i, j) / T.get(i, j))
                }
            }
        }
        val Cm = Matrix(N.numRows, N.numCols)
        for (i in 0..<Cm.numRows) {
            for (j in 0..<Cm.numCols) {
                Cm.set(i, j, N.get(i, j) / X.get(i, j) - Z!!.get(i, j))
            }
        }
        val lG = kotlin.Double.Companion.NaN
        val result = MVAResult()
        result.QN = Q
        result.UN = U
        result.RN = R
        result.TN = T
        result.CN = Cm
        result.XN = X
        result.logNormConstAggr = lG
        result.iter = totiter
        result.method = options.method
        if (snHasClassSwitching(sn)) {
            val ret1 = snDeaggregateChainResults(sn, Lchain, null, STchain, Vchain, alpha, null, null, R, T, null, X)
            result.QN = ret1.Q
            result.UN = ret1.U
            result.RN = ret1.R
            result.TN = ret1.T
            result.CN = ret1.C
            result.XN = ret1.X
        }
        return result
    } else {
        when (options.config.multiserver) {
            "conway", "erlang", "krzesinski" -> options.config.multiserver = "default"
        }
        return solver_amvald(sn, options)
    }
}
