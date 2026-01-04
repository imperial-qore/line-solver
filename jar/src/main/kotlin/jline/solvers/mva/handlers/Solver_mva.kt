/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.mva.handlers

import jline.GlobalConstants
import jline.GlobalConstants.Inf
import jline.api.pfqn.mva.pfqn_mvams
import jline.api.sn.snDeaggregateChainResults
import jline.api.sn.snGetDemandsChain
import jline.api.sn.snHasProductForm
import jline.lang.NetworkStruct
import jline.lang.constant.NodeType
import jline.lang.constant.SchedStrategy
import jline.solvers.SolverOptions
import jline.solvers.mva.MVAResult
import jline.util.Maths
import jline.util.Utils
import jline.util.matrix.Matrix
import java.lang.Double
import kotlin.Int
import kotlin.RuntimeException


/**
 * Handler for the solver_mva function.
 */
fun solver_mva(sn: NetworkStruct, options: SolverOptions): MVAResult {
    val ret = snGetDemandsChain(sn)
    val Lchain = ret.Dchain
    val STchain = ret.STchain
    val Vchain = ret.Vchain
    val alpha = ret.alpha
    val Nchain = ret.Nchain
    val refstatchain = ret.refstatchain
    val nservers = sn.nservers
    // We're using schedStrategies instead of schedid, as schedid is not available in NetworkStruct
    val M = sn.nstations
    val K = sn.nchains

    // Check for LCFS scheduling - not supported in this version
    for (i in 0..<M) {
        when (sn.sched.get(sn.stations.get(i))) {
            SchedStrategy.LCFS, SchedStrategy.LCFSPR ->
                throw RuntimeException("LCFS queueing networks are not supported in this version.")
            else -> {}
        }
    }

    // For non-LCFS models, check product-form requirement
    val infSET = ArrayList<Int?>() // set of infinite server stations
    val qSET = ArrayList<Int?>() // set of other product-form stations
    if (!snHasProductForm(sn)){
        throw RuntimeException("Unsupported exact MVA analysis, the model does not have a product form")
    }

    for (i in 0..<M) {
        when (sn.sched.get(sn.stations.get(i))) {
            SchedStrategy.EXT -> {}
            SchedStrategy.INF -> infSET.add(i)
            SchedStrategy.PS, SchedStrategy.LCFSPR -> qSET.add(i)
            SchedStrategy.FCFS, SchedStrategy.SIRO -> qSET.add(i)
            else -> throw RuntimeException("Unsupported exact MVA analysis for " + SchedStrategy.toText(sn.sched.get(sn.stations.get(
                i))) + " scheduling")
        }
    }
    val Uchain = Matrix(M, K)
    val Tchain = Matrix(M, K)
    val C = Matrix(1, K)
    val Wchain = Matrix(M, K)
    val Qchain = Matrix(M, K)
    val lambda = Matrix(1, K)

    val ocl = ArrayList<Int?>()
    for (i in 0..<Nchain.numCols) {
        if (Utils.isInf(Nchain.get(i))) {
            ocl.add(i)
        }
    }
    for (r in ocl) { // open classes
        lambda.set(0, r!!, 1.0 / STchain.get(refstatchain.get(r).toInt(), r))
        Qchain.set(refstatchain.get(r).toInt(), r, Inf)
    }
    val rset = ArrayList<Int?>()
    for (i in 0..<K) {
        if (Nchain.get(i) != 0.0) {
            rset.add(i)
        }
    }
    val Lp = Matrix(qSET.size, STchain.numCols)
    for (i in qSET.indices) {
        for (j in 0..<Lp.numCols) {
            Lp.set(i, j, STchain.get(qSET.get(i)!!, j) * Vchain.get(qSET.get(i)!!, j))
        }
    }
    val Zp = Matrix(infSET.size, STchain.numCols)
    for (i in infSET.indices) {
        for (j in 0..<Zp.numCols) {
            Zp.set(i, j, STchain.get(infSET.get(i)!!, j) * Vchain.get(infSET.get(i)!!, j))
        }
    }
    val nserversp = Matrix(qSET.size, 1)
    for (i in qSET.indices) {
        nserversp.set(i, nservers.get(qSET.get(i)!!))
    }
    val ret1 = pfqn_mvams(lambda, Lp, Nchain, Zp, Matrix.ones(qSET.size, 1), nserversp)
    val Xchain = ret1.X
    val Qpf = ret1.Q/* Uchain = ret1.UN; -> not needed */
    val lG = ret1.lGN
    for (i in qSET.indices) {
        for (j in 0..<Qchain.numCols) {
            Qchain.set(qSET.get(i)!!, j, Qpf.get(i, j))
        }
    }

    val Q2 = Matrix(infSET.size, Qchain.numCols)
    val Xchainrep = Xchain.repmat(infSET.size, 1)
    for (i in infSET.indices) {
        for (j in 0..<Q2.numCols) {
            Q2.set(i, j, Xchainrep.get(i, j) * STchain.get(infSET.get(i)!!, j) * Vchain.get(infSET.get(i)!!, j))
        }
    }

    for (i in infSET.indices) {
        for (j in 0..<Qchain.numCols) {
            Qchain.set(infSET.get(i)!!, j, Q2.get(i, j))
        }
    }

    val ccl = ArrayList<Int?>()
    for (i in 0..<Nchain.numCols) {
        if (Double.isFinite(Nchain.get(i))) {
            ccl.add(i)
        }
    }
    for (r in rset) {
        for (k in infSET) {
            Wchain.set(k!!, r!!, STchain.get(k, r))
        }
        for (k in qSET) {
            if (Utils.isInf(nservers.get(k!!))) { // Infinite server
                Wchain.set(k, r!!, STchain.get(k, r))
            } else {
                if (Vchain.get(k, r!!) == 0.0 || Xchain.get(r) == 0.0) {
                    Wchain.set(k, r, 0)
                } else {
                    Wchain.set(k, r, Qchain.get(k, r) / (Xchain.get(r) * Vchain.get(k, r)))
                }
            }
        }
    }
    for (r in rset) {
        if (Matrix.extractColumn(Wchain, r!!, null).elementSum() == 0.0) {
            Xchain.set(r, 0.0)
        } else {
            if (Utils.isInf(Nchain.get(r))) {
                val vt = Matrix.extractColumn(Vchain, r, null)
                val wt = Matrix.extractColumn(Wchain, r, null)
                C.set(r, vt.transpose().mult(wt).get(0))
            } else if (Nchain.get(r) == 0.0) {
                Xchain.set(r, 0.0)
                C.set(r, 0.0)
            } else {
                val vt = Matrix.extractColumn(Vchain, r, null)
                val wt = Matrix.extractColumn(Wchain, r, null)
                C.set(r, vt.transpose().mult(wt).get(0))
                Xchain.set(r, Nchain.get(r) / C.get(r))
            }
        }

        for (k in 0..<M) {
            Qchain.set(k, r, Xchain.get(r) * Vchain.get(k, r) * Wchain.get(k, r))
            Tchain.set(k, r, Xchain.get(r) * Vchain.get(k, r))
        }
    }
    for (k in 0..<M) {
        for (r in rset) {
            if (Utils.isInf(nservers.get(k))) { // infinite server
                Uchain.set(k, r!!, Vchain.get(k, r) * STchain.get(k, r) * Xchain.get(r))
            } else {
                Uchain.set(k, r!!, Vchain.get(k, r) * STchain.get(k, r) * Xchain.get(r) / nservers.get(k))
            }
        }
    }
    for (k in 0..<M) {
        for (r in 0..<K) {
            if (Vchain.get(k, r) * STchain.get(k, r) > options.tol) {
                when (sn.sched.get(sn.stations.get(k))) {
                    SchedStrategy.FCFS, SchedStrategy.PS -> {
                        val Urow = Matrix.extractRows(Uchain, k, k + 1, null)
                        val UrowSum = Urow.elementSum()
                        if (UrowSum > 1 + options.tol) {
                            val Vrow = Matrix.extractRows(Vchain, k, k + 1, null)
                            val STrow = Matrix.extractRows(STchain, k, k + 1, null)
                            Uchain.set(k,
                                r,
                                Maths.min(1.0, UrowSum) * Vchain.get(k, r) * STchain.get(k,
                                    r) * Xchain.get(r) / Vrow.elementMult(STrow, null).mult(Xchain.columnMajorOrder())
                                    .get(0))
                        }
                    }
                    else -> {}
                }
            }
        }
    }
    var Vs: Matrix? = null
    for (key in sn.nodevisits.keys) {
        if (Vs == null) {
            Vs = sn.nodevisits.get(key)
        } else {
            Vs = Vs.add(1.0, sn.nodevisits.get(key))
        }
    }
    val sinks = ArrayList<Int?>()
    for (i in sn.nodetype.indices) {
        if (sn.nodetype.get(i) == NodeType.Sink) {
            sinks.add(i)
        }
    }
//    val Vsink = Matrix(sinks.size, Vs!!.getNumCols())
//    for (i in sinks.indices) {
//        for (j in 0..<Vsink.getNumCols()) {
//            Vsink.set(i, j, Vs.get(sinks.get(i)!!, j))
//        }
//    }
//    for (r in 0..<Nchain.getNumCols()) {
//        if (Utils.isInf(Nchain.get(r))) { // open classes
//            Xchain.set(r, Vsink.get(r) / STchain.get(refstatchain.get(r).toInt(), r))
//        }
//    }

    val Rchain = Matrix(Qchain.numRows, Qchain.numCols)
    for (i in 0..<Rchain.numRows) {
        for (j in 0..<Rchain.numCols) {
            if (Qchain.get(i, j) < GlobalConstants.Zero) {
                Rchain.set(i, j, 0.0)
            } else {
                Rchain.set(i, j, Qchain.get(i, j) / Tchain.get(i, j))
            }
        }
    }

    for (i in 0..<Xchain.numRows) {
        for (j in 0..<Xchain.numCols) {
            if (!Double.isFinite(Xchain.get(i, j))) Xchain.set(i, j, 0)
        }
    }

    for (i in 0..<Uchain.numRows) {
        for (j in 0..<Uchain.numCols) {
            if (!Double.isFinite(Uchain.get(i, j))) Uchain.set(i, j, 0)
        }
    }

    for (i in 0..<Qchain.numRows) {
        for (j in 0..<Qchain.numCols) {
            if (!Double.isFinite(Qchain.get(i, j))) Qchain.set(i, j, 0)
        }
    }

    for (i in 0..<Rchain.numRows) {
        for (j in 0..<Rchain.numCols) {
            if (!Double.isFinite(Rchain.get(i, j))) Rchain.set(i, j, 0)
        }
    }

    val Nzero = ArrayList<Int?>()
    for (i in 0..<Nchain.numCols) {
        if (Nchain.get(i) == 0.0) Nzero.add(i)
    }
    for (j in Nzero) {
        Xchain.set(j!!, 0.0)
        for (i in 0..<Uchain.numRows) {
            Uchain.set(i, j, 0.0)
        }
        for (i in 0..<Qchain.numRows) {
            Qchain.set(i, j, 0.0)
        }
        for (i in 0..<Rchain.numRows) {
            Rchain.set(i, j, 0.0)
        }
        for (i in 0..<Tchain.numRows) {
            Tchain.set(i, j, 0.0)
        }
        for (i in 0..<Wchain.numRows) {
            Wchain.set(i, j, 0.0)
        }
    }

    val ret2 =
        snDeaggregateChainResults(sn, Lchain, null, STchain, Vchain, alpha, null, null, Rchain, Tchain, null, Xchain)

    val res = MVAResult()
    res.QN = ret2.Q
    res.UN = ret2.U
    res.RN = ret2.R
    res.TN = ret2.T
    res.CN = ret2.C
    res.XN = ret2.X
    res.logNormConstAggr = lG
    return res
}

