/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.mva.handlers

import jline.api.pfqn.ld.pfqn_mvaldmx
import jline.api.sn.snDeaggregateChainResults
import jline.api.sn.snGetDemandsChain
import jline.lang.NetworkStruct
import jline.solvers.SolverOptions
import jline.solvers.mva.MVAResult
import jline.util.Utils
import jline.util.matrix.Matrix
import java.lang.Double


/**
 * Handler for the solver_mvald function.
 */
fun solver_mvald(sn: NetworkStruct, options: SolverOptions): MVAResult {
    val chainReturn = snGetDemandsChain(sn)
    val Lchain = chainReturn.Dchain
    val STchain = chainReturn.STchain
    val Vchain = chainReturn.Vchain
    val alpha = chainReturn.alpha
    val Nchain = chainReturn.Nchain
    val ST = Matrix(sn.rates.numRows, sn.rates.numCols)
    for (i in 0..<ST.numRows) {
        for (j in 0..<ST.numCols) {
            val serviceTime = 1.0 / sn.rates.get(i, j)
            if (!Double.isNaN(serviceTime)) {
                ST.set(i, j, serviceTime)
            }
        }
    }
    val M = STchain.numRows
    val C = sn.nchains
    val S = sn.nservers
    val N = sn.njobs
    var NchainSum = 0.0
    for (i in 0..<Nchain.numCols) {
        if (Double.isFinite(Nchain.get(i))) {
            NchainSum += Nchain.get(i)
        }
    }
    val mu_chain = Matrix.ones(M, NchainSum.toInt())
    for (i in 0..<M) {
        if (Utils.isInf(S.get(i))) {
            var j = 0
            while (j < NchainSum) {
                mu_chain.set(i, j, j + 1)
                j++
            }
        } else if (sn.lldscaling!=null && !sn.lldscaling.isEmpty) {
            var j = 0
            while (j < NchainSum) {
                mu_chain.set(i, j, sn.lldscaling.get(i, j))
                j++
            }
        } /* else {
            mu_chain(i,1:sum(Nchain)) = ones(1,sum(Nchain(isfinite(Nchain))));
            } is redundant, mu_chain is already set to Matrix.ones
            */
    }
    val lambda = Matrix(1, C)
    for (c in 0..<sn.nchains) {
        for (r in 0..<N.numCols) {
            if (Utils.isInf(N.get(r)) && sn.chains.get(c, r) == 1.0) {
                Nchain.set(0, c, kotlin.Double.Companion.POSITIVE_INFINITY)
                lambda.set(0, c, lambda.get(0, c) + 1.0 / ST.get(sn.refstat.get(r).toInt(), r))
            }
        }
    }
    val ret = pfqn_mvaldmx(lambda, Lchain, Nchain, Matrix(Nchain.numRows, Nchain.numCols), mu_chain, S)
    val Xchain = ret.X
    val Qchain = ret.Q
    val Uchain = ret.U
    val Tchain = Xchain.repmat(M, 1)
    for (i in 0..<Tchain.numRows) {
        for (j in 0..<Tchain.numCols) {
            Tchain.set(i, j, Tchain.get(i, j) * Vchain.get(i, j))
        }
    }
    val Rchain = Xchain.repmat(M, 1)
    for (i in 0..<Rchain.numRows) {
        for (j in 0..<Rchain.numCols) {
            Rchain.set(i, j, Qchain.get(i, j) / Rchain.get(i, j))
        }
    }
    val lG = kotlin.Double.Companion.NaN

    // This is likely wrong as it uses Little's law for the utilization computation
    val deAggregateReturn =
        snDeaggregateChainResults(sn, Lchain, null, STchain, Vchain, alpha, null, Uchain, Rchain, Tchain, null, Xchain)
    val iter = 1
    val res = MVAResult()
    res.QN = deAggregateReturn.Q
    res.UN = deAggregateReturn.U
    res.RN = deAggregateReturn.R
    res.TN = deAggregateReturn.T
    res.CN = deAggregateReturn.C
    res.XN = deAggregateReturn.X
    res.logNormConstAggr = lG
    res.iter = iter
    res.method = options.method
    return res
}