/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.mva.analyzers

import jline.lang.NetworkStruct
import jline.solvers.SolverOptions
import jline.solvers.mva.MVAResult
import jline.solvers.mva.handlers.solver_amva
import jline.solvers.mva.handlers.solver_mvald

/**
 * MVALD Analyzer
 */
fun solver_mvald_analyzer(sn: NetworkStruct, options: SolverOptions): MVAResult {
    var res = MVAResult()
    val startTime = System.nanoTime()
    val method = options.method
    var ret: MVAResult? = null
    when (method) {
        "exact", "mva" -> {
            if (sn.cdscaling != null && !sn.cdscaling.isEmpty()) {
                throw RuntimeException("Exact class-dependent solver not available in MVA.")
            }
            ret = solver_mvald(sn, options)
        }

        "default", "amva", "qd", "lin", "qdlin" -> ret = solver_amva(sn, options)
        else -> throw RuntimeException("The " + method + " method is not supported by the load-dependent MVA solver.")
    }
    val endTime = System.nanoTime()
    res.QN = ret.QN
    res.UN = ret.UN
    res.RN = ret.RN
    res.TN = ret.TN
    res.CN = ret.CN
    res.XN = ret.XN
    res.AN = ret.AN
    res.WN = ret.WN
    res.logNormConstAggr = ret.logNormConstAggr
    res.runtime = (endTime - startTime) / 1000000000.0
    res.iter = ret.iter
    res.method = method
    return res
}

