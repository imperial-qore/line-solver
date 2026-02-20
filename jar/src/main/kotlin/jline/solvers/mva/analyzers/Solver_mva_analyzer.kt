/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.mva.analyzers

import jline.api.sn.snHasFractionalPopulations
import jline.api.sn.snHasProductForm
import jline.io.line_warning
import jline.io.mfilename
import jline.lang.NetworkStruct
import jline.VerboseLevel
import jline.solvers.SolverOptions
import jline.solvers.mva.MVAResult
import jline.solvers.mva.handlers.solver_amva
import jline.solvers.mva.handlers.solver_mva
import jline.solvers.mva.handlers.solver_qna

/**
 * MVA Analyzer class
 */
fun solver_mva_analyzer(sn: NetworkStruct, options: SolverOptions): MVAResult {
    var res = MVAResult()
    var method = options.method
    method = method.replaceFirst("^amva\\.".toRegex(), "")

    val startTime = System.nanoTime()
    var ret: MVAResult? = null
    when (method) {
        "exact", "mva" -> {
            ret = solver_mva(sn, options)
            ret.iter = 0
        }

        "qna" -> if (options.verbose != VerboseLevel.SILENT) ret = solver_qna(sn, options)
        "default" -> {
            if (sn.nchains <= 4 && sn.njobs.sumRows()
                .toDouble() <= 20 && snHasProductForm(sn) && !snHasFractionalPopulations(sn)) {
                // The parameters above take in the worst case a handful of ms
                ret = solver_mva(sn, options)
                method = "exact"
            } else {
                ret = solver_amva(sn, options)
                method = ret.method
            }
        }

        "amva", "bs", "qd", "qli", "fli", "lin", "qdlin", "sqni", "gflin", "egflin",
        "schmidt", "schmidt-ext", "ab" -> {
            ret = solver_amva(sn, options)
            method = ret.method
        }

        else -> {
            res = MVAResult()
            res.QN = null
            res.UN = null
            res.RN = null
            res.TN = null
            res.CN = null
            res.XN = null
            res.AN = null
            res.WN = null
            res.logNormConstAggr = null!!.toDouble()
            res.iter = null!!.toInt()
            res.runtime = (System.nanoTime() - startTime) / 1000000000.0
            if (options.verbose != VerboseLevel.SILENT) line_warning(mfilename(object : Any() {}),
                "Unsupported SolverMVA method.")
            return res
        }
    }
    val endTime = System.nanoTime()
    res.QN = ret!!.QN
    res.UN = ret.UN
    res.RN = ret.RN
    res.TN = ret.TN
    res.CN = ret.CN
    res.XN = ret.XN
    res.AN = ret.AN
    res.WN = ret.WN
    res.logNormConstAggr = ret.logNormConstAggr
    res.iter = ret.iter
    res.method = method
    res.runtime = (endTime - startTime) / 1000000000.0
    return res
}
