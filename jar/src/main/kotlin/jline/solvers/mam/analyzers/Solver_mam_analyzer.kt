package jline.solvers.mam.analyzers

import jline.io.line_debug
import jline.lang.NetworkStruct
import jline.lang.constant.SchedStrategy
import jline.solvers.SolverOptions
import jline.solvers.mam.MAMResult
import jline.solvers.mam.handlers.solver_mam
import jline.solvers.mam.handlers.solver_mam_basic
import jline.solvers.mam.handlers.solver_mam_ldqbd
import jline.solvers.mam.handlers.solver_mam_ag
import jline.solvers.mam.handlers.solver_mna_open
import jline.solvers.mam.handlers.solver_mna_closed
import jline.api.sn.snIsOpenModel
import jline.api.sn.snIsClosedModel
import jline.api.sn.snNonmarkovToPh
import kotlin.RuntimeException

fun solver_mam_analyzer(snInput: NetworkStruct, options: SolverOptions): MAMResult {
    val start = System.nanoTime()

    // Convert non-Markovian distributions to PH
    val sn = snNonmarkovToPh(snInput, options)

    // Check if the model is mixed (has both open and closed classes)
    val isOpen = snIsOpenModel(sn)
    val isClosed = snIsClosedModel(sn)
    val isMixed = isOpen && isClosed
    
    // Mixed models are supported by the dec.source method
    
    options.config.merge = "super"
    options.config.compress = "mixture.order1"
    options.config.space_max = 128
    var result = MAMResult()

    if (options.method == "dec.mmap") {
        line_debug(options.verbose, "Using dec.mmap method, calling solver_mam")
        result = solver_mam(sn, options)
        result.method = "dec.mmap"
    } else if (options.method == "default" || options.method == "dec.source") {
        line_debug(options.verbose, "Using default/dec.source method, calling solver_mam_basic")
        result = solver_mam_basic(sn, options)
        result.method = "dec.source"
    } else if (options.method == "dec.poisson") {
        line_debug(options.verbose, "Using dec.poisson method with space_max=1, calling solver_mam_basic")
        options.config.space_max = 1
        result = solver_mam_basic(sn, options)
        result.method = "dec.poisson"
    } else if (options.method == "mna") {
        if (snIsClosedModel(sn)) {
            line_debug(options.verbose, "MNA method (closed)")
            result = solver_mna_closed(sn, options)
            result.method = "mna"
        } else if (snIsOpenModel(sn)) {
            line_debug(options.verbose, "MNA method (open)")
            result = solver_mna_open(sn, options) as MAMResult
            result.method = "mna"
        } else {
            throw RuntimeException("The mna method in SolverMAM does not support mixed models.")
        }
    } else if (options.method == "inap" || options.method == "inapplus" || options.method == "exact") {
        line_debug(options.verbose, "Using RCAT method: ${options.method}")
        result = solver_mam_ag(sn, options)
        result.method = options.method
    } else if (options.method == "ldqbd") {
        line_debug(options.verbose, "Using LDQBD method for single-class closed network")
        result = solver_mam_ldqbd(sn, options)
        result.method = "ldqbd"
    } else {
        throw RuntimeException("Unknown method")
    }

    for (i in 0..<sn.nstations) {
        if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.EXT) {
            for (j in 0..<result.TN.numCols) {
                result.TN.set(i, j, sn.rates.get(i, j))
            }
        }
    }

    result.QN.setNaNToZero()
    result.CN.setNaNToZero()
    result.RN.setNaNToZero()
    result.UN.setNaNToZero()
    result.XN.setNaNToZero()
    result.TN.setNaNToZero()
    val finish = System.nanoTime()
    result.runtime = (finish - start).toDouble() / 1000000000.0

    return result
}

