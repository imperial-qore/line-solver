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
import jline.solvers.mam.handlers.solver_mam_retrial
import jline.solvers.mam.handlers.solver_mna_open
import jline.solvers.mam.handlers.solver_mna_closed
import jline.api.qsys.qsys_is_retrial
import jline.api.qsys.RetrialInfo
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
        // Check if network is a valid BMAP/PH/N/N bufferless retrial queue
        val retInfo = try { qsys_is_retrial(sn) } catch (e: Exception) { RetrialInfo() }

        if (retInfo.isRetrial) {
            // Use BMAP/PH/N/N retrial solver
            if (options.method == "default") {
                line_debug(options.verbose, "Default method: using retrial for BMAP/PH/N/N bufferless topology")
            }
            line_debug(options.verbose, "Detected BMAP/PH/N/N retrial topology, using retrial method")
            result = solver_mam_retrial(sn, options)
            result.method = "retrial"
        } else {
            // arrival process per chain rescaled by visits at each node
            if (options.method == "default") {
                line_debug(options.verbose, "Default method: using dec.source")
            }
            line_debug(options.verbose, "Using default/dec.source method, calling solver_mam_basic")
            result = solver_mam_basic(sn, options)
            result.method = "dec.source"
        }
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
    } else if (options.method == "retrial") {
        line_debug(options.verbose, "Using retrial method for BMAP/PH/N/N bufferless topology")
        result = solver_mam_retrial(sn, options)
        result.method = "retrial"
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

    // Handle self-looping classes: override metrics if the method produced
    // incorrect values (Inf/0). Mirrors MATLAB solver_mam_ag.m lines 675-725.
    val M = sn.nstations
    val K = sn.nclasses
    if (sn.isslc != null) {
        var hasSlc = false
        for (r in 0..<K) {
            if (sn.isslc.get(r) == 1.0) {
                hasSlc = true
                break
            }
        }
        if (hasSlc) {
            for (r in 0..<K) {
                if (sn.isslc.get(r) == 1.0) {
                    val refst = sn.refstat.get(r).toInt()
                    if (refst in 0..<M) {
                        val qVal = result.QN.get(refst, r)
                        // Override if QN is 0, Inf, or NaN (method didn't compute correctly)
                        if (qVal == 0.0 || java.lang.Double.isInfinite(qVal) || java.lang.Double.isNaN(qVal)) {
                            // Clear SLC metrics at all stations first
                            for (i in 0..<M) {
                                result.QN.set(i, r, 0.0)
                                result.UN.set(i, r, 0.0)
                                result.TN.set(i, r, 0.0)
                                result.RN.set(i, r, 0.0)
                            }
                            // All jobs stay at reference station
                            result.QN.set(refst, r, sn.njobs.get(r))
                            val muIr = sn.rates.get(refst, r)
                            if (!java.lang.Double.isNaN(muIr) && muIr > 0) {
                                val nservers = sn.nservers.get(refst)
                                if (java.lang.Double.isInfinite(nservers)) {
                                    // Delay (infinite server)
                                    result.UN.set(refst, r, result.QN.get(refst, r))
                                    result.TN.set(refst, r, muIr * result.QN.get(refst, r))
                                } else {
                                    // Queue (finite server): use remaining capacity
                                    var otherUtil = 0.0
                                    for (s in 0..<K) {
                                        if (s != r && sn.isslc.get(s) != 1.0) {
                                            otherUtil += result.UN.get(refst, s)
                                        }
                                    }
                                    val remainingCapacity = Math.max(0.0, 1.0 - otherUtil)
                                    val slcDemand = result.QN.get(refst, r) / muIr
                                    result.UN.set(refst, r, Math.min(slcDemand, remainingCapacity))
                                    result.TN.set(refst, r, muIr * result.UN.get(refst, r))
                                }
                            }
                        }
                    }
                }
            }
            // Recompute response times from Little's law for all entries
            for (i in 0..<M) {
                for (r in 0..<K) {
                    if (result.TN.get(i, r) > 0) {
                        result.RN.set(i, r, result.QN.get(i, r) / result.TN.get(i, r))
                    } else {
                        result.RN.set(i, r, 0.0)
                    }
                }
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

