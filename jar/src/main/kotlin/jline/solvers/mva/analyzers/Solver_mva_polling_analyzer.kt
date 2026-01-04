/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.mva.analyzers

import jline.api.polling.polling_qsys_exhaustive
import jline.api.polling.polling_qsys_gated
import jline.api.polling.polling_qsys_1limited
import jline.api.mam.map_exponential
import jline.lang.NetworkStruct
import jline.lang.constant.NodeType
import jline.lang.constant.PollingType
import jline.lang.nodeparam.QueueNodeParam
import jline.lang.processes.Exp
import jline.lang.processes.Immediate
import jline.lang.processes.Distribution
import jline.solvers.SolverOptions
import jline.solvers.mva.MVAResult
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Converts a Distribution to a MAP representation.
 * Currently supports Exponential and Immediate distributions.
 */
private fun distributionToMAP(dist: Distribution?): MatrixCell {
    return when (dist) {
        is Exp -> map_exponential(dist.mean)
        is Immediate -> {
            val MAP = MatrixCell()
            val D0 = Matrix(1, 1)
            D0.set(0, 0, 0.0)
            val D1 = Matrix(1, 1)
            D1.set(0, 0, 0.0)
            MAP[0] = D0
            MAP[1] = D1
            MAP
        }
        null -> {
            val MAP = MatrixCell()
            val D0 = Matrix(1, 1)
            D0.set(0, 0, 0.0)
            val D1 = Matrix(1, 1)
            D1.set(0, 0, 0.0)
            MAP[0] = D0
            MAP[1] = D1
            MAP
        }
        else -> map_exponential(dist.mean)
    }
}

/**
 * MVA Polling System analyzer
 */
fun solver_mva_polling_analyzer(sn: NetworkStruct, options: SolverOptions): MVAResult {
    val startTime = System.nanoTime()
    var method = options.method
    val totiter = 1

    val QN = Matrix(sn.nstations, sn.nclasses)
    val UN = Matrix(sn.nstations, sn.nclasses)
    val RN = Matrix(sn.nstations, sn.nclasses)
    val TN = Matrix(sn.nstations, sn.nclasses)
    val CN = Matrix(sn.nstations, sn.nclasses)
    val AN = Matrix(sn.nstations, sn.nclasses)
    val WN = Matrix(sn.nstations, sn.nclasses)
    val XN = Matrix(sn.nstations, sn.nclasses)
    val lG = 0.0

    var source_ist = -1
    var queue_ist = -1
    var queue_node = -1
    for (i in sn.nodetype.indices) {
        if (sn.nodetype.get(i) == NodeType.Source) {
            source_ist = sn.nodeToStation.get(i).toInt()
        } else if (sn.nodetype.get(i) == NodeType.Queue) {
            queue_ist = sn.nodeToStation.get(i).toInt()
            queue_node = i
        }
    }

    if (source_ist == -1 || queue_ist == -1) {
        throw RuntimeException("Polling analyzer requires a Source and a Queue node")
    }

    if (sn.visits == null || sn.visits.size <= source_ist || sn.visits.get(source_ist) == null) {
        throw RuntimeException("Invalid visits matrix: source station $source_ist not found in visits matrix")
    }
    if (sn.stationToStateful == null || sn.stationToStateful.length() <= queue_ist) {
        throw RuntimeException("Invalid stationToStateful mapping: queue station $queue_ist not found")
    }

    val statefulIndex = sn.stationToStateful.get(queue_ist).toInt()
    val visitsMatrix = sn.visits.get(source_ist)
    if (visitsMatrix == null || visitsMatrix.length() <= statefulIndex) {
        throw RuntimeException("Invalid visits matrix: stateful index $statefulIndex not found in visits for station $source_ist")
    }

    val lambda = Matrix(1, sn.nclasses)
    val mu = Matrix(1, sn.nclasses)
    val k = sn.nservers.get(queue_ist).toInt()

    for (r in 0 until sn.nclasses) {
        lambda.set(0, r, sn.rates.get(source_ist) * visitsMatrix.get(statefulIndex, r))
        mu.set(0, r, sn.rates.get(queue_ist, r))
    }

    val sourceStation = sn.stations[source_ist]
    val queueStation = sn.stations[queue_ist]

    if (sourceStation == null || queueStation == null) {
        throw RuntimeException("Source or queue station not found in stations map")
    }

    val arvMAPs = Array<MatrixCell>(sn.nclasses) { r ->
        val jobClass = sn.jobclasses[r]
        sn.proc[sourceStation]?.get(jobClass) ?: throw RuntimeException("No arrival MAP for class $r")
    }

    val svcMAPs = Array<MatrixCell>(sn.nclasses) { r ->
        val jobClass = sn.jobclasses[r]
        sn.proc[queueStation]?.get(jobClass) ?: throw RuntimeException("No service MAP for class $r")
    }

    val switchMAPs = Array<MatrixCell>(sn.nclasses) { r ->
        val jobClass = sn.jobclasses[r]
        val queueNode = sn.nodes[queue_node]
        val nodeParam = sn.nodeparam[queueNode]

        if (nodeParam is QueueNodeParam && nodeParam.switchoverTime != null && nodeParam.switchoverTime.containsKey(jobClass)) {
            val switchoverDist = nodeParam.switchoverTime[jobClass]
            distributionToMAP(switchoverDist)
        } else {
            distributionToMAP(null)
        }
    }

    val queueNode = sn.nodes[queue_node]
    val nodeParam = sn.nodeparam[queueNode]
    var pollingType = PollingType.EXHAUSTIVE
    var pollingK = 1

    if (nodeParam is QueueNodeParam && nodeParam.pollingType != null) {
        pollingType = nodeParam.pollingType
        if (pollingType == PollingType.KLIMITED && nodeParam.pollingPar != null) {
            pollingK = nodeParam.pollingPar.toInt()
        }
    }

    if (method.equals("exact", ignoreCase = true)) {
        val ca = Matrix(1, sn.nclasses)
        for (r in 0 until sn.nclasses) {
            ca.set(0, r, Math.sqrt(sn.scv.get(source_ist, r)))
        }
        if (k == 1 && ca.elementMax() == 1.0 && ca.elementMin() == 1.0 &&
            (pollingType == PollingType.EXHAUSTIVE || pollingType == PollingType.GATED)) {
            method = "stationtime"
        } else {
            throw RuntimeException("MVA exact method unavailable for this model.")
        }
    }

    if (method == "default") {
        method = "stationtime"
    }

    val W = Matrix(1, sn.nclasses)
    val R = Matrix(1, sn.nclasses)

    when (method) {
        "stationtime" -> {
            val waitTimes = when (pollingType) {
                PollingType.EXHAUSTIVE -> polling_qsys_exhaustive(arvMAPs, svcMAPs, switchMAPs)
                PollingType.GATED -> polling_qsys_gated(arvMAPs, svcMAPs, switchMAPs)
                PollingType.KLIMITED -> {
                    if (pollingK == 1) {
                        polling_qsys_1limited(arvMAPs, svcMAPs, switchMAPs)
                    } else {
                        throw RuntimeException("MVA method unavailable for K-limited polling with K>1.")
                    }
                }
                else -> throw RuntimeException("Unsupported polling type: $pollingType")
            }
            for (r in 0 until sn.nclasses) {
                W.set(0, r, waitTimes[r])
                R.set(0, r, waitTimes[r] + 1.0 / mu.get(0, r))
            }
        }
        else -> throw RuntimeException("Unsupported polling solution method: $method")
    }

    for (r in 0 until sn.nclasses) {
        val visitRatio = visitsMatrix.get(statefulIndex, r)
        RN.set(queue_ist, r, R.get(0, r) * visitRatio)
        CN.set(queue_ist, r, R.get(0, r))
        XN.set(queue_ist, r, lambda.get(0, r))
        UN.set(queue_ist, r, lambda.get(0, r) / mu.get(0, r) / k)
        TN.set(source_ist, r, lambda.get(0, r))
        TN.set(queue_ist, r, lambda.get(0, r))
        QN.set(queue_ist, r, XN.get(queue_ist, r) * RN.get(queue_ist, r))
    }

    val res = MVAResult()
    res.method = method
    res.QN = QN
    res.RN = RN
    res.XN = XN
    res.UN = UN
    res.TN = TN
    res.CN = CN
    res.AN = AN
    res.WN = WN
    res.runtime = (System.nanoTime() - startTime) / 1000000000.0
    res.iter = totiter
    res.logNormConstAggr = lG

    return res
}
