/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.mva.analyzers

import jline.api.qsys.*
import jline.io.Ret
import jline.lang.NetworkStruct
import jline.lang.constant.NodeType
import jline.lang.constant.ProcessType
import jline.solvers.SolverOptions
import jline.solvers.mva.MVAResult
import jline.util.SerializableFunction
import jline.util.matrix.Matrix
import org.apache.commons.math3.analysis.UnivariateFunction
import org.apache.commons.math3.analysis.solvers.BrentSolver
import org.apache.commons.math3.analysis.solvers.UnivariateSolver
import org.apache.commons.math3.util.FastMath

/**
 * MVA Query System analyzer
 */
fun solver_mva_qsys_analyzer(sn: NetworkStruct, options: SolverOptions): MVAResult {
    var res = MVAResult()
    val startTime = System.nanoTime()
    var method = options.method
    val QN = Matrix(sn.nstations, sn.nclasses)
    val UN = Matrix(sn.nstations, sn.nclasses)
    val RN = Matrix(sn.nstations, sn.nclasses)
    val TN = Matrix(sn.nstations, sn.nclasses)
    val CN = Matrix(sn.nstations, sn.nclasses)
    val AN = Matrix(sn.nstations, sn.nclasses)
    val WN = Matrix(sn.nstations, sn.nclasses)
    val XN = Matrix(sn.nstations, sn.nclasses)
    var lG = Double.NaN
    val it = 1

    var source_ist = -1
    var queue_ist = -1
    for (i in sn.nodetype.indices) {
        if (sn.nodetype.get(i) == NodeType.Source) {
            source_ist = sn.nodeToStation.get(i).toInt()
        } else if (sn.nodetype.get(i) == NodeType.Queue) {
            queue_ist = sn.nodeToStation.get(i).toInt()
        }
    }

    // For single-chain qsys models, use chain index 0
    val chainIdx = 0

    // Check if visits matrix exists and has proper dimensions
    if (sn.visits == null || sn.visits.isEmpty() || sn.visits.get(chainIdx) == null) {
        throw RuntimeException("Invalid visits matrix: chain $chainIdx not found in visits matrix")
    }
    if (sn.stationToStateful == null || sn.stationToStateful.length() <= queue_ist) {
        throw RuntimeException("Invalid stationToStateful mapping: queue station $queue_ist not found")
    }

    val statefulIndex = sn.stationToStateful.get(queue_ist).toInt()
    val visitsMatrix = sn.visits.get(chainIdx)
    if (visitsMatrix == null || visitsMatrix.length() <= statefulIndex) {
        throw RuntimeException("Invalid visits matrix: stateful index $statefulIndex not found in visits for chain $chainIdx")
    }

    val lambda = sn.rates.get(source_ist) * visitsMatrix.get(statefulIndex)
    val k = sn.nservers.get(queue_ist).toInt()
    val mu = sn.rates.get(queue_ist)
    val ca = FastMath.sqrt(sn.scv.get(source_ist))
    val cs = FastMath.sqrt(sn.scv.get(queue_ist))

    // Check if this is a BMAP arrival process (MX/M/1)
    val sourceStation = sn.stations[source_ist]
    val jobClass = sn.jobclasses[0]  // Single class model
    val sourceProcType = sn.procid?.get(sourceStation)?.get(jobClass)
    val isBMAP = sourceProcType == ProcessType.BMAP

    if (method == "exact") {
        if (isBMAP && cs == 1.0 && k == 1) {
            // MX/M/1: Batch Poisson arrivals with exponential service
            method = "mxm1"
        } else if (ca == 1.0 && cs == 1.0 && k == 1) {
            method = "mm1"
        } else if (ca == 1.0 && cs == 1.0 && k > 1) {
            method = "mmk"
        } else if (ca == 1.0 && k == 1) {
            method = "mg1"
        } else if (cs == 1.0 && k == 1) {
            method = "gm1"
        } else {
            throw RuntimeException("MVA exact method unavailable for this model.")
        }
    }
    if (method == "default") {
        if (isBMAP && cs == 1.0 && k == 1) {
            // Default for BMAP with exponential service is MX/M/1
            method = "mxm1"
        } else if (k > 1) {
            method = "gigk"
        } else {
            method = "gig1.klb"
        }
    }
    var R = 0.0
    var lambdaEffective = lambda  // Effective job arrival rate (may differ from batch rate for BMAP)
    when (method) {
        "mm1" -> {
            qsys_mm1(lambda, mu)
            R = Ret.qsys.W
        }

        "mxm1" -> {
            // MX/M/1: Batch Poisson arrivals with exponential service
            // Extract batch size distribution from BMAP process
            val proc = sn.proc?.get(sourceStation)?.get(jobClass)
            if (proc == null) {
                throw RuntimeException("BMAP process not found for source station")
            }

            // BMAP structure: proc[0]=D0, proc[1]=D1_total, proc[2+k]=D_{k+1} (batch size k+1)
            val D0 = proc.get(0)
            val maxBatchSize = proc.size() - 2

            // For single-phase BMAP, batch arrival rate is -D0[0,0]
            val lambdaBatch = -D0.get(0, 0)

            // Calculate E[X] and E[X²] from batch matrices
            // For single-phase: rate_k proportional to D_k[0,0]
            var meanBatchSize = 0.0
            var secondMomentBatchSize = 0.0
            var totalRate = 0.0

            for (idx in 2 until proc.size()) {
                val batchSize = idx - 1  // D_k at index k+1 means batch size k
                val Dk = proc.get(idx)
                val rate_k = Dk.get(0, 0)  // Single phase: just the (0,0) element
                totalRate += rate_k
                meanBatchSize += batchSize * rate_k
                secondMomentBatchSize += batchSize * batchSize * rate_k
            }

            // Normalize by total rate to get moments
            if (totalRate > 0) {
                meanBatchSize /= totalRate
                secondMomentBatchSize /= totalRate
            }

            qsys_mxm1(lambdaBatch, mu, meanBatchSize, secondMomentBatchSize)
            R = Ret.qsys.W

            // Set effective job arrival rate for result metrics
            lambdaEffective = lambdaBatch * meanBatchSize

            if (options.method == "default") {
                method = "default [mxm1]"
            }
        }

        "mmk" -> {
            qsys_mmk(lambda, mu, k)
            R = Ret.qsys.W
        }

        "mg1", "mgi1" -> {
            qsys_mg1(lambda, mu, cs)
            R = Ret.qsys.W
        }

        "gigk" -> {
            qsys_gigk_approx(lambda, mu, ca, cs, k)
            R = Ret.qsys.W
        }

        "gigk.kingman_approx" -> {
            qsys_gigk_approx_kingman(lambda, mu, ca, cs, k)
            R = Ret.qsys.W
        }

        "gig1", "gig1.kingman" -> {
            qsys_gig1_ubnd_kingman(lambda, mu, ca, cs)
            R = Ret.qsys.W
        }

        "gig1.heyman" -> {
            qsys_gig1_approx_heyman(lambda, mu, ca, cs)
            R = Ret.qsys.W
        }

        "gig1.allen" -> {
            qsys_gig1_approx_allencunneen(lambda, mu, ca, cs)
            R = Ret.qsys.W
        }

        "gig1.kobayashi" -> {
            qsys_gig1_approx_kobayashi(lambda, mu, ca, cs)
            R = Ret.qsys.W
        }

        "gig1.klb" -> {
            qsys_gig1_approx_klb(lambda, mu, ca, cs)
            R = Ret.qsys.W
            if (options.method == "default") {
                method = "default [gig1.klb]"
            }
        }

        "gig1.marchal" -> {
            qsys_gig1_approx_marchal(lambda, mu, ca, cs)
            R = Ret.qsys.W
        }

        "gm1", "gim1" -> {
            // sigma = Load at arrival instants (Laplace transform of the inter-arrival times)
            @Suppress("UNCHECKED_CAST")
            val F = sn.lst[sn.stations[source_ist]] as SerializableFunction<Double, Double>
            val LA = UnivariateFunction { s: Double -> F.apply(s) }
            // Define the function for fzero
            val func = UnivariateFunction { x: Double -> LA.value(mu - mu * x) - x }
            // Use BrentSolver to find the root of the function
            val solver: UnivariateSolver = BrentSolver()
            val sigma = solver.solve(1000, func, 0.0, 1.0) // Adjust the max evaluations and range as needed
            qsys_gm1(sigma, mu)
            R = Ret.qsys.W
        }

        else -> throw RuntimeException("Unsupported method for a model with 1 station and 1 class.")
    }
    // Follow MATLAB implementation pattern
    // Use lambdaEffective for job-level metrics (differs from batch rate for BMAP)
    for (r in 0 until sn.nclasses) {
        RN.set(queue_ist, r, R * visitsMatrix.get(statefulIndex))
        CN.set(queue_ist, r, RN.get(queue_ist, r))
        XN.set(queue_ist, r, lambdaEffective)
        UN.set(queue_ist, r, lambdaEffective / mu / k)
        TN.set(source_ist, r, lambdaEffective)
        TN.set(queue_ist, r, lambdaEffective)
        QN.set(queue_ist, r, XN.get(queue_ist, r) * RN.get(queue_ist, r))
    }
    lG = 0.0
    val endTime = System.nanoTime()
    res.QN = QN
    res.UN = UN
    res.RN = RN
    res.TN = TN
    res.CN = CN
    res.XN = XN
    res.AN = AN
    res.WN = WN
    res.logNormConstAggr = lG
    res.runtime = (endTime - startTime) / 1000000000.0
    res.iter = it
    res.method = method
    return res
}