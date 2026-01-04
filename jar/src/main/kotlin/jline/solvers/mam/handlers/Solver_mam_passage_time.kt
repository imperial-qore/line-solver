package jline.solvers.mam.handlers

import jline.api.mam.*
import jline.api.map.MAPM1PSCdfRespT
import jline.lang.NetworkStruct
import jline.GlobalConstants
import jline.lang.constant.SchedStrategy
import jline.lib.butools.MMAPPH1FCFS
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import kotlin.math.sqrt
import jline.io.mfilename
import jline.io.line_warning
import jline.io.line_error

/**
 * Analyzes class priorities and determines if the configuration is supported.
 *
 * @param sn NetworkStruct containing class priorities
 * @return PriorityAnalysis result indicating support level
 */
data class PriorityAnalysis(
    val isIdentical: Boolean,
    val isAllDistinct: Boolean,
    val isSupported: Boolean,
    val message: String
)

fun analyzePriorities(sn: NetworkStruct, context: String = "basic"): PriorityAnalysis {
    val priorities = sn.classprio
    val firstPriority = priorities.get(0)
    var identical = true

    for (i in 1..<priorities.length()) {
        if (priorities.get(i) != firstPriority) {
            identical = false
            break
        }
    }

    if (identical) {
        return PriorityAnalysis(
            isIdentical = true,
            isAllDistinct = false,
            isSupported = true,
            message = "Identical priorities - fully supported"
        )
    }

    // Check if all priorities are distinct
    val uniquePriorities = mutableSetOf<Double>()
    for (i in 0..<priorities.length()) {
        uniquePriorities.add(priorities.get(i))
    }

    val allDistinct = uniquePriorities.size == priorities.length()

    return if (allDistinct) {
        if (context == "passage_time") {
            PriorityAnalysis(
                isIdentical = false,
                isAllDistinct = true,
                isSupported = false,
                message = "Response time distribution in priority models not yet supported"
            )
        } else {
            PriorityAnalysis(
                isIdentical = false,
                isAllDistinct = true,
                isSupported = true,
                message = "All distinct priorities - using NPPR analysis"
            )
        }
    } else {
        PriorityAnalysis(
            isIdentical = false,
            isAllDistinct = false,
            isSupported = false,
            message = "Mixed priority configuration not supported - requires identical or all distinct priorities"
        )
    }
}

fun solver_mam_passage_time(sn: NetworkStruct,
                            PH: Map<*, Map<*, MatrixCell>>,
                            options: jline.solvers.SolverOptions): Map<Int, MatrixCell> {
    val RD: MutableMap<Int, MatrixCell> = HashMap<Int, MatrixCell>()
    val M = sn.nstations
    val K = sn.nclasses
    val N = sn.njobs.transpose()

    // Check if model is open (all jobs are infinite)
    var open = true
    for (i in 0..<N.length()) {
        if (java.lang.Double.isFinite(N.get(i))) {
            open = false
            break
        }
    }

    var A: MatrixCell? = MatrixCell()
    var idx_arv = 0
    var idx_q = 0
    var is_ps = false

    if (M == 2 && open) {
        val pie: MutableMap<Int, Matrix> = HashMap<Int, Matrix>()
        val S: MutableMap<Int, Matrix> = HashMap<Int, Matrix>()

        for (i in 0..<M) {
            val station = sn.stations.get(i)

            if (sn.sched.get(station) == SchedStrategy.EXT) {
                // Process external (arrival) station
                val stationProc = PH.get(station) as? Map<*, MatrixCell>
                if (stationProc != null) {
                    val jobClass0 = sn.jobclasses.get(0)
                    val proc0 = stationProc.get(jobClass0) as? MatrixCell
                    if (proc0 != null) {
                        A!!.set(0, proc0.get(0)) // D0 matrix
                        A.set(1, proc0.get(1))   // D1 matrix
                        A.set(2, proc0.get(1))   // D1 matrix (duplicate for superposition)

                        for (k in 1..<K) {
                            val jobClassK = sn.jobclasses.get(k)
                            val procK = stationProc.get(jobClassK) as? MatrixCell
                            if (procK != null) {
                                val B = MatrixCell()
                                B.set(0, procK.get(0)) // D0 matrix
                                B.set(1, procK.get(1)) // D1 matrix
                                B.set(2, procK.get(1)) // D1 matrix
                                A = mmap_super(A!!, B)
                            }
                        }
                    }
                }
                idx_arv = i

            } else if (sn.sched.get(station) == SchedStrategy.FCFS || sn.sched.get(station) == SchedStrategy.HOL) {
                // Process service station
                val stationProc = PH.get(station) as? Map<*, MatrixCell>
                if (stationProc != null) {
                    for (k in 0..<K) {
                        val jobClass = sn.jobclasses.get(k)
                        val procK = stationProc.get(jobClass) as? MatrixCell
                        if (procK != null) {
                            // Scale by number of servers like MATLAB version
                            val scaledMean = map_mean(procK.get(0), procK.get(1)) / sn.nservers.get(i)
                            val scaledProc = map_scale(procK.get(0), procK.get(1), scaledMean)
                            pie.put(k, map_pie(scaledProc.get(0), scaledProc.get(1)))
                            S.put(k, scaledProc.get(0))
                        }
                    }
                }
                idx_q = i
                is_ps = false

            } else if (sn.sched.get(station) == SchedStrategy.PS) {
                // Processor-sharing queue
                val stationProc = PH.get(station) as? Map<*, MatrixCell>
                if (stationProc != null) {
                    for (k in 0..<K) {
                        val jobClass = sn.jobclasses.get(k)
                        val procK = stationProc.get(jobClass) as? MatrixCell
                        if (procK != null) {
                            // For PS, service times remain exponential (no scaling needed for distribution)
                            S.put(k, procK.get(0))  // Service process subgenerator
                            pie.put(k, map_pie(procK.get(0), procK.get(1)))
                        }
                    }
                }
                idx_q = i
                is_ps = true

            } else {
                throw RuntimeException("Unsupported scheduling strategy")
            }
        }

        val priorityAnalysis = analyzePriorities(sn, "passage_time")

        if (!priorityAnalysis.isSupported) {
            if (priorityAnalysis.isAllDistinct) {
                // Future NPPR support - warn but don't throw exception
                line_warning(mfilename(object {}), priorityAnalysis.message)
                return RD // Return empty result for graceful degradation
            } else {
                // Invalid configuration - throw exception
                throw RuntimeException(priorityAnalysis.message)
            }
        } else if (is_ps) {
            // MAP/M/1-PS queue: use specialized sojourn time distribution algorithm
            // Check that service processes are exponential (M)
            for (k in 0..<K) {
                if (S.get(k)!!.getNumRows() != 1) {
                    line_error(mfilename(object {}), "PS queue requires exponential (Markovian) service times")
                }
            }

            // For single-class case, directly compute sojourn time CDF
            if (K == 1) {
                // Extract MAP arrival parameters
                val C_map = A!!.get(0)  // MAP C matrix
                val D_map = A.get(1)    // MAP D matrix

                // Extract exponential service rate
                val mu = -S.get(0)!!.get(0, 0)  // Service rate (S is negative generator)

                // Estimate CDF range
                val lambda = map_lambda(C_map, D_map)
                val rho = lambda / mu
                val approx_mean = 1.0 / (mu * (1.0 - rho))  // Approximate M/M/1-PS mean

                // Generate x points (use similar range as FCFS case)
                val n_pts = options.config.num_cdf_pts
                val x_max = approx_mean * 10.0  // Conservative upper bound
                val x_vals = DoubleArray(n_pts)
                for (i in 0..<n_pts) {
                    x_vals[i] = (x_max * i) / (n_pts - 1).toDouble()
                }

                // Call MAP/M/1-PS sojourn time CDF function
                val W_bar = MAPM1PSCdfRespT.computeCdf(C_map, D_map, mu, x_vals)

                // Convert to CDF format (W_bar is complementary CDF)
                val F = Matrix(n_pts, 1)
                val X = Matrix(n_pts, 1)
                for (i in 0..<n_pts) {
                    F.set(i, 0, 1.0 - W_bar[i])
                    X.set(i, 0, x_vals[i])
                }

                // Store results
                if (!RD.containsKey(idx_arv)) {
                    RD.put(idx_arv, MatrixCell())
                }
                RD.get(idx_arv)!!.set(0, Matrix(0, 0))

                if (!RD.containsKey(idx_q)) {
                    RD.put(idx_q, MatrixCell())
                }
                RD.get(idx_q)!!.set(0, Matrix.concatColumns(F, X, null))

            } else {
                // Multi-class PS: aggregate arrivals into single MAP
                // Extract individual service rates
                val mu_vec = DoubleArray(K)
                for (k in 0..<K) {
                    mu_vec[k] = -S.get(k)!!.get(0, 0)
                }

                // Check if all service rates are equal (required for current implementation)
                var allEqual = true
                for (k in 1..<K) {
                    if (Math.abs(mu_vec[k] - mu_vec[0]) > GlobalConstants.FineTol) {
                        allEqual = false
                        break
                    }
                }
                if (!allEqual) {
                    line_error(mfilename(object {}), "Multi-class PS currently requires identical service rates")
                }
                val mu = mu_vec[0]

                // Estimate CDF range using aggregated arrival rate
                val C_map = A!!.get(0)
                var D_map_sum = A.get(1).copy()
                for (i in 2..<A.size()) {
                    D_map_sum = D_map_sum.add(1.0, A.get(i))
                }

                val lambda = map_lambda(C_map, D_map_sum)
                val rho = lambda / mu
                val approx_mean = 1.0 / (mu * (1.0 - rho))

                val n_pts = options.config.num_cdf_pts
                val x_max = approx_mean * 10.0
                val x_vals = DoubleArray(n_pts)
                for (i in 0..<n_pts) {
                    x_vals[i] = (x_max * i) / (n_pts - 1).toDouble()
                }

                // Compute sojourn time for aggregated system
                val W_bar = MAPM1PSCdfRespT.computeCdf(C_map, D_map_sum, mu, x_vals)

                // Convert to CDF
                val F = Matrix(n_pts, 1)
                val X = Matrix(n_pts, 1)
                for (i in 0..<n_pts) {
                    F.set(i, 0, 1.0 - W_bar[i])
                    X.set(i, 0, x_vals[i])
                }

                // Assign same CDF to all classes (approximation)
                for (k in 0..<K) {
                    if (!RD.containsKey(idx_arv)) {
                        RD.put(idx_arv, MatrixCell())
                    }
                    RD.get(idx_arv)!!.set(k, Matrix(0, 0))

                    if (!RD.containsKey(idx_q)) {
                        RD.put(idx_q, MatrixCell())
                    }
                    RD.get(idx_q)!!.set(k, Matrix.concatColumns(F, X, null))
                }
            }
        } else {
            // FCFS/HOL queue: use existing MMAPPH1FCFS method
            // Prepare arrival process for MMAPPH1FCFS like MATLAB: {A{[1,3:end]}}
            // This means take elements [0, 2, 3, ..., n-1] (skip element 1)
            if (A != null && A.size() > 2) {
                val newA = MatrixCell()
                newA.set(0, A.get(0))  // Keep first element (D0)
                for (i in 2..<A.size()) {
                    newA.set(i - 1, A.get(i))  // Skip element 1, shift rest down
                }
                A = newA
            }

            // Call MMAPPH1FCFS to get phase-type response time distributions
            val mmapResult = MMAPPH1FCFS(A!!,
                pie.mapKeys { it.key as Int? }.mapValues { it.value },
                S.mapKeys { it.key as Int? }.mapValues { it.value },
                null, null, null, null, false, true, null, null)

            val alpha: Map<Int, Matrix> = mmapResult.get("stDistrPH_alpha")?.let { distr ->
                @Suppress("UNCHECKED_CAST")
                distr as? Map<Int, Matrix>
            } ?: emptyMap()

            val D0: Map<Int, Matrix> = mmapResult.get("stDistrPH_A")?.let { d0Map ->
                @Suppress("UNCHECKED_CAST")
                d0Map as? Map<Int, Matrix>
            } ?: emptyMap()

            // Process each class to generate CDF
            for (k in 0..<K) {
                if (alpha.containsKey(k) && D0.containsKey(k)) {
                    val alphaK = alpha.get(k)!!
                    val D0K = D0.get(k)!!

                    // Create phase-type representation: {D0, D1}
                    val negD0K = D0K.copy()
                    negD0K.scaleEq(-1.0)
                    val D1K = negD0K.mult(Matrix.ones(alphaK.length(), 1)).mult(alphaK.transpose())

                    // Calculate variance and mean for CDF range determination
                    val variance = map_var(D0K, D1K)
                    val meanResp = map_mean(D0K, D1K)
                    val sigma = sqrt(variance)

                    // Determine CDF range (like MATLAB: while map_cdf < 1-FineTol)
                    var n = 5
                    var maxTime = meanResp + n * sigma
                    while (map_cdf(D0K, D1K, Matrix.singleton(maxTime)).get(0) < 1 - GlobalConstants.FineTol) {
                        n++
                        maxTime = meanResp + n * sigma
                    }

                    // Generate CDF points
                    val n_pts = options.config.num_cdf_pts
                    val F = Matrix(n_pts, 1)
                    val X = Matrix(n_pts, 1)

                    for (i in 0..<n_pts) {
                        val t = (maxTime * i) / (n_pts - 1).toDouble()
                        X.set(i, 0, t)
                        F.set(i, 0, map_cdf(D0K, D1K, Matrix.singleton(t)).get(0))
                    }

                    // Store results matching MATLAB format: RD{idx_arv,k} = [], RD{idx_q,k} = [F,X]
                    if (!RD.containsKey(idx_arv)) {
                        RD.put(idx_arv, MatrixCell())
                    }
                    RD.get(idx_arv)!!.set(k, Matrix(0, 0)) // Empty matrix for arrival station

                    if (!RD.containsKey(idx_q)) {
                        RD.put(idx_q, MatrixCell())
                    }
                    RD.get(idx_q)!!.set(k, Matrix.concatColumns(F, X, null)) // [F, X] for queue station
                }
            }
        }
    } else {
        // Match MATLAB behavior: warn and return with no result for unsupported models
        line_warning(mfilename(object {}),
            "This model is not supported by SolverMAM yet. Returning with no result.")
    }

    return RD
}