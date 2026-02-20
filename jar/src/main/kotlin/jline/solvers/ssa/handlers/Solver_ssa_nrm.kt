package jline.solvers.ssa.handlers

import jline.lang.NetworkStruct
import jline.GlobalConstants
import jline.GlobalConstants.Inf
import jline.GlobalConstants.NegInf
import jline.lang.constant.NodeType
import jline.lang.constant.SchedStrategy
import jline.VerboseLevel
import jline.lang.state.State
import jline.lang.state.ToMarginal.toMarginalAggr
import jline.solvers.SolverOptions
import jline.util.Maths
import jline.util.RandomManager
import jline.util.matrix.Matrix
import kotlin.math.min

data class SolverSSAResultNRM(val QN: Matrix,
                              val UN: Matrix,
                              val RN: Matrix,
                              val TN: Matrix,
                              val CN: Matrix,
                              val XN: Matrix,
                              val sn: NetworkStruct)

fun solver_ssa_nrm(sn_in: NetworkStruct, options: SolverOptions): SolverSSAResultNRM {

    // Set master seed for reproducible SSA simulation
    RandomManager.setMasterSeed(options.seed)

    // ---------------------------------------------------------------------
    // Parameters & shorthands
    // ---------------------------------------------------------------------
    val sn: NetworkStruct = sn_in.copy()
    val samples = options.samples
    val R = sn.nclasses.toInt()
    val I = sn.nnodes.toInt()
    val state = sn.state

    // ---------------------------------------------------------------------
    // Stoichiometry & reaction mapping (self-loops included)
    // ---------------------------------------------------------------------
    var S = Matrix(0, I * R)  // will transpose at the end
    val fromIdx = mutableListOf<Int>()
    val toIdx = mutableListOf<MutableList<Int>>()
    val fromIR = mutableListOf<IntArray>()
    val probIR = mutableListOf<MutableList<Double>>()

    // this is currently M^2*R^2, it can be lowered to M*R decoupling the routing
    var k = 0
    for (ind in 0 until I) {
        for (r in 0 until R) {
            k++
            fromIR.add(intArrayOf(ind, r))
            fromIdx.add(ind * R + r)
            probIR.add(mutableListOf())
            toIdx.add(mutableListOf())
            val Srow = DoubleArray(I * R) // build stoichiometry row
            if (sn.isslc.get(r) != 0.0) {
                Srow[fromIdx[k - 1]] = NegInf
            } else {
                Srow[fromIdx[k - 1]] = -1.0
                for (jnd in 0 until I) {
                    for (s in 0 until R) {
                        if (sn.rtnodes.get(ind * R + r, jnd * R + s) > 0) {
                            toIdx[k - 1].add(jnd * R + s)
                            val p = sn.rtnodes.get(ind * R + r, jnd * R + s)
                            probIR[k - 1].add(p)
                            Srow[jnd * R + s] = Srow[jnd * R + s] + p
                        }
                    }
                }
            }
            // Expand S matrix if needed
            if (k > S.numRows) {
                val newS = Matrix(k, I * R)
                for (i in 0 until S.numRows) {
                    for (j in 0 until S.numCols) {
                        newS.set(i, j, S.get(i, j))
                    }
                }
                S = newS
            }
            // Set the row
            for (j in Srow.indices) {
                S.set(k - 1, j, Srow[j])
            }
        }
    }
    S = S.transpose() // states × reactions

    // ---------------------------------------------------------------------
    // Initial state vector
    // ---------------------------------------------------------------------
    val nvec0 = Matrix(I * R, 1) // initial state (aggregate state)
    for (ind in 0 until I) {
        if (sn.isstateful.get(ind) != 0.0) {
            val aggr = toMarginalAggr(sn,
                ind,
                state[sn.stateful.get(sn.nodeToStateful.get(ind).toInt())],
                null,
                null,
                null,
                null,
                null)
            for (r in 0 until R) {
                var nir = aggr.nir[r]
                if (nir.isInfinite()) {
                    if (sn.nodetype.get(ind) == NodeType.Source) {
                        nir = 1.0
                    } else {
                        error("Infinite population error.")
                    }
                }
                nvec0.set(ind * R + r, 0, nir)
            }
        }
    }

    // ---------------------------------------------------------------------
    // Rates and servers setup
    // ---------------------------------------------------------------------
    val mi = DoubleArray(I)
    val rates = Array(I) { DoubleArray(R) }
    for (ind in 0 until I) {
        if (sn.isstation.get(ind) != 0.0) {
            for (r in 0 until R) {
                val ist = sn.nodeToStation.get(ind).toInt()
                val muir = sn.rates.get(ist, r)
                if (!muir.isNaN()) {
                    rates[ind][r] = muir
                }
                mi[ind] = sn.nservers.get(ist)
            }
        } else {
            for (r in 0 until R) {
                rates[ind][r] = GlobalConstants.Immediate
                mi[ind] = GlobalConstants.MaxInt.toDouble()
            }
        }
        if (mi[ind].isInfinite()) {
            mi[ind] = GlobalConstants.MaxInt.toDouble()
        }
    }

    // ---------------------------------------------------------------------
    // Propensity functions
    // ---------------------------------------------------------------------
    val epstol = GlobalConstants.Zero
    val a = Array<(Matrix) -> Double>(fromIdx.size) { j ->
        { X: Matrix ->
            val ind = fromIR[j][0]
            val r = fromIR[j][1]
            val base = rates[ind][r]
            when {
                sn.isstation.get(ind) == 0.0 -> base * min(1.0, X.get(fromIdx[j], 0))
                else -> {
                    when (sn.sched.get(sn.stations.get(sn.nodeToStation.get(ind).toInt()))) {
                        SchedStrategy.EXT -> base
                        SchedStrategy.INF -> base * X.get(fromIdx[j], 0)
                        SchedStrategy.PS, SchedStrategy.DPS, SchedStrategy.GPS -> {
                            if (R == 1) { // single class
                                base * min(mi[ind], X.get(fromIdx[j], 0))
                            } else {
                                val total = (0 until R).sumOf { X.get(ind * R + it, 0) } + epstol
                                base * (X.get(fromIdx[j], 0) / total) * min(mi[ind], total)
                            }
                        }
                        SchedStrategy.PSPRIO, SchedStrategy.DPSPRIO, SchedStrategy.GPSPRIO -> {
                            if (R == 1) { // single class - no priority effect
                                base * min(mi[ind], X.get(fromIdx[j], 0))
                            } else {
                                val total = (0 until R).sumOf { X.get(ind * R + it, 0) } + epstol
                                // Find minimum priority among present classes
                                var minPrio = Int.MAX_VALUE
                                for (rr in 0 until R) {
                                    if (X.get(ind * R + rr, 0) > 0) {
                                        val classPrio = sn.classprio.get(rr).toInt()
                                        if (classPrio < minPrio) {
                                            minPrio = classPrio
                                        }
                                    }
                                }
                                val classPrio = sn.classprio.get(r).toInt()
                                // If total <= servers OR this class has highest priority
                                if (total <= mi[ind] || classPrio == minPrio) {
                                    base * (X.get(fromIdx[j], 0) / total) * min(mi[ind], total)
                                } else {
                                    0.0 // Not in highest priority group - no service
                                }
                            }
                        }
                        else -> base * X.get(fromIdx[j], 0)
                    }
                }
            }
        }
    }

    // ---------------------------------------------------------------------
    // Propensity function dependencies
    // ---------------------------------------------------------------------
    val D = Array(S.numCols) { mutableListOf<Int>() }
    for (k in 0 until D.size) {
        val J = mutableListOf<Int>() // set of state variables affected by reaction k
        for (i in 0 until S.numRows) {
            if (S.get(i, k) != 0.0) {
                J.add(i)
            }
        }
        val vecd = mutableListOf<Int>()
        for (j in J.indices) {
            // (ind-1)*R + r -> ind*R + r in 0-based indexing
            val pos = J[j]
            val r = pos % R
            val ind = (pos - r) / R
            for (rr in 0 until R) {
                vecd.add(ind * R + rr)
            }
        }
        // vecd now contains all state variables affected by the firing of
        // reaction k. We now find the propensity functions that depend
        // on those variables
        if (vecd.isNotEmpty()) {
            val vecdUnique = vecd.distinct()
            val vecs = mutableListOf<Int>()
            for (j in vecdUnique.indices) {
                for (kk in 0 until S.numCols) {
                    if (S.get(vecdUnique[j], kk) < 0.0) {
                        vecs.add(kk)
                    }
                }
            }
            D[k].addAll(vecs.distinct())
        }
    }

    // Having accounted for them in D, we can now remove self-loops markings
    for (i in 0 until S.numRows) {
        for (j in 0 until S.numCols) {
            if (S.get(i, j).isInfinite()) {
                S.set(i, j, 0.0)
            }
        }
    }

    // ---------------------------------------------------------------------
    // Run SSA/NRM and compute metrics directly
    // ---------------------------------------------------------------------
    val M = sn.nstations
    val K = sn.nclasses
    sn.njobs.transpose()

    // Initialize metric matrices
    val QN = Matrix(M, K).apply { fill(0.0) }
    val UN = Matrix(M, K).apply { fill(0.0) }
    val RN = Matrix(M, K).apply { fill(0.0) }
    val TN = Matrix(M, K).apply { fill(0.0) }
    val CN = Matrix(1, K).apply { fill(0.0) }
    val XN = Matrix(1, K).apply { fill(0.0) }

    // Run simulation with direct metric computation
    next_reaction_method_direct(S, D, a, nvec0, samples, options, QN, UN, RN, TN, CN, XN, sn, fromIdx, fromIR)

    return SolverSSAResultNRM(QN, UN, RN, TN, CN, XN, sn)
}

// ======================================================================
// Next-Reaction Method with direct metric computation
// ======================================================================
fun next_reaction_method_direct(S: Matrix,
                                D: Array<MutableList<Int>>,
                                a: Array<(Matrix) -> Double>,
                                nvec0: Matrix,
                                samples: Int,
                                options: SolverOptions,
                                QN: Matrix,
                                UN: Matrix,
                                RN: Matrix,
                                TN: Matrix,
                                CN: Matrix,
                                XN: Matrix,
                                sn: NetworkStruct,
                                fromIdx: List<Int>,
                                fromIR: List<IntArray>) {
    val numReactions = S.numCols
    val Ak = DoubleArray(numReactions) { a[it](nvec0) }
    val Pk = DoubleArray(numReactions) { -kotlin.math.ln(Maths.rand()) }
    val Tk = DoubleArray(numReactions)
    val times = mutableListOf(0.0)
    val states = mutableListOf<DoubleArray>()
    var nvec = nvec0.copy()

    val M = sn.nstations
    val K = sn.nclasses
    val R = sn.nclasses
    sn.nnodes
    val NK = sn.njobs.transpose()

    // Accumulate metrics during simulation
    var totalTime = 0.0
    val PH = sn.proc
    val rates = sn.rates
    val S_servers = sn.nservers

    var n = 0
    while (n < samples) {
        val tau = DoubleArray(numReactions) { i ->
            if (Ak[i] > 0) (Pk[i] - Tk[i]) / Ak[i] else Inf
        }
        val kfire = tau.withIndex().minByOrNull { it.value }!!.index
        val dt = tau[kfire]
        totalTime += dt

        // Accumulate state-dependent metrics during this time interval
        for (ist in 0 until M) {
            val ind = sn.stationToNode.get(ist).toInt()
            for (k in 0 until K) {
                val currentPop = nvec.get(ind * R + k, 0)

                // Accumulate queue length (QN)
                QN.set(ist, k, QN.get(ist, k) + currentPop * dt)

                // Compute throughput contribution (TN) from departures
                val irIndex = ind * R + k
                val depRate = fromIdx.withIndex().find { it.value == irIndex }?.let { Ak[it.index] } ?: 0.0
                TN.set(ist, k, TN.get(ist, k) + depRate * dt)

                // Compute utilization based on scheduling policy
                when (sn.sched.get(sn.stations.get(ist))) {
                    SchedStrategy.INF, SchedStrategy.EXT -> {
                        UN.set(ist, k, UN.get(ist, k) + currentPop * dt)
                    }
                    SchedStrategy.PS, SchedStrategy.DPS, SchedStrategy.GPS -> {
                        val phEntry = PH.get(sn.stations.get(ist))?.get(sn.jobclasses.get(k))
                        if (phEntry != null && !phEntry.isEmpty) {
                            val servers = S_servers.get(ist)
                            val totalPop = (0 until R).sumOf { nvec.get(ind * R + it, 0) }
                            val utilization = if (totalPop > 0) {
                                (currentPop / totalPop) * min(servers, totalPop) / servers
                            } else 0.0
                            UN.set(ist, k, UN.get(ist, k) + utilization * dt)
                        }
                    }
                    SchedStrategy.PSPRIO, SchedStrategy.DPSPRIO, SchedStrategy.GPSPRIO -> {
                        val phEntry = PH.get(sn.stations.get(ist))?.get(sn.jobclasses.get(k))
                        if (phEntry != null && !phEntry.isEmpty) {
                            val servers = S_servers.get(ist)
                            val totalPop = (0 until R).sumOf { nvec.get(ind * R + it, 0) }
                            if (totalPop > 0) {
                                // Find minimum priority among present classes
                                var minPrio = Int.MAX_VALUE
                                for (rr in 0 until R) {
                                    if (nvec.get(ind * R + rr, 0) > 0) {
                                        val classPrio = sn.classprio.get(rr).toInt()
                                        if (classPrio < minPrio) {
                                            minPrio = classPrio
                                        }
                                    }
                                }
                                val classPrio = sn.classprio.get(k).toInt()
                                val utilization = if (totalPop <= servers || classPrio == minPrio) {
                                    // If total <= servers, all jobs get service proportionally
                                    // Or this class has highest priority
                                    (currentPop / totalPop) * min(servers, totalPop) / servers
                                } else {
                                    0.0 // Not in highest priority group - no utilization
                                }
                                UN.set(ist, k, UN.get(ist, k) + utilization * dt)
                            }
                        }
                    }
                    else -> {}
                }
            }
        }

        times += times.last() + dt

        // Update state
        nvec.addEq(S.getColumn(kfire))

        // Update times and propensities
        for (i in 0 until numReactions) Tk[i] += Ak[i] * dt
        D[kfire].forEach { i -> Ak[i] = a[i](nvec) }
        Pk[kfire] -= kotlin.math.ln(Maths.rand())

        states += DoubleArray(nvec.numRows) { r -> nvec.get(r, 0) }
        n++
        printProgress(options, n)
    }

    // Normalize metrics by total time
    if (totalTime > 0) {
        for (ist in 0 until M) {
            for (k in 0 until K) {
                QN.set(ist, k, QN.get(ist, k) / totalTime)
                UN.set(ist, k, UN.get(ist, k) / totalTime)
                TN.set(ist, k, TN.get(ist, k) / totalTime)
            }
        }
    }

    // Compute derived metrics
    for (k in 0 until K) {
        // Throughput at reference station
        sn.stationToNode.get(sn.refstat.get(k).toInt()).toInt()
        XN.set(0, k, TN.get(sn.refstat.get(k).toInt(), k))

        // Response times and cycle times
        for (ist in 0 until M) {
            if (TN.get(ist, k) > 0) {
                RN.set(ist, k, QN.get(ist, k) / TN.get(ist, k))
            } else {
                RN.set(ist, k, 0.0)
            }
        }

        if (XN.get(0, k) > 0) {
            CN.set(0, k, NK.get(k) / XN.get(0, k))
        }
    }

    // Handle NaN values
    QN.apply(Double.NaN, 0.0, "equal")
    UN.apply(Double.NaN, 0.0, "equal")
    RN.apply(Double.NaN, 0.0, "equal")
    XN.apply(Double.NaN, 0.0, "equal")
    TN.apply(Double.NaN, 0.0, "equal")
    CN.apply(Double.NaN, 0.0, "equal")
}

private fun printProgress(options: SolverOptions, samples_collected: Int) {
    if (options.method != "parallel" && (options.verbose == VerboseLevel.STD || options.verbose == VerboseLevel.DEBUG)) {
        if (samples_collected == 2) {
            System.out.printf("SSA samples: %9d ", samples_collected)
            System.out.flush()
        } else if (options.verbose == VerboseLevel.DEBUG) {
            System.out.printf("\b\b\b\b\b\b\b\b\b\b %9d", samples_collected)
            System.out.flush()
        } else if (samples_collected % 1000 == 0) {
            System.out.printf("\b\b\b\b\b\b\b\b\b\b %9d", samples_collected)
            System.out.flush()
        }
        if (samples_collected == options.samples) {
            println()
        }
    }
}