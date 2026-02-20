package jline.solvers.ssa.handlers

import jline.lang.Event
import jline.lang.NetworkStruct
import jline.lang.constant.EventType
import jline.lang.constant.NodeType
import jline.VerboseLevel
import jline.lang.nodes.StatefulNode
import jline.lang.state.AfterGlobalEvent
import jline.lang.state.EventCache
import jline.lang.state.State
import jline.lang.state.ToMarginal.toMarginal
import jline.solvers.SolverOptions
import jline.solvers.ssa.SSAValues
import jline.solvers.ssa.SolverSSA
import jline.util.Maths
import jline.util.RandomManager
import jline.util.Utils
import jline.util.matrix.Matrix
import jline.streaming.Collector
import org.apache.commons.math3.util.FastMath
import org.ejml.data.DMatrixRMaj
import org.ejml.data.DMatrixSparseCSC
import org.ejml.data.DMatrixSparseTriplet
import org.ejml.ops.DConvertMatrixStruct
import java.util.stream.Collectors

fun solver_ssa(sn_in: NetworkStruct,
               eventCache: EventCache,
               init_state: MutableMap<StatefulNode?, Matrix?>,
               options: SolverOptions,
               solverSSA: SolverSSA): SSAValues {
    val options = solverSSA.getOptions()
    val sn: NetworkStruct = sn_in.copy()

    // Set master seed for reproducible SSA simulation
    RandomManager.setMasterSeed(options.seed)


    // generate local state spaces
    val nstateful = sn.nstateful
    val R = sn.nclasses
    val N = sn.njobs.transpose()
    val sync = sn.sync
    val csmask = sn.csmask

    // Get properly dimensioned cutoff matrix
    val cutoff = options.getCutoffMatrix(sn.nstations, sn.nclasses)

    val Np = N.transpose()
    val capacityc = Matrix(sn.nnodes, sn.nclasses)
    capacityc.zero()

    for (ind in 0..<sn.nnodes) {
        if (sn.isstation.get(ind) == 1.0) {
            val ist = sn.nodeToStation.get(ind).toInt()
            for (r in 0..<sn.nclasses) {
                var c = 0
                for (i in 0..<sn.chains.numRows) {
                    if (sn.chains.get(i, r) == 1.0) {
                        c = i
                    }
                }

                val proc_m = sn.proc.get(sn.stations.get(ist))!!.get(sn.jobclasses.get(r))!!.get(0)
                var disabled = false
                for (row in 0..<proc_m.numRows) {
                    for (col in 0..<proc_m.numCols) {
                        if (java.lang.Double.isNaN(proc_m.get(row, col)) && sn.nodetype.get(ind) != NodeType.Place) {
                            disabled = true
                        }
                    }
                }

                if (!sn.visits.get(c)!!.isEmpty && sn.visits.get(c)!!.get(ist, r) == 0.0) {
                    capacityc.set(ind, r, 0)
                } else if (!sn.proc.isEmpty() && !sn.proc.get(sn.stations.get(ist))!!.get(sn.jobclasses.get(r))!!
                        .isEmpty && disabled) {
                    capacityc.set(ind, r, 0)
                } else {
                    if (Utils.isInf(N.get(r))) {
                        capacityc.set(ind, r, Maths.min(cutoff.get(ist, r), sn.classcap.get(ist, r)))
                    } else {
                        // sum values in sn,njobs at indices where the c-th row of sn.chains is true
                        var njobs_sum = 0
                        for (i in 0..<sn.njobs.numCols) {
                            if (sn.chains.get(c, i) == 1.0) {
                                njobs_sum += sn.njobs.get(i).toInt()
                            }
                        }
                        capacityc.set(ind, r, njobs_sum)
                    }
                }
            }
            val capacity_sum = capacityc.sumRows().get(ind).toInt()
            if (Utils.isInf(sn.nservers.get(ist))) {
                sn.nservers.set(ist, capacity_sum.toDouble())
            }
            for (col in 0..<sn.cap.numCols) {
                sn.cap.set(ist, col, capacity_sum)
            }
            for (col in 0..<sn.classcap.numCols) {
                sn.classcap.set(ist, col, capacityc.get(ind, col))
            }
        }
    }

    if (Np.hasInfinite()) {
        // set all elements of Np where they are Infinity to 0
        for (col in 0..<Np.numCols) {
            if (Utils.isInf(Np.get(0, col))) {
                Np.set(0, col, 0)
            }
        }
    }

    val init_state_hashed = Matrix(1, nstateful)
    init_state_hashed.zero() // pick first state in space{i}

    val arvRatesSamples: MutableMap<Int?, Matrix?> = HashMap<Int?, Matrix?>()
    val depRatesSamples: MutableMap<Int?, Matrix?> = HashMap<Int?, Matrix?>()
    // Swap map entries with rows of contained matrices compared with matlab
    for (r in 0..<options.samples) {
        val m = Matrix(R, nstateful)
        m.zero()
        arvRatesSamples.put(r, m)
        depRatesSamples.put(r, m.copy())
    }
    val A = sync.size
    var samples_collected = 1
    var state = init_state_hashed.copy()
    var cur_state: MutableMap<Int?, Matrix> = HashMap<Int?, Matrix>()
    val nir: MutableMap<Int?, Matrix> = HashMap<Int?, Matrix>()

    for (ind in 0..<sn.nnodes) {
        if (sn.isstateful.get(ind) == 1.0) {
            val isf = sn.nodeToStateful.get(ind).toInt()
            val state_space = Matrix.extractRows(init_state.get(solverSSA.model.getStatefulNodes().get(isf)),
                state.get(isf).toInt(),
                state.get(isf).toInt() + 1,
                null)
            cur_state.put(isf, state_space)

            if (sn.isstation.get(ind) == 1.0) {
                val ist = sn.nodeToStation.get(ind).toInt()
                // transpose to get nir as column i
                nir.put(ist, toMarginal(sn, ind, state_space, null, null, null, null, null).nir.transpose())
            }
        }
    }


    state = Matrix(0, 0)
    val statelen = Matrix(cur_state.size, 1)

    for (ind in 0..<cur_state.size) {
        if (cur_state.containsKey(ind)) {
            val row: Matrix = cur_state.get(ind)!!
            if (state.isEmpty) {
                state = row
            } else {
                state = Matrix.concatColumns(state, row, null)
            }
            statelen.set(ind, row.numElements.toDouble())
        } else {
            statelen.set(ind, 0.0)
        }
    }

    val nSamples = options.samples - 1

    // Use dense DoubleArray for tranSync (avoids sparse Matrix.set() O(nnz) per call)
    val tranSyncData = DoubleArray(nSamples)

    // Compute initial tranState dimensions: 1 (time) + state length
    val z = Matrix(1, 1)
    z.zero()
    var tranStateInit = Matrix.concatColumns(z, state, null).transpose()
    var tranStateRows = tranStateInit.numRows  // state_size + 1

    // Use dense DMatrixRMaj for tranState: rows = state_size+1, cols = nSamples
    // This avoids O(nnz) per element with sparse DMatrixSparseCSC
    var tranStateD = DMatrixRMaj(tranStateRows, nSamples)
    // Copy initial state to first column
    for (row in 0..<tranStateRows) {
        tranStateD.set(row, 0, tranStateInit.get(row, 0))
    }

    samples_collected = 1

    // Compute SSq dimensions
    var SSqRows = 0
    for (ist in 0..<sn.nstations) {
        if (nir.containsKey(ist)) {
            SSqRows += nir.get(ist)!!.numRows
        }
    }

    // Use dense DMatrixRMaj for SSq
    val SSqD = DMatrixRMaj(SSqRows, nSamples)
    // Copy initial nir to first column
    var rowOff = 0
    for (ist in 0..<sn.nstations) {
        if (nir.containsKey(ist)) {
            val col = nir.get(ist)!!
            for (i in 0..<col.numRows) {
                SSqD.set(rowOff + i, 0, col.get(i))
            }
            rowOff += col.numRows
        }
    }

    val local = sn.nnodes

    val node_a: MutableMap<Int?, Int?> = HashMap<Int?, Int?>()
    val node_p: MutableMap<Int?, Int?> = HashMap<Int?, Int?>()
    val class_a: MutableMap<Int?, Int?> = HashMap<Int?, Int?>()
    val class_p: MutableMap<Int?, Int?> = HashMap<Int?, Int?>()
    val event_a: MutableMap<Int?, EventType?> = HashMap<Int?, EventType?>()
    val event_p: MutableMap<Int?, EventType?> = HashMap<Int?, EventType?>()
    val outprob_a: MutableMap<Int?, Double?> = HashMap<Int?, Double?>()
    val outprob_p: MutableMap<Int?, Double?> = HashMap<Int?, Double?>()
    for (act in 0..<A) {
        val active: Event = sync.get(act)!!.active.get(0)!!
        val passive: Event = sync.get(act)!!.passive.get(0)!!
        node_a.put(act, active.getNode())
        node_p.put(act, passive.getNode())
        class_a.put(act, active.jobClass)
        class_p.put(act, passive.jobClass)
        event_a.put(act, active.getEvent())
        event_p.put(act, passive.getEvent())
    }

    // next state if a give action occurs
    val next_state: MutableMap<Int?, MutableMap<Int?, Matrix>> = HashMap<Int?, MutableMap<Int?, Matrix>>()
    val isSimulation = true // allow state vector to grow, e.g., for FCFS buffers
    var cur_time = 0.0
    val enabled_rates: MutableMap<Int?, Double?> = HashMap<Int?, Double?>()
    val enabled_sync: MutableMap<Int?, Int?> = HashMap<Int?, Int?>()
    var cur_state_1: MutableMap<Int?, Matrix?> = HashMap<Int?, Matrix?>()
    while (samples_collected < options.samples && cur_time <= options.timespan[1]) {
        if (samples_collected == 1) {
            if (options.method == "parallel") {
                //System.out.printf("SSA samples: %6d\n", samples_collected);
                //System.out.println("Running SSA with " + this.DEFAULT_THREADS + " threads\n");
            }
        }

        val node_a_sf: MutableMap<Int?, Int?> = HashMap<Int?, Int?>()
        val node_p_sf: MutableMap<Int?, Int?> = HashMap<Int?, Int?>()
        val prob_sync_p: MutableMap<Int?, Double?> = HashMap<Int?, Double?>()
        enabled_rates.clear()
        enabled_sync.clear()

        // Find enabled pairwise synchronisations
        solver_ssa_findenabled(sn,
            eventCache,
            A,
            node_a,
            next_state,
            cur_state,
            event_a,
            class_a,
            isSimulation,
            outprob_a,
            node_p,
            local,
            event_p,
            class_p,
            outprob_p,
            prob_sync_p,
            sync,
            node_a_sf,
            node_p_sf,
            depRatesSamples,
            samples_collected,
            arvRatesSamples,
            csmask,
            enabled_rates,
            enabled_sync,
            solverSSA)

        // Handle global synchronizations (for SPN Transitions)
        val gsync = sn.gsync
        if (gsync != null && gsync.isNotEmpty()) {
            val G = gsync.size
            var ctr = enabled_rates.size  // Continue from where regular syncs left off

            for (gact in 0 until G) {
                val gSync = gsync[gact] ?: continue
                if (gSync.active.isEmpty()) continue

                val gind = gSync.active[0].node  // node index for global event

                // Convert cur_state map to List<Matrix> for afterGlobalEvent
                val glspace = ArrayList<Matrix>()
                for (isf in 0 until sn.nstateful) {
                    glspace.add(cur_state[isf]?.copy() ?: Matrix(1, 1))
                }

                val result = AfterGlobalEvent.afterGlobalEvent(sn, gind, glspace, gSync, isSimulation)
                val outrate = result.outrate
                val outprob = result.outprob
                val outglspace = result.outglspace

                // Find non-zero rate*prob combinations
                for (ia in 0 until outrate.numRows) {
                    val rate = outrate.get(ia, 0)
                    val prob = if (outprob.numRows > ia) outprob.get(ia, 0) else 1.0

                    if (!rate.isNaN() && rate > 0 && !prob.isNaN() && prob > 0) {
                        val combinedRate = rate * prob
                        enabled_rates[ctr] = combinedRate
                        enabled_sync[ctr] = A + gact  // Global syncs indexed after regular syncs

                        // Store the next state for this global sync
                        val nextStateForGsync = HashMap<Int?, Matrix>()
                        for (isf in 0 until sn.nstateful) {
                            if (outglspace.size > isf && outglspace[isf] != null) {
                                nextStateForGsync[isf] = outglspace[isf].copy()
                            } else if (cur_state.containsKey(isf)) {
                                nextStateForGsync[isf] = cur_state[isf]!!.copy()
                            }
                        }
                        next_state[A + gact] = nextStateForGsync

                        // Record departure/arrival rates for FIRE events
                        if (gSync.active[0].event == EventType.FIRE) {
                            for (pev in gSync.passive) {
                                val pevNode = pev.node
                                if (pevNode >= 0 && pevNode < sn.nnodes &&
                                    !sn.nodeToStateful[pevNode].isNaN() &&
                                    sn.nodeToStateful[pevNode] >= 0) {
                                    val pevIsf = sn.nodeToStateful[pevNode].toInt()

                                    // For class, we use mode as a proxy since passive events
                                    // affect all classes in that mode
                                    val pevClass = pev.mode % R

                                    if (pev.event == EventType.PRE) {
                                        // Departure from input Place (consuming tokens)
                                        val depMatrix = depRatesSamples[samples_collected]
                                        if (depMatrix != null) {
                                            val currentVal = depMatrix.get(pevClass, pevIsf)
                                            depMatrix.set(pevClass, pevIsf, currentVal + combinedRate)
                                        }
                                    } else if (pev.event == EventType.POST) {
                                        // Arrival at output Place (producing tokens)
                                        val arvMatrix = arvRatesSamples[samples_collected]
                                        if (arvMatrix != null) {
                                            val currentVal = arvMatrix.get(pevClass, pevIsf)
                                            arvMatrix.set(pevClass, pevIsf, currentVal + combinedRate)
                                        }
                                    }
                                }
                            }
                        }

                        ctr++
                    }
                }
            }
        }

        val enabled_rates_m = Matrix(1, enabled_rates.size)
        for (i in 0..<enabled_rates.size) {
            enabled_rates_m.set(0, i, enabled_rates.get(i)!!)
        }
        val tot_rate = enabled_rates_m.elementSum()
        val cum_sum = enabled_rates_m.cumsumViaRow()
        val cum_rate = Matrix.scaleMult(cum_sum, 1.0 / tot_rate)

        val rand = Maths.rand()
        var firing_ctr = -1
        for (i in 0..<cum_rate.numElements) {
            if (rand > cum_rate.get(i)) {
                firing_ctr = i
            } else {
                break
            }
        }
        firing_ctr++
        if (enabled_sync.isEmpty()) {
            throw RuntimeException("SSA simulation entered a deadlock before collecting all samples, " + "no synchronization is enabled.")
        }

        // This part is needed to ensure that when the state vector grows the
        // padding of zeros, this is done on the left (e.g., for FCFS buffers)
        tranStateD = update_paddings_dense(sn, cur_state, statelen, tranStateD)
        tranStateRows = tranStateD.numRows

        // Apply the time increment
        val dt = -(FastMath.log(Maths.rand()) / tot_rate)
        cur_time += dt

        // Update sample data using dense arrays - O(1) per element
        save_log_dense(dt,
            cur_state,
            tranStateD,
            samples_collected,
            tranSyncData,
            enabled_sync,
            firing_ctr,
            sn,
            nir,
            cur_state,
            SSqD,
            solverSSA.getStreamingCollector(),
            cur_time)

        cur_state = next_state.get(enabled_sync.get(firing_ctr))!!

        // Display progress to the user
        samples_collected++
        print_progress(options, samples_collected)
    }

    // Copy final state ONCE after loop
    cur_state_1 = cur_state.entries.stream().collect(Collectors.toMap({ it.key },
        { it.value!!.copy() }
    ))

    // --- Post-processing using dense DMatrixRMaj for O(1) element access ---

    // Cumulative sum of dwell times using DMatrixRMaj
    val timesCumSumD = DMatrixRMaj(nSamples, 1)
    var cumSum = 0.0
    for (i in 0 until nSamples) {
        cumSum += tranStateD.get(0, i)
        timesCumSumD.set(i, 0, cumSum)
    }
    val timesCumSum = denseToSparseMatrix(timesCumSumD)

    // Find unique state rows (columns of tranStateD, excluding row 0 which is time)
    val rowKeyMap = HashMap<String, MutableList<Int>>()
    val uniqueKeysOrdered = ArrayList<String>()
    val estimatedCapacity = (tranStateRows - 1) * 20

    for (col in 0 until nSamples) {
        val sb = StringBuilder(estimatedCapacity)
        for (row in 1 until tranStateRows) {
            if (row > 1) sb.append(',')
            sb.append(tranStateD.get(row, col))
        }
        val key = sb.toString()
        val indices = rowKeyMap[key]
        if (indices == null) {
            val newList = ArrayList<Int>()
            newList.add(col)
            rowKeyMap[key] = newList
            uniqueKeysOrdered.add(key)
        } else {
            indices.add(col)
        }
    }

    // Sort unique keys lexicographically
    uniqueKeysOrdered.sort()

    val numUniqueStates = uniqueKeysOrdered.size
    val ui = IntArray(numUniqueStates)
    val uj = Array<List<Int>>(numUniqueStates) { emptyList() }

    for (i in 0 until numUniqueStates) {
        val key = uniqueKeysOrdered[i]
        val indices = rowKeyMap[key]!!
        ui[i] = indices[0]
        uj[i] = indices
    }

    // Build tranSysState
    val stateSizeCount = cur_state_1.size
    val statesz = IntArray(stateSizeCount)
    for (i in 0 until stateSizeCount) {
        statesz[i] = cur_state_1[i]!!.numElements
    }

    val tranSysState = HashMap<Int, Matrix?>()
    tranSysState[0] = timesCumSum

    var start_index = 1
    for (j in 0 until stateSizeCount) {
        val size = statesz[j]
        val end_index = start_index + size
        val tmpD = DMatrixRMaj(nSamples, size)
        for (samp in 0 until nSamples) {
            for (k in start_index until end_index) {
                tmpD.set(samp, k - start_index, tranStateD.get(k, samp))
            }
        }
        tranSysState[j + 1] = denseToSparseMatrix(tmpD)
        start_index = end_index
    }

    // Build arvRates and depRates using DMatrixRMaj, then wrap
    val arvRatesD = Array(R) { DMatrixRMaj(numUniqueStates, sn.nstateful) }
    val depRatesD = Array(R) { DMatrixRMaj(numUniqueStates, sn.nstateful) }

    // Compute pi using DMatrixRMaj
    val piD = DMatrixRMaj(1, numUniqueStates)
    for (s in 0 until numUniqueStates) {
        val stateIndexes = uj[s]
        var dwellSum = 0.0
        for (idx in stateIndexes) {
            dwellSum += tranStateD.get(0, idx)
        }
        piD.set(0, s, dwellSum)
    }

    // Build SSq from unique states using DMatrixRMaj
    val SSqNewD = DMatrixRMaj(numUniqueStates, SSqRows)
    for (s in 0 until numUniqueStates) {
        val srcCol = ui[s]
        for (row in 0 until SSqRows) {
            SSqNewD.set(s, row, SSqD.get(row, srcCol))
        }
    }
    var SSq = denseToSparseMatrix(SSqNewD)

    // Fill arvRates/depRates from samples into DMatrixRMaj
    for (ind in 0 until sn.nnodes) {
        if (sn.isstateful.get(ind) == 1.0) {
            val isf = sn.nodeToStateful.get(ind).toInt()
            for (s in 0 until numUniqueStates) {
                val uis = ui[s]
                for (r in 0 until R) {
                    arvRatesD[r].set(s, isf, arvRatesSamples[uis]!!.get(r, isf))
                    depRatesD[r].set(s, isf, depRatesSamples[uis]!!.get(r, isf))
                }
            }
        }
    }

    // Wrap in Matrix after all data is written
    val arvRates = HashMap<Int, Matrix?>()
    val depRates = HashMap<Int, Matrix?>()
    repeat(R) { i ->
        arvRates[i] = denseToSparseMatrix(arvRatesD[i])
        depRates[i] = denseToSparseMatrix(depRatesD[i])
    }

    // Normalize pi
    val piSum = piD.data.sum()
    for (i in 0 until numUniqueStates) {
        piD.set(0, i, piD.get(0, i) / piSum)
    }
    var pi = denseToSparseMatrix(piD)

    // Convert tranSync to Matrix using dense path
    val tranSyncRMaj = DMatrixRMaj(nSamples, 1)
    for (i in 0 until nSamples) {
        tranSyncRMaj.set(i, 0, tranSyncData[i])
    }
    val tranSync = denseToSparseMatrix(tranSyncRMaj)

    return SSAValues(pi, SSq, arvRates, depRates, tranSysState, tranSync, sn)
}

/**
 * Dense version of save_log using DMatrixRMaj for O(1) element access.
 */
fun save_log_dense(dt: Double,
                   cur_state: MutableMap<Int?, Matrix>,
                   tranStateD: DMatrixRMaj,
                   samples_collected: Int,
                   tranSyncData: DoubleArray,
                   enabled_sync: MutableMap<Int?, Int?>,
                   firing_ctr: Int,
                   sn: NetworkStruct,
                   nir: MutableMap<Int?, Matrix>,
                   stateCell: MutableMap<Int?, Matrix>,
                   SSqD: DMatrixRMaj,
                   streamingCollector: Collector? = null,
                   curTime: Double = 0.0) {
    // Set dt as first element, then copy state data directly to tranStateD - O(1) per element
    val colIdx = samples_collected - 1
    tranStateD.set(0, colIdx, dt)
    var offset = 1
    for (ind in 0..<sn.nnodes) {
        if (cur_state.containsKey(ind)) {
            val row: Matrix = cur_state.get(ind)!!
            for (i in 0..<row.numElements) {
                tranStateD.set(offset + i, colIdx, row.get(i))
            }
            offset += row.numElements
        }
    }

    // +1 to convert 0-based Java indexing to 1-based indexing for LINE compatibility
    tranSyncData[samples_collected - 1] = (enabled_sync.get(firing_ctr)!! + 1).toDouble()

    // Compute marginals for each station
    for (ind in 0..<sn.nnodes) {
        if (sn.isstation.get(ind) == 1.0) {
            val isf = sn.nodeToStateful.get(ind).toInt()
            val ist = sn.nodeToStation.get(ind).toInt()
            nir.put(ist, toMarginal(sn, ind, stateCell.get(isf), null, null, null, null, null).nir.transpose())
        }
    }

    // Copy nir values directly to SSqD - O(1) per element
    // IMPORTANT: iterate by station index to ensure correct ordering
    var rowOffset = 0
    for (ist in 0..<sn.nstations) {
        if (nir.containsKey(ist)) {
            val col: Matrix = nir.get(ist)!!
            for (i in 0..<col.numRows) {
                SSqD.set(rowOffset + i, colIdx, col.get(i))
            }
            rowOffset += col.numRows
        }
    }

    // Notify streaming collector if active
    if (streamingCollector != null) {
        // Build nir_col only when streaming collector is active (rare case)
        var nir_col = Matrix(0, 0)
        for (ist in 0..<sn.nstations) {
            if (nir.containsKey(ist)) {
                val col: Matrix = nir.get(ist)!!
                if (nir_col.isEmpty) {
                    nir_col = col
                } else {
                    nir_col = Matrix.concatRows(nir_col, col, null)
                }
            }
        }
        streamingCollector.recordState(curTime, dt, nir_col, null, null)
    }
}

fun print_progress(options: SolverOptions, samples_collected: Int) {
    if (options.method != "parallel" && (options.verbose == VerboseLevel.STD || options.verbose == VerboseLevel.DEBUG)) {
        if (samples_collected == 2) {
            System.out.printf("SSA samples: %6d ", samples_collected)
            System.out.flush()
        } else if (options.verbose == VerboseLevel.DEBUG) {
            System.out.printf("\b\b\b\b\b\b\b %6d", samples_collected)
            System.out.flush()
        } else if (samples_collected % 100 == 0) {
            System.out.printf("\b\b\b\b\b\b\b %6d", samples_collected)
            System.out.flush()
        }
        if (samples_collected == options.samples) {
            println()
        }
    }
}

/**
 * Convert a DMatrixRMaj to a sparse-backed Matrix via triplet format.
 * Triplet -> CSC conversion is O(nnz) and avoids the O(n^2) pathology
 * of direct DMatrixRMaj -> DMatrixSparseCSC conversion.
 */
private fun denseToSparseMatrix(d: DMatrixRMaj): Matrix {
    val rows = d.numRows
    val cols = d.numCols
    val data = d.data

    // Count non-zeros first for exact triplet allocation
    var nnz = 0
    for (i in 0 until data.size) {
        if (data[i] != 0.0) nnz++
    }

    val triplet = DMatrixSparseTriplet(rows, cols, nnz)
    // Add entries in row-major order (data is row-major in DMatrixRMaj)
    for (i in 0 until rows) {
        val rowOff = i * cols
        for (j in 0 until cols) {
            val v = data[rowOff + j]
            if (v != 0.0) {
                triplet.addItem(i, j, v)
            }
        }
    }

    val csc = DConvertMatrixStruct.convert(triplet, null as DMatrixSparseCSC?)
    return Matrix(csc as org.ejml.data.DMatrix)
}

/**
 * Dense version of update_paddings using DMatrixRMaj.
 * When state vectors grow (e.g., FCFS buffers), pad with zeros by creating a new larger matrix.
 */
fun update_paddings_dense(sn: NetworkStruct,
                          stateCell: MutableMap<Int?, Matrix>,
                          statelen: Matrix,
                          tranStateD: DMatrixRMaj): DMatrixRMaj {
    var result = tranStateD
    for (ind in 0..<sn.nnodes) {
        if (sn.isstation.get(ind) == 1.0) {
            val isf = sn.nodeToStateful.get(ind).toInt()
            val deltalen = (stateCell.get(isf)!!.numElements > statelen.get(isf))
            if (deltalen) {
                statelen.set(isf, stateCell.get(isf)!!.numElements.toDouble())
                // padding: insert a row of zeros at position shift+1
                var shift = 0
                if (ind > 0) {
                    for (col in 0..<isf) {
                        shift += statelen.get(col).toInt()
                    }
                }
                val oldRows = result.numRows
                val nCols = result.numCols
                val newResult = DMatrixRMaj(oldRows + 1, nCols)
                // Copy rows 0..shift (inclusive)
                for (r in 0..shift) {
                    for (c in 0 until nCols) {
                        newResult.set(r, c, result.get(r, c))
                    }
                }
                // Row shift+1 is zeros (default in DMatrixRMaj)
                // Copy rows shift+1..oldRows-1 to shift+2..oldRows
                for (r in shift + 1 until oldRows) {
                    for (c in 0 until nCols) {
                        newResult.set(r + 1, c, result.get(r, c))
                    }
                }
                result = newResult
            }
        }
    }
    return result
}
