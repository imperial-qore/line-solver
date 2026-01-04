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

    val tranSync = Matrix(options.samples - 1, 1)
    tranSync.zero()
    val z = Matrix(1, 1)
    z.zero()
    var tranState = Matrix.concatColumns(z, state, null).transpose()

    val tranState_s = Matrix(tranState.numRows, options.samples - 1)
    for (row in 0..<tranState.numRows) {
        for (col in 0..<tranState.numCols) {
            tranState_s.set(row, col, tranState.get(row, col))
        }
    }
    tranState = tranState_s.copy()

    samples_collected = 1

    var SSq = Matrix(0, 0)
    // Iterate by station index to ensure correct ordering
    for (ist in 0..<sn.nstations) {
        if (nir.containsKey(ist)) {
            val col: Matrix = nir.get(ist)!!
            if (SSq.isEmpty) {
                SSq = col
            } else {
                SSq = Matrix.concatRows(SSq, col, null)
            }
        }
    }

    val SSq_s = Matrix(SSq.numRows, options.samples - 1)
    // copy all o SSq into SSq_new
    for (row in 0..<SSq.numRows) {
        for (col in 0..<SSq.numCols) {
            SSq_s.set(row, col, SSq.get(row, col))
        }
    }
    SSq = SSq_s.copy()

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
        tranState = update_paddings(sn, cur_state, statelen, tranState)

        // Apply the time increment
        val dt = -(FastMath.log(Maths.rand()) / tot_rate)
        cur_time += dt

        // Update sample data
        save_log(dt,
            cur_state,
            tranState,
            samples_collected,
            tranSync,
            enabled_sync,
            firing_ctr,
            sn,
            nir,
            cur_state,
            SSq,
            solverSSA.getStreamingCollector(),
            cur_time)

        cur_state = next_state.get(enabled_sync.get(firing_ctr))!!

        // Display progress to the user
        samples_collected++
        print_progress(options, samples_collected)
    }

    // Copy final state ONCE after loop - was incorrectly inside loop causing O(n²) overhead
    cur_state_1 = cur_state.entries.stream().collect(Collectors.toMap({ it.key },
        { it.value!!.copy() }
    ))

    tranState = tranState.transpose()
    val tranStateTimes = Matrix.extractColumn(tranState, 0, null)
    val timesCumSum = tranStateTimes.cumsumViaCol()

    val uniqueRows = Matrix.uniqueRowIndexesFromColumn(tranState, 1)
    val ui = uniqueRows.vi
    val uj = uniqueRows.vj

    val stateSizeCount = cur_state_1.size
    val statesz = Matrix(1, stateSizeCount)
    for (i in 0 until stateSizeCount) {
        statesz.set(i, cur_state_1[i]!!.numElements.toDouble())
    }

    val tranSysState = HashMap<Int, Matrix?>()
    tranSysState[0] = timesCumSum

    var start_index = 1
    for (j in 0 until stateSizeCount) {
        val size = statesz.get(j).toInt()
        val end_index = start_index + size
        val tmp = Matrix.extract(tranState, 0, tranState.numRows, start_index, end_index)
        tranSysState[j + 1] = tmp
        start_index = end_index
    }

    val arvRates = HashMap<Int, Matrix?>()
    val depRates = HashMap<Int, Matrix?>()
    repeat(R) { i ->
        val rate = Matrix(ui.numRows, sn.nstateful)
        arvRates[i] = rate
        depRates[i] = rate.copy()
    }

    var pi = Matrix(1, ui.numRows)
    for (s in 0 until ui.numRows) {
        val stateIndexes = uj[s]!!
        var sum = 0.0
        for (i in stateIndexes) {
            sum += tranState.get(i, 0)
        }
        pi.set(s, sum)
    }

    // Pre-allocate result matrix to avoid O(n²) concatenation overhead
    val numUniqueStates = ui.numElements
    val SSq_new = Matrix(SSq.numRows, numUniqueStates)
    for (col_ind in 0 until numUniqueStates) {
        val srcCol = ui.get(col_ind).toInt()
        for (row in 0 until SSq.numRows) {
            SSq_new.set(row, col_ind, SSq.get(row, srcCol))
        }
    }
    SSq = SSq_new.transpose()

    for (ind in 0 until sn.nnodes) {
        if (sn.isstateful.get(ind) == 1.0) {
            val isf = sn.nodeToStateful.get(ind).toInt()
            for (s in 0 until ui.numRows) {
                val uis = ui.get(s).toInt()
                for (r in 0 until R) {
                    val arvRate = arvRates[r]!!
                    val depRate = depRates[r]!!
                    arvRate.set(s, isf, arvRatesSamples[uis]!!.get(r, isf))
                    depRate.set(s, isf, depRatesSamples[uis]!!.get(r, isf))
                }
            }
        }
    }

    pi = Matrix.scaleMult(pi, 1.0 / pi.elementSum())

    //if (options.method.equals("para") || options.method.equals ("parallel")) {
    //    System.out.printf("SSA samples: %6d\n", samples_collected);
    //}
    return SSAValues(pi, SSq, arvRates, depRates, tranSysState, tranSync, sn)
}

fun save_log(dt: Double,
             cur_state: MutableMap<Int?, Matrix>,
             tranState: Matrix,
             samples_collected: Int,
             tranSync: Matrix,
             enabled_sync: MutableMap<Int?, Int?>,
             firing_ctr: Int,
             sn: NetworkStruct,
             nir: MutableMap<Int?, Matrix>,
             stateCell: MutableMap<Int?, Matrix>,
             SSq: Matrix,
             streamingCollector: Collector? = null,
             curTime: Double = 0.0) {
    // Build state vector directly without matrix concatenation
    // First compute total size
    var totalSize = 0
    for (ind in 0..<sn.nnodes) {
        if (cur_state.containsKey(ind)) {
            totalSize += cur_state.get(ind)!!.numElements
        }
    }

    // Set dt as first element, then copy state data directly to tranState
    val colIdx = samples_collected - 1
    tranState.set(0, colIdx, dt)
    var offset = 1
    for (ind in 0..<sn.nnodes) {
        if (cur_state.containsKey(ind)) {
            val row: Matrix = cur_state.get(ind)!!
            for (i in 0..<row.numElements) {
                tranState.set(offset + i, colIdx, row.get(i))
            }
            offset += row.numElements
        }
    }

    // add new row to tranSync if needed
//            if (samples_collected - 1 >= tranSync.getNumRows()) {
//                Matrix tranSync_new_row = new Matrix(1, tranSync.getNumCols());
//                tranSync_new_row.zero();
//                tranSync = Matrix.concatRows(tranSync, tranSync_new_row, null);
//            }
    // +1 to convert 0-based Java indexing to 1-based indexing for LINE compatibility
    tranSync.set(samples_collected - 1, 0, enabled_sync.get(firing_ctr)!! + 1)

    for (ind in 0..<sn.nnodes) {
        if (sn.isstation.get(ind) == 1.0) {
            val isf = sn.nodeToStateful.get(ind).toInt()
            val ist = sn.nodeToStation.get(ind).toInt()
            nir.put(ist, toMarginal(sn, ind, stateCell.get(isf), null, null, null, null, null).nir.transpose())
        }
    }

    // Copy nir values directly to SSq column without matrix concatenation
    // IMPORTANT: iterate by station index to ensure correct ordering in SSq
    val colIdx2 = samples_collected - 1
    var rowOffset = 0
    for (ist in 0..<sn.nstations) {
        if (nir.containsKey(ist)) {
            val col: Matrix = nir.get(ist)!!
            for (i in 0..<col.numRows) {
                SSq.set(rowOffset + i, colIdx2, col.get(i))
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

fun update_paddings(sn: NetworkStruct,
                    stateCell: MutableMap<Int?, Matrix>,
                    statelen: Matrix,
                    tranState: Matrix): Matrix {
    var tranState = tranState
    for (ind in 0..<sn.nnodes) {
        if (sn.isstation.get(ind) == 1.0) {
            val isf = sn.nodeToStateful.get(ind).toInt()
            val deltalen = (stateCell.get(isf)!!.numElements > statelen.get(isf))
            if (deltalen) {
                statelen.set(isf, stateCell.get(isf)!!.numElements.toDouble())
                // padding
                var shift = 0
                if (ind > 0) {
                    for (col in 0..<isf) {
                        shift += statelen.get(col).toInt()
                    }
                }
                val pad = Matrix(1, tranState.numCols)

                val top = Matrix.extractRows(tranState, 0, shift + 1, null)
                val bottom = Matrix.extractRows(tranState, shift + 1, tranState.numRows, null)
                val tmp = Matrix.concatRows(top, pad, null)
                tranState = Matrix.concatRows(tmp, bottom, null)
            }
        }
    }
    return tranState
}