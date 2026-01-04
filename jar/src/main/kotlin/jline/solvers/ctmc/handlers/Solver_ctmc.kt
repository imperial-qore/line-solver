package jline.solvers.ctmc.handlers

import jline.api.mc.ctmc_makeinfgen
import jline.api.mc.ctmc_ssg
import jline.api.mc.ctmc_stochcomp
import jline.lang.NetworkStruct
import jline.lang.constant.EventType
import jline.GlobalConstants
import jline.lang.constant.NodeType
import jline.VerboseLevel
import jline.api.mc.ctmc_ssg_reachability
import jline.lang.nodes.Node
import jline.lang.state.State
import jline.solvers.SolverOptions
import jline.solvers.ctmc.ResultCTMC
import jline.solvers.ctmc.SolverCTMC
import jline.util.Pair
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import java.lang.Double
import kotlin.Array
import kotlin.DoubleArray
import kotlin.IllegalArgumentException
import kotlin.Int
import kotlin.checkNotNull

class Solver_ctmc(private val solverCTMC: SolverCTMC?) {
    companion object {
        @JvmStatic
        fun solver_ctmc(sn: NetworkStruct, options: SolverOptions): ResultCTMC {
            var sn = sn
            val nstateful = sn.nstateful
            val nclasses = sn.nclasses
            val sync = sn.sync
            val A = sync.size
            val csmask = sn.csmask
            if (options.config.state_space_gen == null) {
                options.config.state_space_gen = "default"
            }

            var stateSpace: Matrix
            var stateSpaceAggr: Matrix
            var stateSpaceHashed: Matrix
            when (options.config.state_space_gen) {
                "reachable" -> {
                    // Note: reachable state space generation does not handle open models yet (no cutoff)
                    val ssgResult = ctmc_ssg_reachability(sn, options)
                    stateSpace = ssgResult.stateSpace
                    stateSpaceAggr = ssgResult.stateSpaceAggr
                    stateSpaceHashed = ssgResult.stateSpaceHashed
                    sn = ssgResult.sn
                }
                "default", "full" -> {
                    val ssgResult = ctmc_ssg(sn, options)
                    stateSpace = ssgResult.getStateSpace()
                    stateSpaceAggr = ssgResult.getStateSpaceAggr()
                    stateSpaceHashed = ssgResult.getStateSpaceHashed()
                    sn = ssgResult.getSn()
                }

                else -> throw IllegalArgumentException("Unknown state space generation method: " + options.config.state_space_gen)
            }

            //    Matrix Q = sparse()
            val size = stateSpaceHashed.getNumRows()
            var Q = Matrix.eye(size)
            val Dfilt = MatrixCell()
            for (a in 0..<A) {
                Dfilt.set(a, Matrix(size, size))
            }
            val local = sn.nnodes + 1
            for (a in 0..<A) {
                val stateCell = MatrixCell()
                for (s in 0..<stateSpaceHashed.getNumRows()) {
                    val state = stateSpaceHashed.getRow(s)

                    for (ind in 0..<sn.nnodes) {
                        if (sn.isstateful.get(ind, 0) == 1.0) {
                            val isf = sn.nodeToStateful.get(ind).toInt()
                            val stateIndex = state.get(isf).toInt()

                            val spaceMatrix: Matrix = sn.space.get(sn.stateful.get(isf))!!
                            val stateRow = spaceMatrix.getRow(stateIndex)
                            stateCell.set(isf, stateRow)
                        }
                    }

                    val node_a = sync.get(a)!!.active.get(0)!!.getNode()
                    val state_a = state.get(sn.nodeToStateful.get(node_a).toInt())
                    val class_a = sync.get(a)!!.active.get(0)!!.getJobClass()
                    val event_a = sync.get(a)!!.active.get(0)!!.getEvent()

                    val eventResult = State.afterEventHashed(sn, node_a, state_a, event_a, class_a)
                    val new_state_a = eventResult.outspace
                    val rate_a = eventResult.outrate

                    if (new_state_a.get(0) == -1.0) {
                        continue
                    }

                    for (ia in 0..<new_state_a.length()) {
                        if (rate_a.get(ia) > 0) {
                            val node_p = sync.get(a)!!.passive.get(0)!!.getNode()
                            if (node_p + 1 != local) {
                                val state_p = state.get(sn.nodeToStateful.get(node_p).toInt()).toInt()
                                val class_p = sync.get(a)!!.passive.get(0)!!.getJobClass()
                                val event_p = sync.get(a)!!.passive.get(0)!!.getEvent()


                                var afterEventResult: jline.io.Ret.EventResult? = null
                                if (node_p == node_a) {
                                    if (new_state_a.get(ia) != -1.0) {
                                        afterEventResult =
                                            State.afterEventHashed(sn, node_p, new_state_a.get(ia), event_p, class_p)
                                    }
                                } else {
                                    if (new_state_a.get(ia) != -1.0) {
                                        afterEventResult =
                                            State.afterEventHashed(sn, node_p, state_p.toDouble(), event_p, class_p)
                                    }
                                }

                                if (afterEventResult == null) {
                                    continue;
                                }

                                val new_state_p = afterEventResult.outspace
                                val outprob_p = afterEventResult.outprob

                                for (ip in 0..<new_state_p.getNumRows()) {
                                    var prob_sync_p = 0.0

                                    // Skip invalid states (hash == -1 means state not found in state space)
                                    if (new_state_p.get(ip) == -1.0) {
                                        continue
                                    }

                                    if (node_p + 1 != local) {
                                        if (new_state_p.get(ip) != -1.0) {
                                            if (sn.isstatedep.get(node_a, 2) != 0.0) {
                                                val newStateCell = MatrixCell(stateCell)
                                                val statefulNodeA = sn.nodeToStateful.get(node_a).toInt()
                                                val statefulNodeP = sn.nodeToStateful.get(node_p).toInt()

                                                val spaceMatrixA = sn.space.get(sn.stateful.get(statefulNodeA))!!
                                                    .getRow(new_state_a.get(ia).toInt())
                                                val spaceMatrixP = sn.space.get(sn.stateful.get(statefulNodeP))!!
                                                    .getRow(new_state_p.get(ip).toInt())

                                                newStateCell.set(statefulNodeA, spaceMatrixA)

                                                newStateCell.set(statefulNodeP, spaceMatrixP)

                                                val stateCell_node: MutableMap<Node?, Matrix?> =
                                                    HashMap<Node?, Matrix?>()
                                                for (entry in stateCell.toMap().entries) {
                                                    val node_index = entry.key
                                                    val matrix = entry.value
                                                    val node = sn.nodes.get(node_index)
                                                    if (node != null) {
                                                        stateCell_node.putIfAbsent(node, matrix)
                                                    }
                                                }

                                                val newStateCell_node: MutableMap<Node?, Matrix?> =
                                                    HashMap<Node?, Matrix?>()
                                                for (entry in newStateCell.toMap().entries) {
                                                    val node_index = entry.key
                                                    val matrix = entry.value
                                                    val node = sn.nodes.get(node_index)
                                                    if (node != null) {
                                                        newStateCell_node.putIfAbsent(node, matrix)
                                                    }
                                                }

                                                val nodePairs =
                                                    Pair<MutableMap<Node?, Matrix?>?, MutableMap<Node?, Matrix?>?>(
                                                        stateCell_node,
                                                        newStateCell_node
                                                    )

                                                prob_sync_p = sync.get(a)!!.passive.get(0)!!
                                                    .getProb(nodePairs) * outprob_p.get(ip)
                                            } else {
                                                prob_sync_p =
                                                    sync.get(a)!!.passive.get(0)!!.getProb() * outprob_p.get(ip)
                                            }
                                        } else {
                                            prob_sync_p = 0.0
                                        }
                                    }
                                    //               Matlab line 149
                                    var new_state: Matrix? = null
                                    if (!Double.isNaN(new_state_a.get(ia))) {
                                        if (node_p + 1 == local) {
                                            new_state = state.copy()
                                            new_state.set(sn.nodeToStateful.get(node_a).toInt(), new_state_a.get(ia))
                                            prob_sync_p = outprob_p.get(ip)
                                        } else if (!new_state_p.isEmpty()) {
                                            new_state = state.copy()
                                            new_state.set(sn.nodeToStateful.get(node_a).toInt(), new_state_a.get(ia))
                                            new_state.set(sn.nodeToStateful.get(node_p).toInt(), new_state_p.get(ip))
                                        }
                                        //                  MATLAB line 192
                                        checkNotNull(new_state)
                                        val ns = Matrix.matchrow(stateSpaceHashed, new_state)


                                        if (ns >= 0) {
                                            if (!rate_a.isEmpty()) {
                                                if (node_p + 1 < local && csmask.get(
                                                        class_a,
                                                        class_p
                                                    ) == 0.0 && rate_a.get(ia) * prob_sync_p > 0 && sn.nodetype.get(
                                                        node_p
                                                    ) != NodeType.Source
                                                ) {
                                                    System.err.printf(
                                                        "Error: state-dependent routing at node %d (%s) violates the class switching mask (node %s -> node %s, class %s -> class %s).",
                                                        node_a,
                                                        sn.nodenames.get(node_a),
                                                        sn.nodenames.get(node_a),
                                                        sn.nodenames.get(node_p),
                                                        sn.classnames.get(class_a),
                                                        sn.classnames.get(class_p)
                                                    )
                                                }

                                                val finalRate = rate_a.get(ia) * prob_sync_p

                                                if (Dfilt.get(a).getNumRows() >= s + 1 && Dfilt.get(a)
                                                        .getNumCols() >= ns + 1
                                                ) {
                                                    Dfilt.get(a).set(
                                                        s,
                                                        ns,
                                                        Dfilt.get(a).get(s, ns) + finalRate
                                                    )
                                                } else {
                                                    Dfilt.get(a).set(s, ns, finalRate)
                                                }
                                            }
                                        }
                                    }
                                }
                            } else {
                                if (!Double.isNaN(new_state_a.get(ia))) {
                                    val new_state = state.copy()
                                    new_state.set(sn.nodeToStateful.get(node_a).toInt(), new_state_a.get(ia))
                                    val prob_sync_p = 1.0
                                    val ns = Matrix.matchrow(stateSpaceHashed, new_state)
                                    if (ns > -1) {
                                        if (!rate_a.hasNaN()) {
                                            if (Dfilt.get(a).getNumRows() >= s && Dfilt.get(a).getNumCols() >= ns) {
                                                Dfilt.get(a)
                                                    .set(s, ns, Dfilt.get(a).get(s, ns) + rate_a.get(ia) * prob_sync_p)
                                            } else {
                                                Dfilt.get(a).set(s, ns, rate_a.get(ia) * prob_sync_p)
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            for (a in 0..<A) {
                Q = Q.add(1.0, Dfilt.get(a))
            }

            var diag_Q: Matrix? = Matrix(Q)
            Matrix.extractDiag(Q, diag_Q)
            val colMatrix = diag_Q!!.getColumn(0)
            diag_Q = Matrix.diag(*colMatrix.toArray1D())
            Q = Q.sub(diag_Q)

            var arvRates: Array<Array<DoubleArray?>?>? = Array<Array<DoubleArray?>?>(stateSpaceHashed.getNumRows()) {
                Array<DoubleArray?>(nstateful) {
                    DoubleArray(nclasses)
                }
            }
            var depRates: Array<Array<DoubleArray?>?>? = Array<Array<DoubleArray?>?>(stateSpaceHashed.getNumRows()) {
                Array<DoubleArray?>(nstateful) {
                    DoubleArray(nclasses)
                }
            }
            for (a in 0..<A) {
                val node_a = sync.get(a)!!.active.get(0)!!.getNode()
                val class_a = sync.get(a)!!.active.get(0)!!.getJobClass()
                val event_a = sync.get(a)!!.active.get(0)!!.getEvent()

                val node_p = sync.get(a)!!.passive.get(0)!!.getNode()
                val class_p = sync.get(a)!!.passive.get(0)!!.getJobClass()
                if (event_a == EventType.DEP) {
                    val node_a_sf = sn.nodeToStateful.get(node_a).toInt()
                    val node_p_sf = sn.nodeToStateful.get(node_p).toInt()
                    for (s in 0..<stateSpaceHashed.getNumRows()) {
                        val value = depRates!![s]!![node_a_sf]!![class_a] + Dfilt.get(a).sumRows(s)
                        depRates[s]!![node_a_sf]!![class_a] =
                            depRates[s]!![node_a_sf]!![class_a] + Dfilt.get(a).sumRows(s)
                        arvRates!![s]!![node_p_sf]!![class_p] =
                            arvRates[s]!![node_p_sf]!![class_p] + Dfilt.get(a).sumRows(s)
                    }
                }
            }

            val zero_row = Matrix.findIndexWithZeroSum(Q, true)
            val zero_col = Matrix.findIndexWithZeroSum(Q, false)

            Q.expandMatrixToSquare()
            val negativeIdentity_row = Matrix.eye(zero_row.size)
            negativeIdentity_row.mulByMinusOne()

            for ((rowIdx, row) in zero_row.withIndex()) {
                for ((colIdx, col) in zero_row.withIndex()) {
                    Q.set(row!!, col!!, negativeIdentity_row.get(rowIdx, colIdx))
                }
            }

            val negativeIdentity_col = Matrix.eye(zero_col.size)
            negativeIdentity_col.mulByMinusOne()

            for ((rowIdx, row) in zero_col.withIndex()) {
                for ((colIdx, col) in zero_col.withIndex()) {
                    Q.set(row!!, col!!, negativeIdentity_col.get(rowIdx, colIdx))
                }
            }

            //    MATLAB line 267
            for (a in 0..<A) {
                val colBound = Dfilt.get(a).getNumRows() - 1
                for (col in Dfilt.get(a).getNumCols()..<colBound) {
                    for (row in 0..<Dfilt.get(a).getNumRows()) {
                        Dfilt.get(a).set(row, col, 0)
                    }
                }
            }

            Q = ctmc_makeinfgen(Q)

            if (options.config.hide_immediate) {
                // Find non-station stateful nodes, excluding Cache and Router nodes
                // Cache nodes need their immediate transitions to compute hit/miss rates
                // Router nodes must preserve immediate transitions for correct throughput computation
                val imm_list: MutableList<kotlin.Double?> = ArrayList<kotlin.Double?>()

                for (ind in 0..<sn.nnodes) {
                    // Only consider non-station stateful nodes that are not Caches or Routers
                    if (sn.isstateful.get(ind) != 0.0 && sn.isstation.get(ind) == 0.0 && sn.nodetype.get(ind) != NodeType.Cache && sn.nodetype.get(ind) != NodeType.Router) {
                        val isf = sn.nodeToStateful.get(ind).toInt()

                        val space_slice = Matrix.extract(
                            sn.space.get(sn.stateful.get(isf)),
                            0,
                            sn.space.get(sn.stateful.get(isf))!!.getNumRows(),
                            0,
                            nclasses
                        )
                        val colSum = space_slice.sumCols()
                        val imm_st: MutableList<Int?> = ArrayList<Int?>()
                        for (row in 0..<colSum.getNumRows()) {
                            if (colSum.get(row, 0) > 0) {
                                imm_st.add(row)
                            }
                        }

                        for (s in 0..<stateSpaceHashed.getNumRows()) {
                            val finalIsf = isf
                            val finalS = s
                            val anyMatch = imm_st.stream().anyMatch { immStIndex: Int? ->
                                stateSpaceHashed.get(finalS, finalIsf) == immStIndex!!.toDouble()
                            }
                            if (anyMatch) {
                                imm_list.add(s.toDouble())
                            }
                        }
                    }
                }

                // Convert imm_list to unique sorted values
                val imm_unique = imm_list.filterNotNull().distinct().sorted()
                val imm = Matrix(if (imm_unique.isEmpty()) 0 else imm_unique.size, 1)
                for (i in imm_unique.indices) {
                    imm.set(i, 0, imm_unique[i])
                }

                val nonimm: MutableList<kotlin.Double?> = ArrayList<kotlin.Double?>()
                val allStates = (0..<Q.getNumRows()).map { it.toDouble() }

                if (imm_unique.isNotEmpty()) {
                    for (state in allStates) {
                        if (!imm_unique.contains(state)) {
                            nonimm.add(state)
                        }
                    }
                } else {
                    for (state in allStates) {
                        nonimm.add(state)
                    }
                }

                // Remove immediate state rows from stateSpace, stateSpaceAggr and stateSpaceHashed
                if (imm_unique.isNotEmpty()) {
                    // Create new matrices without the immediate state rows
                    val newStateSpace = Matrix(nonimm.size, stateSpace.getNumCols())
                    val newStateSpaceAggr = Matrix(nonimm.size, stateSpaceAggr.getNumCols())
                    val newStateSpaceHashed = Matrix(nonimm.size, stateSpaceHashed.getNumCols())
                    
                    for (i in nonimm.indices) {
                        val origRow = nonimm[i]!!.toInt()
                        for (col in 0..<stateSpace.getNumCols()) {
                            newStateSpace.set(i, col, stateSpace.get(origRow, col))
                        }
                        for (col in 0..<stateSpaceAggr.getNumCols()) {
                            newStateSpaceAggr.set(i, col, stateSpaceAggr.get(origRow, col))
                        }
                        for (col in 0..<stateSpaceHashed.getNumCols()) {
                            newStateSpaceHashed.set(i, col, stateSpaceHashed.get(origRow, col))
                        }
                    }
                    
                    stateSpace = newStateSpace
                    stateSpaceAggr = newStateSpaceAggr
                    stateSpaceHashed = newStateSpaceHashed

                    val stochcompResult = ctmc_stochcomp(Q, nonimm)
                    Q = stochcompResult.S
                    val Q12 = stochcompResult.Q12
                    val Q22 = stochcompResult.Q22

                    // Apply stochastic complement to event filters
                    for (a in 0..<A) {
                        val Q21a = Matrix(imm_unique.size, nonimm.size)
                        for (row in imm_unique.indices) {
                            val rowIndex = imm_unique[row].toInt()
                            for (col in nonimm.indices) {
                                val colIndex = nonimm[col]!!.toInt()
                                Q21a.set(row, col, Dfilt.get(a).get(rowIndex, colIndex))
                            }
                        }

                        Q22.changeSign()
                        val Q22_inv = Q22.inv()
                        var Ta = Q22_inv.mult(Q21a)
                        Ta = Q12.mult(Ta)
                        val dfilt_value = Matrix(nonimm.size, nonimm.size)

                        for (row in nonimm.indices) {
                            val rowIndex = nonimm[row]!!.toInt()
                            for (col in nonimm.indices) {
                                val colIndex = nonimm[col]!!.toInt()
                                dfilt_value.set(row, col, Dfilt.get(a).get(rowIndex, colIndex))
                            }
                        }
                        
                        if (Ta.getNumNonZeros() > 0) {
                            dfilt_value.add(Ta)
                        }
                        
                        Dfilt.set(a, dfilt_value)
                    }
                    
                    // Remove rows for immediate states from depRates and arvRates
                    val newDepRates: Array<Array<DoubleArray?>?> = Array<Array<DoubleArray?>?>(nonimm.size) {
                        Array<DoubleArray?>(nstateful) {
                            DoubleArray(nclasses)
                        }
                    }
                    val newArvRates: Array<Array<DoubleArray?>?> = Array<Array<DoubleArray?>?>(nonimm.size) {
                        Array<DoubleArray?>(nstateful) {
                            DoubleArray(nclasses)
                        }
                    }
                    
                    for (i in nonimm.indices) {
                        val origRow = nonimm[i]!!.toInt()
                        for (j in 0..<nstateful) {
                            for (k in 0..<nclasses) {
                                newDepRates[i]!![j]!![k] = depRates!![origRow]!![j]!![k]
                                newArvRates[i]!![j]!![k] = arvRates!![origRow]!![j]!![k]
                            }
                        }
                    }
                    
                    depRates = newDepRates
                    arvRates = newArvRates
                } else {
                    // No immediate states, keep the original matrices
                }
            }

            return ResultCTMC(Q, stateSpace, stateSpaceAggr, Dfilt, arvRates, depRates, sn)
        }
    }
}
