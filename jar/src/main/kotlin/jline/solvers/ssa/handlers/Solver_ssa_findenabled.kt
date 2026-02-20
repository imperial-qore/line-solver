package jline.solvers.ssa.handlers

import jline.lang.NetworkStruct
import jline.lang.Sync
import jline.lang.constant.EventType
import jline.GlobalConstants
import jline.lang.constant.NodeType
import jline.lang.state.EventCache
import jline.lang.state.State
import jline.solvers.ssa.SolverSSA
import jline.util.matrix.Matrix

fun solver_ssa_findenabled(sn: NetworkStruct,
                           eventCache: EventCache,
                           A: Int,
                           node_a: MutableMap<Int?, Int?>,
                           next_state: MutableMap<Int?, MutableMap<Int?, Matrix>>,
                           stateCell: MutableMap<Int?, Matrix>,
                           event_a: MutableMap<Int?, EventType?>,
                           class_a: MutableMap<Int?, Int?>,
                           isSimulation: Boolean,
                           outprob_a: MutableMap<Int?, Double?>,
                           node_p: MutableMap<Int?, Int?>,
                           local: Int,
                           event_p: MutableMap<Int?, EventType?>,
                           class_p: MutableMap<Int?, Int?>,
                           outprob_p: MutableMap<Int?, Double?>,
                           prob_sync_p: MutableMap<Int?, Double?>,
                           sync: MutableMap<Int?, Sync?>,
                           node_a_sf: MutableMap<Int?, Int?>,
                           node_p_sf: MutableMap<Int?, Int?>,
                           depRatesSamples: MutableMap<Int?, Matrix?>,
                           samples_collected: Int,
                           arvRatesSamples: MutableMap<Int?, Matrix?>,
                           csmask: Matrix,
                           enabled_rates: MutableMap<Int?, Double?>,
                           enabled_sync: MutableMap<Int?, Int?>,
                           solverSSA: SolverSSA) {
    for (act in 0..<A) {
        val rate_a: MutableMap<Int?, Matrix?> = HashMap<Int?, Matrix?>()
        val isf_a = sn.nodeToStateful.get(node_a.get(act)!!).toInt()
        var isf_p: Int
        next_state[act] = stateCell.mapValues { it.value!!.copy() }.toMutableMap()

        solverSSA.run {
            val eventResult = State.afterEvent(sn,
                node_a.get(act)!!,
                stateCell.get(isf_a),
                event_a.get(act),
                class_a.get(act)!!,
                isSimulation,
                eventCache)
            if (!eventResult.outspace.isEmpty) {
                next_state.get(act)!!.put(sn.nodeToStateful.get(node_a.get(act)!!).toInt(), eventResult.outspace)
            } else {
                next_state.get(act)!!.remove(sn.nodeToStateful.get(node_a.get(act)!!).toInt())
            }
            if (!eventResult.outrate.isEmpty) {
                rate_a.put(act, eventResult.outrate)
            } else {
                rate_a.remove(act)
            }
            if (!eventResult.outprob.isEmpty) {
                outprob_a.put(act, eventResult.outprob.toDouble())
            } else {
                outprob_a.remove(act)
            }
        }

        if (!next_state.get(act)!!.containsKey(isf_a) || !rate_a.containsKey(act)) {
            continue
        }

        for (ia in 0..<next_state.get(act)!!.get(isf_a)!!.numRows) {
            if (java.lang.Double.isNaN(rate_a.get(act)!!.get(ia)) || rate_a.get(act)!!.get(ia) == 0.0) {
                // handling degenerate rate values
                rate_a.get(act)!!.set(ia, GlobalConstants.Zero)
            }

            val hash_check: Matrix = next_state.get(act)!!.get(isf_a)!!
            var hash_found = false
            for (col in 0..<hash_check.numCols) {
                if (hash_check.get(ia, col) != -1.0) {
                    hash_found = true
                    break
                }
            }

            if (!hash_found) {
                continue
            }

            //boolean update_cond = true;
            if (rate_a.get(act)!!.get(ia) > 0) {
                if (node_p.get(act) != local) {
                    isf_p = sn.nodeToStateful.get(node_p.get(act)!!).toInt()
                    if (node_p.get(act) == node_a.get(act)) {
                        // self-loop

                        val eventResult = State.afterEvent(sn,
                            node_p.get(act)!!,
                            next_state.get(act)!!.get(isf_a),
                            event_p.get(act),
                            class_p.get(act)!!,
                            isSimulation,
                            eventCache)
                        if (!eventResult.outspace.isEmpty) {
                            next_state.get(act)!!.put(isf_p, eventResult.outspace)
                        } else {
                            next_state.get(act)!!.remove(isf_p)
                        }
                        if (!eventResult.outprob.isEmpty) {
                            outprob_p.put(act, eventResult.outprob.toDouble())
                        }
                    } else {
                        // departure
                        val eventResult = State.afterEvent(sn,
                            node_p.get(act)!!,
                            next_state.get(act)!!.get(isf_p),
                            event_p.get(act),
                            class_p.get(act)!!,
                            isSimulation,
                            eventCache)

                        if (!eventResult.outspace.isEmpty) {
                            next_state.get(act)!!.put(isf_p, eventResult.outspace)
                        } else {
                            next_state.get(act)!!.remove(isf_p)
                        }
                        if (!eventResult.outprob.isEmpty) {
                            outprob_p.put(act, eventResult.outprob.toDouble())
                        }
                    }

                    if (next_state.get(act)!!.containsKey(isf_p)) {
                        // Check if the source node (node_a) has state-dependent routing
                        if (node_a.get(act)!! < sn.nnodes && sn.isstatedep.get(node_a.get(act)!!, 2) == 1.0) {
                            // State-dependent routing: probability depends on current and next state
                            // Convert stateCell and next_state maps from stateful index -> Matrix to Node -> Matrix
                            // IMPORTANT: Use sn.stateful (not model.getStatefulNodes()) to get the same Node
                            // objects that are captured by the routing lambda in Network.refreshRouting()
                            val stateCell_node: MutableMap<jline.lang.nodes.Node?, jline.util.matrix.Matrix?> = HashMap()
                            for (entry in stateCell.entries) {
                                val stateful_index = entry.key
                                val matrix = entry.value
                                if (stateful_index != null && stateful_index < sn.stateful.size) {
                                    val node = sn.stateful.get(stateful_index)
                                    if (node != null) {
                                        stateCell_node.put(node, matrix)
                                    }
                                }
                            }
                            val nextState_node: MutableMap<jline.lang.nodes.Node?, jline.util.matrix.Matrix?> = HashMap()
                            for (entry in next_state.get(act)!!.entries) {
                                val stateful_index = entry.key
                                val matrix = entry.value
                                if (stateful_index != null && stateful_index < sn.stateful.size) {
                                    val node = sn.stateful.get(stateful_index)
                                    if (node != null) {
                                        nextState_node.put(node, matrix)
                                    }
                                }
                            }
                            val nodePairs = jline.util.Pair(stateCell_node, nextState_node)
                            prob_sync_p.put(act, sync.get(act)!!.passive.get(0)!!.getProb(nodePairs))
                        } else {
                            prob_sync_p.put(act, sync.get(act)!!.passive.get(0)!!.getProb())
                        }
                    } else {
                        prob_sync_p.put(act, 0.0)
                    }
                }
                if (next_state.get(act)!!.containsKey(isf_a)) {
                    if (node_p.get(act) == local) {
                        prob_sync_p.put(act, 1.0)
                    }
                    if (!java.lang.Double.isNaN(rate_a.get(act)!!.toDouble())) {
                        if (next_state.get(act)!!.size == stateCell.size) {
                            if (event_a.get(act) == EventType.DEP) {
                                isf_p = sn.nodeToStateful.get(node_p.get(act)!!).toInt()
                                node_a_sf.put(act, isf_a)
                                node_p_sf.put(act, isf_p)

                                //                                        Matrix original_departure = depRatesSamples.get(class_a.get(act));
//                                        Matrix original_arrival = arvRatesSamples.get(class_p.get(act));
                                val added_value = (outprob_a.get(act)!! * outprob_p.get(act)!! * rate_a.get(act)!!
                                    .get(ia) * prob_sync_p.get(act)!!)

                                val a_sf_act: Int = node_a_sf.get(act)!!
                                val p_sf_act: Int = node_p_sf.get(act)!!

                                val dep_value = (depRatesSamples.get(samples_collected - 1)!!
                                    .get(class_a.get(act)!!, a_sf_act) + added_value)
                                val arv_val = (arvRatesSamples.get(samples_collected - 1)!!
                                    .get(class_p.get(act)!!, p_sf_act) + added_value)

                                depRatesSamples.get(samples_collected - 1)!!
                                    .set(class_a.get(act)!!, a_sf_act, dep_value)
                                arvRatesSamples.get(samples_collected - 1)!!.set(class_p.get(act)!!, p_sf_act, arv_val)
                            }
                            if (node_p.get(act)!! < local && csmask.get(class_a.get(act)!!,
                                    class_p.get(act)!!) != 1.0 && sn.nodetype.get(node_p.get(act)!!) != NodeType.Source && (rate_a.get(
                                    act)!!.get(ia) * prob_sync_p.get(act)!! > 0)) {
                                // Error: state-dependent routing violates the class switching mask
                                throw RuntimeException("Error: state-dependent routing at node ${node_a.get(act)} violates the class switching mask (node ${node_a.get(act)} -> node ${node_p.get(act)}, class ${class_a.get(act)} -> class ${class_p.get(act)}).")
                            }

                            val ctr = enabled_rates.size
                            enabled_rates.put(ctr, rate_a.get(act)!!.get(ia) * prob_sync_p.get(act)!!)
                            enabled_sync.put(ctr, act)
                        }
                    }
                }
            }
        }
    }
}
