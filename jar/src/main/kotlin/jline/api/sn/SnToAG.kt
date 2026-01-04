/**
 * @file Convert LINE network structure to Agent (RCAT) format for SolverAG
 *
 * Converts a LINE queueing network to the R, AP format required by
 * the INAP and AUTOCAT algorithms using the sn.sync data structure.
 *
 * @since LINE 3.0
 */
package jline.api.sn

import jline.io.Ret
import jline.lang.Event
import jline.lang.NetworkStruct
import jline.lang.constant.EventType
import jline.lang.constant.NodeType
import jline.lang.constant.SignalType
import jline.util.Utils
import jline.util.matrix.Matrix
import java.util.ArrayList

/**
 * Convert LINE network structure to Agent (RCAT) format for SolverAG
 *
 * @param sn NetworkStruct object for the queueing network model
 * @param maxStates Maximum number of states per process (truncation level), default 100
 * @return Ret.snToAG containing R, AP, processMap, actionMap, N
 */
@JvmOverloads
fun snToAG(sn: NetworkStruct, maxStates: Int = 100): Ret.snToAG {
    val M = sn.nstations
    val K = sn.nclasses
    val sync = sn.sync
    val local = sn.nnodes + 1  // local action indicator

    // Identify station types
    val sourceStations = ArrayList<Int>()
    val sinkStations = ArrayList<Int>()
    val queueStations = ArrayList<Int>()

    for (ist in 0 until M) {
        val nodeIdx = sn.stationToNode[ist, 0].toInt()
        when (sn.nodetype[nodeIdx]) {
            NodeType.Source -> sourceStations.add(ist)
            NodeType.Sink -> sinkStations.add(ist)
            else -> queueStations.add(ist)  // Queue, Delay, or other service stations
        }
    }

    // Create process mapping: each (station, class) pair at queue stations
    // In RCAT, each queue with a specific class is modeled as a separate process
    // Note: Signal classes do NOT form processes - they only trigger actions
    var processIdx = 0
    val processMap = Matrix(M, K)
    for (ist in queueStations) {
        for (r in 0 until K) {
            // Skip signal classes - they don't form processes
            if (sn.issignal != null && sn.issignal[r, 0] != 0.0) {
                continue
            }
            val rate = sn.rates[ist, r]
            if (!java.lang.Double.isNaN(rate) && rate > 0) {
                processIdx++
                processMap[ist, r] = processIdx.toDouble()
            }
        }
    }
    val numProcesses = processIdx

    if (numProcesses == 0) {
        return Ret.snToAG(
            arrayOf(),
            Matrix(1, 2),
            processMap,
            ArrayList(),
            intArrayOf()
        )
    }

    // Determine number of states for each process
    val N = IntArray(numProcesses)
    for (p in 0 until numProcesses) {
        // Find station and class for this process (1-indexed process)
        for (ist in 0 until M) {
            for (r in 0 until K) {
                if (processMap[ist, r].toInt() == p + 1) {
                    N[p] = if (Utils.isInf(sn.njobs[r])) {
                        maxStates  // Open class - truncate at maxStates
                    } else {
                        sn.njobs[r].toInt() + 1  // Closed class - states 0, 1, ..., njobs
                    }
                }
            }
        }
    }

    // Build actions from sn.sync
    // Each sync event with DEP (active) -> ARV (passive) between queue stations is an action
    val actionMapList = ArrayList<Ret.ActionMapEntry>()

    if (sync != null) {
        for (s in sync.keys) {
            val syncEntry = sync[s] ?: continue
            if (syncEntry.active.isEmpty() || syncEntry.passive.isEmpty()) {
                continue
            }

            val activeEvent = syncEntry.active[0] ?: continue
            val passiveEvent = syncEntry.passive[0] ?: continue

            // We only care about DEP -> ARV synchronizations between queue stations
            if (activeEvent.event != EventType.DEP || passiveEvent.event != EventType.ARV) {
                continue
            }

            // Get station indices
            val nodeA = activeEvent.node
            val nodeP = passiveEvent.node

            if (nodeP == local) {
                continue  // Local actions are handled separately
            }

            // Convert node to station
            if (sn.isstation[nodeA, 0] == 0.0 || sn.isstation[nodeP, 0] == 0.0) {
                continue
            }

            val ist = sn.nodeToStation[nodeA, 0].toInt()
            val jst = sn.nodeToStation[nodeP, 0].toInt()

            // Skip if either is source or sink
            if (sourceStations.contains(ist) || sinkStations.contains(ist)) {
                continue
            }
            if (sourceStations.contains(jst) || sinkStations.contains(jst)) {
                continue
            }

            val r = activeEvent.jobClass
            val sClass = passiveEvent.jobClass

            // Check if both processes exist
            if (processMap[ist, r] == 0.0 || processMap[jst, sClass] == 0.0) {
                continue
            }

            // Get routing probability
            val prob = if (passiveEvent.probFun != null) {
                1.0  // State-dependent routing - use 1 as placeholder
            } else {
                val rawProb = passiveEvent.prob
                if (java.lang.Double.isNaN(rawProb)) 1.0 else rawProb
            }

            if (prob <= 0) {
                continue
            }

            // Check if this is a negative signal class
            var isNegativeClass = false
            if (sn.issignal != null && sn.issignal[r, 0] != 0.0) {
                if (sn.signaltype != null && r < sn.signaltype.size) {
                    val signalType = sn.signaltype[r]
                    if (signalType == SignalType.NEGATIVE) {
                        isNegativeClass = true
                    }
                }
            }

            actionMapList.add(Ret.ActionMapEntry(ist, r, jst, sClass, prob, isNegativeClass))
        }
    }
    val numActions = actionMapList.size

    // Initialize R and AP
    val numCols = if (numProcesses < 2) 2 else numProcesses
    val R = Array(numActions + 1) { arrayOfNulls<Matrix>(numCols) }
    val AP = Matrix(if (numActions < 1) 1 else numActions, 2)

    // Build local/hidden rate matrices L for each process (R[numActions][k])
    for (p in 0 until numProcesses) {
        // Find station and class for this process
        var pStation = -1
        var pClass = -1
        outer@ for (ist in 0 until M) {
            for (r in 0 until K) {
                if (processMap[ist, r].toInt() == p + 1) {
                    pStation = ist
                    pClass = r
                    break@outer
                }
            }
        }
        if (pStation >= 0 && pClass >= 0) {
            R[numActions][p] = buildLocalRates(sn, pStation, pClass, N[p], sourceStations, K)
        }
    }

    // Build active and passive matrices for each action
    for (a in 0 until numActions) {
        val am = actionMapList[a]

        // Active process (departure)
        val ist = am.fromStation
        val r = am.fromClass
        val pActive = processMap[ist, r].toInt() - 1  // Convert to 0-indexed
        AP[a, 0] = (pActive + 1).toDouble()  // Store as 1-indexed

        // Get service rate
        val muIr = sn.rates[ist, r]
        val prob = am.prob

        // Active matrix Aa: rate of active transition
        // For a departure: state n -> n-1 with rate mu*prob
        val Aa = Matrix.zeros(N[pActive], N[pActive])
        for (n in 1 until N[pActive]) {
            Aa[n, n - 1] = muIr * prob
        }
        R[a][0] = Aa

        // Passive process (arrival or signal effect)
        val jst = am.toStation
        val sClass = am.toClass
        val pPassive = processMap[jst, sClass].toInt() - 1  // Convert to 0-indexed
        AP[a, 1] = (pPassive + 1).toDouble()  // Store as 1-indexed

        // Passive matrix Pb: probability (not rate!) of passive transition
        val Pb = Matrix.zeros(N[pPassive], N[pPassive])
        if (am.isNegative) {
            // NEGATIVE: Job removal at destination (G-network negative customer)
            // Empty queue: no effect (stay at state 0)
            Pb[0, 0] = 1.0
            // Non-empty queues: decrement (n -> n-1)
            for (n in 1 until N[pPassive]) {
                Pb[n, n - 1] = 1.0
            }
        } else {
            // POSITIVE: Normal job arrival at destination
            // Arrival: state n -> n+1
            for (n in 0 until N[pPassive] - 1) {
                Pb[n, n + 1] = 1.0
            }
            // Boundary: at max capacity, stay (job is lost or blocked)
            Pb[N[pPassive] - 1, N[pPassive] - 1] = 1.0
        }
        R[a][1] = Pb
    }

    // Convert R to Matrix[][] (need non-null matrices)
    @Suppress("UNCHECKED_CAST")
    val RResult = R.map { row -> row.map { it ?: Matrix.zeros(1, 1) }.toTypedArray() }.toTypedArray()

    return Ret.snToAG(RResult, AP, processMap, actionMapList, N)
}

/**
 * Build local/hidden transition matrix for process at station ist, class r
 *
 * Local transitions include:
 * - External arrivals (from source)
 * - Departures to sink (jobs leaving the system)
 *
 * Note: Departures to other queues are NOT local - they are handled as actions
 *
 * In LINE's open network model:
 * - Sink is not a stateful node
 * - Departures to sink are encoded as routing back to Source in sn.rt
 * - So we calculate prob_sink as routing from queue to source
 */
private fun buildLocalRates(
    sn: NetworkStruct,
    ist: Int,
    r: Int,
    Np: Int,
    sourceStations: List<Int>,
    K: Int
): Matrix {
    val L = Matrix.zeros(Np, Np)
    val rt = sn.rt

    // External arrivals from source
    var lambdaIr = 0.0
    for (isrc in sourceStations) {
        for (sSrc in 0 until K) {
            // Arrivals from source station isrc, class sSrc to station ist, class r
            val isfSrc = sn.nodeToStateful[sn.stationToNode[isrc, 0].toInt(), 0].toInt()
            val isfIst = sn.nodeToStateful[sn.stationToNode[ist, 0].toInt(), 0].toInt()

            val probSrc = rt[isfSrc * K + sSrc, isfIst * K + r]
            val srcRate = sn.rates[isrc, sSrc]
            if (probSrc > 0 && !java.lang.Double.isNaN(srcRate)) {
                // Check if source class sSrc is a negative signal - handle separately
                var isNegativeSignal = false
                if (sn.issignal != null && sn.issignal[sSrc, 0] != 0.0) {
                    if (sn.signaltype != null && sSrc < sn.signaltype.size) {
                        if (sn.signaltype[sSrc] == SignalType.NEGATIVE) {
                            isNegativeSignal = true
                        }
                    }
                }
                if (!isNegativeSignal) {
                    lambdaIr += srcRate * probSrc
                }
            }
        }
    }

    // External arrival transitions: n -> n+1 with rate lambda
    if (lambdaIr > 0) {
        for (n in 0 until Np - 1) {
            L[n, n + 1] = lambdaIr
        }
        // At capacity, arrivals are lost (or we could model blocking)
    }

    // Departures to sink (modeled as routing to source in LINE's open network representation)
    val muIr = sn.rates[ist, r]
    if (!java.lang.Double.isNaN(muIr) && muIr > 0) {
        // In LINE, for open networks, departures to sink are encoded as routing to source
        // So prob_sink = sum of routing probabilities from this queue to source stations
        var probSink = 0.0
        val isfIst = sn.nodeToStateful[sn.stationToNode[ist, 0].toInt(), 0].toInt()

        for (jsrc in sourceStations) {
            val jsf = sn.nodeToStateful[sn.stationToNode[jsrc, 0].toInt(), 0].toInt()
            for (s in 0 until K) {
                probSink += rt[isfIst * K + r, jsf * K + s]
            }
        }

        // Departure to sink transitions: n -> n-1 with rate mu*prob_sink
        if (probSink > 0) {
            for (n in 1 until Np) {
                L[n, n - 1] = L[n, n - 1] + muIr * probSink
            }
        }
    }

    return L
}

/**
 * Stochastic network ToAG algorithms
 */
@Suppress("unused")
class SntoagKt {
    companion object {
        // Class documentation marker for Dokka
    }
}
