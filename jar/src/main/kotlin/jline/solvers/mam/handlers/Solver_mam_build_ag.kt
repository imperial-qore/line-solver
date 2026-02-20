package jline.solvers.mam.handlers

import jline.lang.NetworkStruct
import jline.lang.constant.NodeType
import jline.lang.constant.SignalType
import jline.lang.processes.DiscreteDistribution
import jline.util.matrix.Matrix
import kotlin.math.max
import kotlin.math.min

data class ActionInfo(
    val fromStation: Int,
    val fromClass: Int,
    val toStation: Int,
    val toClass: Int,
    val prob: Double,
    val isNegative: Boolean = false,
    val isCatastrophe: Boolean = false,
    val removalDistribution: DiscreteDistribution? = null
)

data class RCATModel(
    val R: Array<Array<Matrix?>>,
    val AP: Matrix,
    val processMap: Matrix,
    val actionMap: List<ActionInfo>,
    val N: IntArray
)

fun solver_mam_build_ag(sn: NetworkStruct, maxStates: Int): RCATModel {
    val M = sn.nstations
    val K = sn.nclasses
    val rt = sn.rt

    val sourceStations = mutableListOf<Int>()
    val queueStations = mutableListOf<Int>()

    // Find sink nodes (using node indices, not station indices)
    // Sinks are nodes, not stations, so we search in nodetype
    val sinkNodes = mutableListOf<Int>()
    for (nodeIdx in 0 until sn.nodetype.size) {
        if (sn.nodetype[nodeIdx] == NodeType.Sink) {
            sinkNodes.add(nodeIdx)
        }
    }

    for (ist in 0 until M) {
        val nodeIdx = sn.stationToNode.get(ist).toInt()
        when (sn.nodetype[nodeIdx]) {
            NodeType.Source -> sourceStations.add(ist)
            else -> queueStations.add(ist)
        }
    }

    var processIdx = 0
    val processMap = Matrix(M, K)
    processMap.fill(-1.0)

    for (ist in queueStations) {
        for (r in 0 until K) {
            // Skip Signal classes - they don't have their own queue state
            if (sn.issignal != null && sn.issignal.get(r, 0) > 0) {
                continue
            }
            val rate = sn.rates.get(ist, r)
            if (!rate.isNaN() && rate > 0) {
                processMap.set(ist, r, processIdx.toDouble())
                processIdx++
            }
        }
    }
    val numProcesses = processIdx

    if (numProcesses == 0) {
        return RCATModel(
            R = Array(1) { arrayOfNulls<Matrix>(1) },
            AP = Matrix(0, 2),
            processMap = processMap,
            actionMap = emptyList(),
            N = IntArray(0)
        )
    }

    val N = IntArray(numProcesses)
    for (p in 0 until numProcesses) {
        outer@ for (ist in 0 until M) {
            for (r in 0 until K) {
                if (processMap.get(ist, r).toInt() == p) {
                    val njobs = sn.njobs.get(r)
                    if (njobs.isInfinite()) {
                        N[p] = maxStates
                    } else {
                        N[p] = njobs.toInt() + 1
                    }
                    break@outer
                }
            }
        }
    }

    val actionMap = mutableListOf<ActionInfo>()

    for (ist in queueStations) {
        for (r in 0 until K) {
            if (processMap.get(ist, r).toInt() >= 0) {
                // Check if class r is a negative signal class
                val isNegativeClass = sn.issignal != null &&
                                      sn.issignal.get(r, 0) > 0 &&
                                      sn.signaltype != null &&
                                      sn.signaltype[r] == SignalType.NEGATIVE

                // Check if class r is a catastrophe signal
                val isCatastropheClass = sn.isCatastrophe != null &&
                                          sn.isCatastrophe.get(r, 0) > 0

                // Get removal distribution for this class (null for single removal or catastrophe)
                val removalDist: DiscreteDistribution? = if (sn.signalRemovalDist != null && r < sn.signalRemovalDist.size) {
                    sn.signalRemovalDist[r]
                } else {
                    null
                }

                for (jst in queueStations) {
                    for (s in 0 until K) {
                        if (processMap.get(jst, s).toInt() >= 0) {
                            val probIJRS = rt.get(ist * K + r, jst * K + s)
                            if (probIJRS > 0 && (ist != jst || r != s)) {
                                actionMap.add(ActionInfo(ist, r, jst, s, probIJRS, isNegativeClass, isCatastropheClass, removalDist))
                            }
                        }
                    }
                }
            }
        }
    }
    val numActions = actionMap.size

    val R = Array(numActions + 1) { arrayOfNulls<Matrix>(max(numProcesses, 2)) }
    val AP = Matrix(max(numActions, 1), 2)

    for (p in 0 until numProcesses) {
        for (ist in 0 until M) {
            for (r in 0 until K) {
                if (processMap.get(ist, r).toInt() == p) {
                    R[numActions][p] = buildLocalRates(sn, ist, r, N[p], rt, sourceStations, sinkNodes, K)
                    break
                }
            }
            if (R[numActions][p] != null) break
        }
    }

    for ((a, am) in actionMap.withIndex()) {
        val ist = am.fromStation
        val r = am.fromClass
        val pActive = processMap.get(ist, r).toInt()
        AP.set(a, 0, pActive.toDouble())

        val muIr = sn.rates.get(ist, r)
        val prob = am.prob

        val Aa = Matrix(N[pActive], N[pActive])
        for (n in 1 until N[pActive]) {
            Aa.set(n, n - 1, muIr * prob)
        }
        if (N[pActive] > 0) {
            Aa.set(N[pActive] - 1, N[pActive] - 1, muIr * prob)
        }
        R[a][0] = Aa

        val jst = am.toStation
        val s = am.toClass
        val pPassive = processMap.get(jst, s).toInt()
        AP.set(a, 1, pPassive.toDouble())

        val Pb = Matrix(N[pPassive], N[pPassive])
        if (am.isNegative) {
            // NEGATIVE: Job removal at destination (G-network negative customer)
            if (am.isCatastrophe) {
                // CATASTROPHE: All jobs are removed - all states transition to 0
                for (n in 0 until N[pPassive]) {
                    Pb.set(n, 0, 1.0)
                }
            } else if (am.removalDistribution != null) {
                // BATCH REMOVAL: Remove a random number of jobs based on distribution
                // P[n, m] = probability of transition from n to m jobs
                // = sum over k of P(remove k) where m = max(0, n - k)
                val dist = am.removalDistribution
                for (n in 0 until N[pPassive]) {
                    if (n == 0) {
                        // Empty queue: no effect
                        Pb.set(0, 0, 1.0)
                    } else {
                        // For each possible resulting state m (from 0 to n)
                        for (m in 0..n) {
                            val k = n - m  // Number of jobs to remove to go from n to m
                            // Probability of removing exactly k jobs when queue has n jobs
                            // = P(removal = k) for k < n, P(removal >= n) for k = n (i.e., m = 0)
                            if (m > 0) {
                                // Remove exactly k jobs: P(removal = k)
                                val prob = dist.evalPMF(k.toDouble())
                                if (prob > 0) {
                                    Pb.set(n, m, Pb.get(n, m) + prob)
                                }
                            } else {
                                // Remove all jobs (m = 0): P(removal >= n)
                                // = sum_{j >= n} P(removal = j)
                                // = 1 - CDF(n-1) = 1 - sum_{j=0}^{n-1} P(removal = j)
                                var cdfNMinus1 = 0.0
                                for (j in 0 until n) {
                                    cdfNMinus1 += dist.evalPMF(j.toDouble())
                                }
                                val probAtLeastN = 1.0 - cdfNMinus1
                                if (probAtLeastN > 0) {
                                    Pb.set(n, 0, Pb.get(n, 0) + probAtLeastN)
                                }
                            }
                        }
                    }
                }
            } else {
                // DEFAULT: Remove exactly 1 job (original behavior)
                // Empty queue: no effect (state 0 stays at state 0)
                Pb.set(0, 0, 1.0)
                // Non-empty queues: decrement (n -> n-1)
                for (n in 1 until N[pPassive] - 1) {
                    Pb.set(n, n - 1, 1.0)
                }
                // Boundary at max capacity: decrement
                if (N[pPassive] > 1) {
                    Pb.set(N[pPassive] - 1, N[pPassive] - 2, 1.0)
                }
            }
        } else {
            // POSITIVE: Normal job arrival at destination
            for (n in 0 until N[pPassive] - 1) {
                Pb.set(n, n + 1, 1.0)
            }
            // Boundary: at max capacity
            if (N[pPassive] > 0) {
                Pb.set(N[pPassive] - 1, N[pPassive] - 1, 1.0)
            }
        }
        R[a][1] = Pb
    }

    return RCATModel(R, AP, processMap, actionMap, N)
}

private fun buildLocalRates(
    sn: NetworkStruct,
    ist: Int,
    r: Int,
    Np: Int,
    rt: Matrix,
    sourceStations: List<Int>,
    sinkNodes: List<Int>,  // Node indices for sink nodes (not station indices)
    K: Int
): Matrix {
    val L = Matrix(Np, Np)

    // External arrivals from source - separate positive, negative (single), batch, and catastrophe
    var lambdaIrPos = 0.0  // Positive arrivals
    var lambdaIrNegSingle = 0.0  // Negative arrivals with single removal (default)
    var lambdaIrCatastrophe = 0.0  // Catastrophe arrivals (remove all)
    // Batch removal arrivals: list of (rate, distribution) pairs
    val batchArrivals = mutableListOf<Pair<Double, DiscreteDistribution>>()

    for (isrc in sourceStations) {
        for (sSrc in 0 until K) {
            // For signals: check routing to ANY class at this station (signals route as Signal->Signal
            // but their effect is on positive customer processes)
            val isSignalSrc = sn.issignal != null && sn.issignal.get(sSrc, 0) > 0
            val probSrc: Double
            if (isSignalSrc) {
                var pSum = 0.0
                for (sDst in 0 until K) {
                    pSum += rt.get(isrc * K + sSrc, ist * K + sDst)
                }
                probSrc = pSum
            } else {
                probSrc = rt.get(isrc * K + sSrc, ist * K + r)
            }
            val srcRate = sn.rates.get(isrc, sSrc)
            if (probSrc > 0 && !srcRate.isNaN()) {
                // Check if source class sSrc is a negative signal
                val isNegSrc = sn.issignal != null &&
                               sn.issignal.get(sSrc, 0) > 0 &&
                               sn.signaltype != null &&
                               sn.signaltype[sSrc] == SignalType.NEGATIVE
                if (isNegSrc) {
                    // Check if it's a catastrophe
                    val isCat = sn.isCatastrophe != null &&
                                sn.isCatastrophe.get(sSrc, 0) > 0
                    if (isCat) {
                        lambdaIrCatastrophe += srcRate * probSrc
                    } else {
                        // Check if it has a removal distribution
                        val removalDist: DiscreteDistribution? = if (sn.signalRemovalDist != null && sSrc < sn.signalRemovalDist.size) {
                            sn.signalRemovalDist[sSrc]
                        } else {
                            null
                        }
                        if (removalDist != null) {
                            batchArrivals.add(Pair(srcRate * probSrc, removalDist))
                        } else {
                            lambdaIrNegSingle += srcRate * probSrc
                        }
                    }
                } else {
                    lambdaIrPos += srcRate * probSrc
                }
            }
        }
    }

    // Positive arrival transitions: n -> n+1
    if (lambdaIrPos > 0) {
        for (n in 0 until Np - 1) {
            L.set(n, n + 1, lambdaIrPos)
        }
    }

    // Catastrophe arrival transitions: n -> 0 (for all n > 0)
    if (lambdaIrCatastrophe > 0) {
        for (n in 1 until Np) {
            L.set(n, 0, L.get(n, 0) + lambdaIrCatastrophe)
        }
    }

    // Batch removal arrival transitions: n -> m at rate λ * P(remove n-m) for m < n
    for ((rate, dist) in batchArrivals) {
        for (n in 1 until Np) {
            for (m in 0 until n) {
                val k = n - m  // Number of jobs to remove
                val prob = if (m > 0) {
                    // Remove exactly k jobs: P(removal = k)
                    dist.evalPMF(k.toDouble())
                } else {
                    // Remove all jobs (m = 0): P(removal >= n)
                    var cdfNMinus1 = 0.0
                    for (j in 0 until n) {
                        cdfNMinus1 += dist.evalPMF(j.toDouble())
                    }
                    1.0 - cdfNMinus1
                }
                if (prob > 0) {
                    L.set(n, m, L.get(n, m) + rate * prob)
                }
            }
        }
    }

    // Single removal negative arrival transitions: n -> n-1 (only if queue non-empty)
    if (lambdaIrNegSingle > 0) {
        for (n in 1 until Np) {
            L.set(n, n - 1, L.get(n, n - 1) + lambdaIrNegSingle)
        }
    }

    val muIr = sn.rates.get(ist, r)
    if (!muIr.isNaN() && muIr > 0) {
        // Departures to sink - use rtnodes with node indices (like MATLAB)
        // Get the node index for this station
        val nodeIdx = sn.stationToNode.get(ist).toInt()

        var probSink = 0.0
        if (sn.rtnodes != null && sn.rtnodes.numRows > 0) {
            val nNodes = sn.nodetype.size
            for (jsnk in sinkNodes) {
                for (s in 0 until K) {
                    // rtnodes indices: nodeIdx * K + classIdx (0-based)
                    val fromIdx = nodeIdx * K + r
                    val toIdx = jsnk * K + s
                    if (fromIdx < sn.rtnodes.numRows && toIdx < sn.rtnodes.numCols) {
                        probSink += sn.rtnodes.get(fromIdx, toIdx)
                    }
                }
            }
        }

        val probSelf = rt.get(ist * K + r, ist * K + r)
        val localDepartureRate = muIr * probSink

        // Departure transitions: n -> n-1 (add to existing negative arrival effects)
        for (n in 1 until Np) {
            L.set(n, n - 1, L.get(n, n - 1) + localDepartureRate)
        }

        // Self-service transitions: n -> n (diagonal, for phase transitions)
        for (n in 1 until Np) {
            L.set(n, n, muIr * probSelf)
        }
    }

    return L
}
