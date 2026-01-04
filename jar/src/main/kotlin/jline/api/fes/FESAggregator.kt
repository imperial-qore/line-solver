package jline.api.fes

import jline.api.mc.dtmc_stochcomp
import jline.api.pfqn.ld.ljd_linearize
import jline.lang.ClosedClass
import jline.lang.Network
import jline.lang.NetworkStruct
import jline.lang.constant.SchedStrategy
import jline.lang.nodes.Delay
import jline.lang.nodes.Queue
import jline.lang.nodes.Station
import jline.lang.processes.Disabled
import jline.lang.processes.Exp
import jline.solvers.mva.SolverMVA
import jline.util.Utils
import jline.util.matrix.Matrix
import java.util.*

/**
 * Flow-Equivalent Server (FES) Aggregation (Internal Implementation)
 *
 * Replaces a subset of stations in a closed product-form queueing network
 * with a single Flow-Equivalent Server having Limited Joint Dependence (LJD)
 * service rates based on throughputs in an isolated subnetwork.
 *
 * @deprecated Use jline.lang.ModelAdapter.aggregateFES() for public API
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
@Deprecated("Use ModelAdapter.aggregateFES() instead")
object FESAggregator {

    /**
     * Aggregate a station subset into a Flow-Equivalent Server
     *
     * @param model Closed product-form Network model
     * @param stationSubset List of Station objects to aggregate
     * @return FESResult containing the new model and deaggregation info
     */
    @JvmStatic
    fun aggregateFES(model: Network, stationSubset: List<Station>): FESResult {
        return aggregateFES(model, stationSubset, FESOptions.defaults())
    }

    /**
     * Aggregate a station subset into a Flow-Equivalent Server with options
     *
     * @param model Closed product-form Network model
     * @param stationSubset List of Station objects to aggregate
     * @param options FES aggregation options
     * @return FESResult containing the new model and deaggregation info
     */
    @JvmStatic
    fun aggregateFES(model: Network, stationSubset: List<Station>, options: FESOptions): FESResult {
        // Validate inputs
        validateInputs(model, stationSubset)

        // Get network structure
        val sn = model.getStruct(true)
        val M = sn.nstations
        val K = sn.nclasses

        // Get station indices for subset
        val modelStations = model.stations
        val subsetIndices = IntArray(stationSubset.size)
        for (i in stationSubset.indices) {
            for (j in 0 until M) {
                if (stationSubset[i] === modelStations[j]) {
                    subsetIndices[i] = j
                    break
                }
            }
        }

        // Get complement indices
        val allIndices = (0 until M).toList()
        val complementIndices = allIndices.filter { it !in subsetIndices.toList() }.toIntArray()

        // Set cutoffs
        val cutoffs = options.cutoffs ?: Matrix(1, K).apply {
            for (k in 0 until K) {
                set(0, k, sn.njobs.get(k).toDouble())
            }
        }

        // Compute stochastic complement routing matrices
        val rt = sn.rt

        // Build index sets for rt matrix
        val subsetRtIndices = mutableListOf<Int>()
        for (i in subsetIndices) {
            val isf = sn.stationToStateful.get(i).toInt()
            for (k in 0 until K) {
                subsetRtIndices.add((isf - 1) * K + k)
            }
        }

        val complementRtIndices = mutableListOf<Int>()
        for (i in complementIndices) {
            val isf = sn.stationToStateful.get(i).toInt()
            for (k in 0 until K) {
                complementRtIndices.add((isf - 1) * K + k)
            }
        }

        // Compute stochastic complements
        val stochCompSubset = dtmc_stochcomp(rt, subsetRtIndices.toMutableList())
        val stochCompComplement = dtmc_stochcomp(rt, complementRtIndices.toMutableList())

        // Build isolated subnetwork
        val isolatedModel = buildIsolatedModel(model, stationSubset, stochCompSubset, sn)

        // Compute throughputs for all population states
        val throughputTable = computeThroughputs(isolatedModel, cutoffs, options)

        // Create the FES model
        val fesModel = Network("${model.name}_FES")

        // Copy complement stations
        val stationMap = HashMap<Int, Station>()
        for (i in complementIndices) {
            val origStation = modelStations[i]
            val newStation = when (origStation) {
                is Queue -> {
                    val q = Queue(fesModel, origStation.name, origStation.schedStrategy)
                    if (!Utils.isInf(origStation.numberOfServers.toDouble())) {
                        q.setNumberOfServers(origStation.numberOfServers)
                    }
                    q
                }
                is Delay -> Delay(fesModel, origStation.name)
                else -> throw IllegalArgumentException("Unsupported station type: ${origStation::class.java}")
            }
            stationMap[i] = newStation
        }

        // Create FES station
        val fesStation = Queue(fesModel, "FES", SchedStrategy.PS)
        fesStation.setNumberOfServers(1)

        // Create job classes
        val refStation = if (complementIndices.isNotEmpty()) {
            stationMap[complementIndices[0]]!!
        } else {
            fesStation
        }

        val newClasses = ArrayList<ClosedClass>()
        for (k in 0 until K) {
            val origClass = model.classes[k]
            val population = sn.njobs.get(k).toInt()
            val newClass = ClosedClass(fesModel, origClass.name, population, refStation, 0)
            newClasses.add(newClass)
        }

        // Set service distributions for complement stations
        for (i in complementIndices) {
            val newStation = stationMap[i]!!
            for (k in 0 until K) {
                val rate = sn.rates.get(i, k)
                if (rate > 0 && !rate.isNaN()) {
                    newStation.setService(newClasses[k], Exp(rate))
                } else {
                    newStation.setService(newClasses[k], Disabled())
                }
            }
        }

        // Set FES service distribution with LJCD (per-class scaling)
        // For multi-class FES, use Limited Joint Class Dependence (LJCD) which allows
        // each class to have its own state-dependent scaling factor.
        // The FES service rate for class c in state (n1,...,nK) equals throughput_c(n).

        // Set base service rate Exp(1.0) for all classes; LJCD scaling provides actual rate
        for (k in 0 until K) {
            fesStation.setService(newClasses[k], Exp(1.0))
        }

        // Handle zeros in scaling tables
        val tableSize = computeTableSize(cutoffs)
        for (k in 0 until K) {
            for (idx in 0 until tableSize) {
                if (throughputTable[k].get(0, idx) < 1e-10) {
                    throughputTable[k].set(0, idx, 1e-10)
                }
            }
        }

        // Build per-class scaling map for LJCD
        val scalingMap = HashMap<jline.lang.JobClass, Matrix>()
        for (k in 0 until K) {
            scalingMap[newClasses[k]] = throughputTable[k]
        }

        // Set per-class joint dependence (LJCD)
        fesStation.setJointClassDependence(scalingMap, cutoffs)

        // Build routing matrix
        val nodes = fesModel.nodes
        val nNodes = nodes.size

        // Find FES node index
        var fesNodeIdx = 0
        for (n in 0 until nNodes) {
            if (nodes[n] === fesStation) {
                fesNodeIdx = n
                break
            }
        }

        // Build node index mapping for complement
        val complementNodeMap = HashMap<Int, Int>()
        for (i in complementIndices) {
            val station = stationMap[i]!!
            for (n in 0 until nNodes) {
                if (nodes[n] === station) {
                    complementNodeMap[i] = n
                    break
                }
            }
        }

        // Set routing probabilities using model.link()
        val P = fesModel.initRoutingMatrix()

        for (k in 0 until K) {
            val Pk = Matrix(nNodes, nNodes)

            // Routes within complement
            for (i in complementIndices) {
                if (!complementNodeMap.containsKey(i)) continue
                val iNode = complementNodeMap[i]!!
                val isf_i = sn.stationToStateful.get(i).toInt()
                val localI = complementIndices.indexOf(i)

                for (j in complementIndices) {
                    if (!complementNodeMap.containsKey(j)) continue
                    val jNode = complementNodeMap[j]!!
                    val localJ = complementIndices.indexOf(j)

                    val rtIdx_i = localI * K + k
                    val rtIdx_j = localJ * K + k

                    if (rtIdx_i < stochCompComplement.numRows && rtIdx_j < stochCompComplement.numCols) {
                        val prob = stochCompComplement.get(rtIdx_i, rtIdx_j)
                        if (prob > 1e-10) {
                            Pk.set(iNode, jNode, prob)
                        }
                    }
                }

                // Routes from complement to FES
                for (j in subsetIndices) {
                    val isf_j = sn.stationToStateful.get(j).toInt()
                    val rtIdx_i = (isf_i - 1) * K + k
                    val rtIdx_j = (isf_j - 1) * K + k

                    if (rtIdx_i < rt.numRows && rtIdx_j < rt.numCols) {
                        val prob = rt.get(rtIdx_i, rtIdx_j)
                        if (prob > 1e-10) {
                            Pk.set(iNode, fesNodeIdx, Pk.get(iNode, fesNodeIdx) + prob)
                        }
                    }
                }
            }

            // Routes from FES to complement
            for (j in complementIndices) {
                if (!complementNodeMap.containsKey(j)) continue
                val jNode = complementNodeMap[j]!!
                val isf_j = sn.stationToStateful.get(j).toInt()

                var probSum = 0.0
                for (i in subsetIndices) {
                    val isf_i = sn.stationToStateful.get(i).toInt()
                    val rtIdx_i = (isf_i - 1) * K + k
                    val rtIdx_j = (isf_j - 1) * K + k

                    if (rtIdx_i < rt.numRows && rtIdx_j < rt.numCols) {
                        probSum += rt.get(rtIdx_i, rtIdx_j)
                    }
                }

                if (probSum > 1e-10) {
                    Pk.set(fesNodeIdx, jNode, probSum / subsetIndices.size)
                }
            }

            // Normalize rows
            for (n in 0 until nNodes) {
                var rowSum = 0.0
                for (m in 0 until nNodes) {
                    rowSum += Pk.get(n, m)
                }
                if (rowSum > 1e-10) {
                    for (m in 0 until nNodes) {
                        Pk.set(n, m, Pk.get(n, m) / rowSum)
                    }
                }
            }

            P.set(newClasses[k], Pk)
        }

        // Link the model
        fesModel.link(P)

        // Build deaggregation info
        val deaggInfo = FESDeaggInfo(
            originalModel = model,
            stationSubset = stationSubset,
            subsetIndices = subsetIndices,
            complementIndices = complementIndices,
            throughputTable = throughputTable,
            cutoffs = cutoffs,
            stochCompSubset = stochCompSubset,
            stochCompComplement = stochCompComplement,
            isolatedModel = isolatedModel,
            fesNodeIdx = fesNodeIdx
        )

        return FESResult(fesModel, fesStation, deaggInfo)
    }

    /**
     * Validate inputs for FES aggregation
     */
    private fun validateInputs(model: Network, stationSubset: List<Station>) {
        if (stationSubset.isEmpty()) {
            throw IllegalArgumentException("Station subset cannot be empty")
        }

        val sn = model.getStruct(true)

        // Check model has only closed classes
        for (k in 0 until sn.nclasses) {
            if (Utils.isInf(sn.njobs.get(k))) {
                throw IllegalArgumentException("FES aggregation only applies to closed queueing networks")
            }
        }

        val modelStations = model.stations
        if (stationSubset.size >= modelStations.size) {
            throw IllegalArgumentException("Cannot aggregate all stations")
        }

        // Validate each station
        for (station in stationSubset) {
            if (station !is Queue && station !is Delay) {
                throw IllegalArgumentException("FES aggregation only supports Queue and Delay stations")
            }

            var found = false
            for (ms in modelStations) {
                if (station === ms) {
                    found = true
                    break
                }
            }
            if (!found) {
                throw IllegalArgumentException("Station ${station.name} does not belong to the model")
            }
        }
    }

    /**
     * Build isolated subnetwork from station subset
     */
    private fun buildIsolatedModel(
        model: Network,
        stationSubset: List<Station>,
        stochCompS: Matrix,
        sn: NetworkStruct
    ): Network {
        val isolatedModel = Network("${model.name}_isolated")
        val K = sn.nclasses
        val M_sub = stationSubset.size

        // Get subset station indices
        val modelStations = model.stations
        val subsetIndices = IntArray(M_sub)
        for (i in stationSubset.indices) {
            for (j in 0 until sn.nstations) {
                if (stationSubset[i] === modelStations[j]) {
                    subsetIndices[i] = j
                    break
                }
            }
        }

        // Create new stations
        val newStations = ArrayList<Station>()
        for (i in 0 until M_sub) {
            val origStation = stationSubset[i]
            val newStation = when (origStation) {
                is Queue -> {
                    val q = Queue(isolatedModel, origStation.name, origStation.schedStrategy)
                    if (!Utils.isInf(origStation.numberOfServers.toDouble())) {
                        q.setNumberOfServers(origStation.numberOfServers)
                    }
                    q
                }
                is Delay -> Delay(isolatedModel, origStation.name)
                else -> throw IllegalArgumentException("Unsupported station type")
            }
            newStations.add(newStation)
        }

        // Create closed classes with population 1 (will be changed during enumeration)
        val refStation = newStations[0]
        val newClasses = ArrayList<ClosedClass>()
        for (k in 0 until K) {
            val origClass = model.classes[k]
            val newClass = ClosedClass(isolatedModel, origClass.name, 1, refStation, 0)
            newClasses.add(newClass)
        }

        // Set service distributions
        for (i in 0 until M_sub) {
            val origIdx = subsetIndices[i]
            val newStation = newStations[i]

            for (k in 0 until K) {
                val rate = sn.rates.get(origIdx, k)
                if (rate > 0 && !rate.isNaN()) {
                    newStation.setService(newClasses[k], Exp(rate))
                } else {
                    newStation.setService(newClasses[k], Disabled())
                }
            }
        }

        // Build routing from stochastic complement
        val P = isolatedModel.initRoutingMatrix()
        val nodes = isolatedModel.nodes
        val nNodes = nodes.size

        val stationToNode = IntArray(M_sub)
        for (i in 0 until M_sub) {
            for (n in 0 until nNodes) {
                if (nodes[n] === newStations[i]) {
                    stationToNode[i] = n
                    break
                }
            }
        }

        for (k in 0 until K) {
            val Pk = Matrix(nNodes, nNodes)

            for (i in 0 until M_sub) {
                val rowIdx = i * K + k
                val iNode = stationToNode[i]

                for (j in 0 until M_sub) {
                    val colIdx = j * K + k
                    val jNode = stationToNode[j]

                    if (rowIdx < stochCompS.numRows && colIdx < stochCompS.numCols) {
                        val prob = stochCompS.get(rowIdx, colIdx)
                        if (prob > 1e-10) {
                            Pk.set(iNode, jNode, prob)
                        }
                    }
                }

                // Normalize row
                var rowSum = 0.0
                for (n in 0 until nNodes) {
                    rowSum += Pk.get(iNode, n)
                }
                if (rowSum > 1e-10) {
                    for (n in 0 until nNodes) {
                        Pk.set(iNode, n, Pk.get(iNode, n) / rowSum)
                    }
                } else {
                    // Self-loop if no routing
                    Pk.set(iNode, iNode, 1.0)
                }
            }

            P.set(newClasses[k], Pk)
        }

        isolatedModel.link(P)
        return isolatedModel
    }

    /**
     * Compute throughputs for all population states
     */
    private fun computeThroughputs(
        isolatedModel: Network,
        cutoffs: Matrix,
        options: FESOptions
    ): List<Matrix> {
        val K = isolatedModel.classes.size
        val tableSize = computeTableSize(cutoffs)

        // Initialize throughput tables
        val throughputTable = ArrayList<Matrix>()
        for (k in 0 until K) {
            throughputTable.add(Matrix(1, tableSize))
        }

        // Enumerate all population states using BFS
        val visited = HashSet<String>()
        val queue = LinkedList<IntArray>()
        queue.add(IntArray(K))

        while (queue.isNotEmpty()) {
            val nvec = queue.poll()
            val stateKey = nvec.contentToString()

            if (visited.contains(stateKey)) continue
            visited.add(stateKey)

            val totalPop = nvec.sum()

            // Convert IntArray to Matrix for ljd_linearize
            val nvecMatrix = Matrix(1, K).apply {
                for (i in 0 until K) set(0, i, nvec[i].toDouble())
            }
            val idx = ljd_linearize(nvecMatrix, cutoffs)

            if (totalPop == 0) {
                // Zero population: throughput is 0
                for (k in 0 until K) {
                    throughputTable[k].set(0, idx, 0.0)
                }
            } else {
                try {
                    // Update populations
                    for (k in 0 until K) {
                        (isolatedModel.classes[k] as ClosedClass).population = nvec[k].toDouble()
                    }

                    // Solve with MVA
                    val solver = SolverMVA(isolatedModel)
                    solver.runAnalyzer()
                    val result = solver.result

                    // Extract throughputs
                    for (k in 0 until K) {
                        if (nvec[k] > 0 && result.XN != null) {
                            throughputTable[k].set(0, idx, result.XN.get(k))
                        } else {
                            throughputTable[k].set(0, idx, 0.0)
                        }
                    }
                } catch (e: Exception) {
                    // If solver fails, set throughput to 0
                    for (k in 0 until K) {
                        throughputTable[k].set(0, idx, 0.0)
                    }
                }
            }

            // Add neighboring states
            for (k in 0 until K) {
                if (nvec[k] < cutoffs.get(0, k).toInt()) {
                    val newState = nvec.copyOf()
                    newState[k] = newState[k] + 1
                    val newKey = newState.contentToString()
                    if (!visited.contains(newKey)) {
                        queue.add(newState)
                    }
                }
            }
        }

        return throughputTable
    }

    /**
     * Compute table size from cutoffs
     */
    private fun computeTableSize(cutoffs: Matrix): Int {
        var size = 1
        for (k in 0 until cutoffs.numCols) {
            size *= (cutoffs.get(0, k).toInt() + 1)
        }
        return size
    }
}
