/**
 * @file MAPQN Parameters Factory
 * 
 * Factory class for creating parameter objects for MAP queueing network
 * (MAPQN) analysis. Handles the instantiation and configuration of various
 * parameter types used in MAP-based queueing network modeling.
 * 
 * @since LINE 3.0
 */
package jline.api.mapqn

import jline.GlobalConstants.Inf
import jline.lang.NetworkStruct
import jline.lang.constant.NodeType
import jline.lang.constant.ProcessType
import jline.util.matrix.Matrix
import java.util.Arrays

/**
 * Factory class for creating Mapqn_parameters from NetworkStruct.
 * This factory analyzes the network structure and creates the appropriate
 * Mapqn_parameters subclass based on the network characteristics.
 */
object Mapqn_parameters_factory {
    
    /**
     * Creates appropriate Mapqn_parameters based on NetworkStruct properties.
     * 
     * @param networkStruct The network structure to convert
     * @return Appropriate Mapqn_parameters instance
     * @throws IllegalArgumentException if the network structure is not supported
     */
    fun createFromNetworkStruct(networkStruct: NetworkStruct): Mapqn_parameters {
        // Validate network structure
        networkStruct.validateStructuralConsistency()
        
        // Extract basic parameters
        val M = networkStruct.nstations
        val N = networkStruct.nclosedjobs
        
        if (M <= 0) {
            throw IllegalArgumentException("Network must have at least one station")
        }
        if (N <= 0) {
            throw IllegalArgumentException("Network must have closed jobs for MAPQN analysis")
        }
        
        // Determine the type of network and create appropriate parameters
        return when {
            hasFiniteCapacity(networkStruct) -> createFiniteCapacityParameters(networkStruct)
            hasMultiplePhases(networkStruct) -> createLinearReductionParameters(networkStruct)
            else -> createBasicParameters(networkStruct)
        }
    }
    
    /**
     * Creates LinearReductionParameters from NetworkStruct.
     * This constructor is used for multi-phase networks without finite capacity.
     */
    fun createLinearReductionParameters(networkStruct: NetworkStruct): LinearReductionParameters {
        val M = networkStruct.nstations
        val N = networkStruct.nclosedjobs
        
        // Extract number of phases for each station
        val K = IntArray(M)
        for (i in 0 until M) {
            // Get station index in nodes array
            val stationIdx = networkStruct.stationToNode?.get(i, 0)?.toInt() ?: i
            // Phases matrix should have nstations rows
            K[i] = networkStruct.phases?.get(i, 0)?.toInt() ?: 1
        }
        
        // Extract service rates
        val mu = Array(M) { i ->
            val stationObj = networkStruct.stations?.get(i)
            val muMap = networkStruct.mu?.get(stationObj)
            
            if (muMap != null && muMap.isNotEmpty()) {
                // Get the first job class service rates
                val firstClass = networkStruct.jobclasses?.firstOrNull()
                val muMatrix = muMap[firstClass]
                
                if (muMatrix != null) {
                    // Convert to phase transition matrix
                    convertToPhaseTransitionMatrix(muMatrix, K[i])
                } else {
                    // Default: single phase with unit rate
                    Matrix.eye(K[i])
                }
            } else {
                // Default: single phase with unit rate
                Matrix.eye(K[i])
            }
        }
        
        // Extract routing probabilities
        val r = extractRoutingMatrix(networkStruct)
        
        // Background transition rates (typically zero for basic networks)
        val v = Array(M) { i ->
            Matrix.zeros(K[i], K[i])
        }
        
        return LinearReductionParameters(M, N, K, mu, r, v)
    }
    
    /**
     * Creates parameters for networks with finite capacity queues.
     * Determines whether to use Mapqn_qr_bounds_bas or Mapqn_qr_bounds_rsrd based on blocking characteristics.
     */
    private fun createFiniteCapacityParameters(networkStruct: NetworkStruct): Mapqn_parameters {
        val M = networkStruct.nstations
        val N = networkStruct.nclosedjobs
        
        // Extract capacities
        val F = IntArray(M)
        for (i in 0 until M) {
            F[i] = networkStruct.cap?.get(i, 0)?.toInt() ?: Int.MAX_VALUE
        }
        
        // Check if we have blocking configurations (would indicate BAS model)
        // For now, we'll default to RSRD model as it's more general
        return createMapqn_qr_bounds_rsrd_parameters(networkStruct, F)
    }
    
    /**
     * Creates Mapqn_qr_bounds_rsrd_parameters from NetworkStruct.
     */
    private fun createMapqn_qr_bounds_rsrd_parameters(networkStruct: NetworkStruct, F: IntArray): Mapqn_qr_bounds_rsrd_parameters {
        val M = networkStruct.nstations
        val N = networkStruct.nclosedjobs
        
        // Extract number of phases
        val K = IntArray(M)
        for (i in 0 until M) {
            K[i] = networkStruct.phases?.get(i, 0)?.toInt() ?: 1
        }
        
        // Extract service rates
        val mu = Array(M) { i ->
            extractServiceRateMatrix(networkStruct, i, K[i])
        }
        
        // Background transition rates
        val v = Array(M) { i ->
            Matrix.zeros(K[i], K[i])
        }
        
        // Load-dependent rates (default to 1.0 for all loads)
        val alpha = Array(M) { i ->
            DoubleArray(N) { 1.0 }
        }
        
        // Extract load-dependent scaling if available
        if (networkStruct.cdscaling != null && networkStruct.stations != null) {
            for (i in 0 until M) {
                val stationObj = networkStruct.stations[i]
                val scalingFunc = networkStruct.cdscaling[stationObj]
                if (scalingFunc != null) {
                    // Apply scaling function for each population level
                    for (n in 1..N) {
                        val stateMatrix = Matrix(1, 1)
                        stateMatrix.set(0, 0, n.toDouble())
                        alpha[i][n-1] = scalingFunc.apply(stateMatrix)
                    }
                }
            }
        }
        
        // Extract routing probabilities
        val r = extractRoutingMatrix(networkStruct)
        
        return Mapqn_qr_bounds_rsrd_parameters(M, N, F, K, mu, v, alpha, r)
    }
    
    /**
     * Creates basic parameters for single-phase networks.
     * This creates LinearReductionParameters with K[i]=1 for all stations.
     */
    private fun createBasicParameters(networkStruct: NetworkStruct): LinearReductionParameters {
        val M = networkStruct.nstations
        val N = networkStruct.nclosedjobs
        
        // All stations have single phase
        val K = IntArray(M) { 1 }
        
        // Extract service rates
        val mu = Array(M) { i ->
            val rate = networkStruct.rates?.get(i, 0) ?: 1.0
            Matrix(arrayOf(doubleArrayOf(rate)))
        }
        
        // Extract routing probabilities
        val r = extractRoutingMatrix(networkStruct)
        
        // No background transitions
        val v = Array(M) { Matrix(arrayOf(doubleArrayOf(0.0))) }
        
        return LinearReductionParameters(M, N, K, mu, r, v)
    }
    
    /**
     * Checks if the network has finite capacity queues.
     */
    private fun hasFiniteCapacity(networkStruct: NetworkStruct): Boolean {
        if (networkStruct.cap == null) return false
        
        for (i in 0 until networkStruct.nstations) {
            val capacity = networkStruct.cap.get(i, 0)
            if (capacity > 0 && capacity < Inf) {
                return true
            }
        }
        return false
    }
    
    /**
     * Checks if the network has multi-phase service.
     */
    private fun hasMultiplePhases(networkStruct: NetworkStruct): Boolean {
        if (networkStruct.phases == null) return false
        
        for (i in 0 until networkStruct.nstations) {
            if (networkStruct.phases.get(i, 0) > 1) {
                return true
            }
        }
        return false
    }
    
    /**
     * Extracts the routing probability matrix from NetworkStruct.
     * Returns an MxM matrix where entry (i,j) is the probability of routing from station i to station j.
     */
    private fun extractRoutingMatrix(networkStruct: NetworkStruct): Matrix {
        val M = networkStruct.nstations
        
        // Check if we have a routing table
        if (networkStruct.rt != null) {
            // Extract station-to-station routing from the full routing table
            // rt is indexed by (station*class), we need to extract station-level routing
            val nclasses = networkStruct.nclasses
            val r = Matrix(M, M)
            
            for (i in 0 until M) {
                for (j in 0 until M) {
                    // Sum routing probabilities across all classes
                    var sum = 0.0
                    for (c in 0 until nclasses) {
                        val fromIdx = i * nclasses + c
                        val toIdx = j * nclasses + c  // Assuming no class switching for now
                        sum += networkStruct.rt.get(fromIdx, toIdx)
                    }
                    r.set(i, j, sum / nclasses)  // Average across classes
                }
            }
            
            // Normalize rows to ensure stochastic matrix
            for (i in 0 until M) {
                val rowSum = (0 until M).sumOf { j -> r.get(i, j) }
                if (rowSum > 0) {
                    for (j in 0 until M) {
                        r.set(i, j, r.get(i, j) / rowSum)
                    }
                }
            }
            
            return r
        }
        
        // Fallback: create default routing based on network topology
        return createDefaultRoutingMatrix(networkStruct)
    }
    
    /**
     * Creates a default routing matrix based on network topology.
     */
    private fun createDefaultRoutingMatrix(networkStruct: NetworkStruct): Matrix {
        val M = networkStruct.nstations
        val r = Matrix(M, M)
        
        // Simple cyclic routing as default
        for (i in 0 until M) {
            for (j in 0 until M) {
                if (j == (i + 1) % M) {
                    r.set(i, j, 1.0)
                } else {
                    r.set(i, j, 0.0)
                }
            }
        }
        
        return r
    }
    
    /**
     * Extracts service rate matrix for a specific station.
     */
    private fun extractServiceRateMatrix(networkStruct: NetworkStruct, stationIdx: Int, numPhases: Int): Matrix {
        // Try to get from mu map first
        if (networkStruct.mu != null && networkStruct.stations != null) {
            val station = networkStruct.stations[stationIdx]
            val muMap = networkStruct.mu[station]
            
            if (muMap != null && muMap.isNotEmpty()) {
                // Get rates for first job class
                val firstClass = networkStruct.jobclasses?.firstOrNull()
                val muMatrix = muMap[firstClass]
                
                if (muMatrix != null) {
                    return convertToPhaseTransitionMatrix(muMatrix, numPhases)
                }
            }
        }
        
        // Fallback to rates matrix
        if (networkStruct.rates != null) {
            val rate = networkStruct.rates.get(stationIdx, 0)
            if (rate > 0) {
                // Create diagonal matrix with service rate
                val mu = Matrix.zeros(numPhases, numPhases)
                for (k in 0 until numPhases) {
                    mu.set(k, k, rate)
                }
                return mu
            }
        }
        
        // Default: identity matrix (unit rates)
        return Matrix.eye(numPhases)
    }
    
    /**
     * Converts a service rate specification to a phase transition matrix.
     */
    private fun convertToPhaseTransitionMatrix(serviceRates: Matrix, numPhases: Int): Matrix {
        // If already the right size, return as-is
        if (serviceRates.numRows == numPhases && serviceRates.numCols == numPhases) {
            return serviceRates.copy()
        }
        
        // If single value, create diagonal matrix
        if (serviceRates.numRows == 1 && serviceRates.numCols == 1) {
            val rate = serviceRates.get(0, 0)
            val mu = Matrix.zeros(numPhases, numPhases)
            for (k in 0 until numPhases) {
                mu.set(k, k, rate)
            }
            return mu
        }
        
        // If vector of phase rates, create diagonal matrix
        if ((serviceRates.numRows == numPhases && serviceRates.numCols == 1) || 
            (serviceRates.numRows == 1 && serviceRates.numCols == numPhases)) {
            val mu = Matrix.zeros(numPhases, numPhases)
            for (k in 0 until numPhases) {
                val rate = if (serviceRates.numRows == 1) {
                    serviceRates.get(0, k)
                } else {
                    serviceRates.get(k, 0)
                }
                mu.set(k, k, rate)
            }
            return mu
        }
        
        // Default: identity matrix
        return Matrix.eye(numPhases)
    }
}