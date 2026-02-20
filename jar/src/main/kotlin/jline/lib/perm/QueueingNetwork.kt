package jline.lib.perm

import jline.GlobalConstants.Inf
import jline.util.matrix.Matrix
import kotlin.math.pow

/**
 * Utility function to make a matrix doubly stochastic.
 * 
 * Uses the Sinkhorn algorithm to iteratively normalize rows and columns
 * until the matrix becomes approximately doubly stochastic.
 * 
 * @param M The input matrix
 * @return Pair of (doubly stochastic matrix, rescaling factor)
 */
fun preprocessingDS(M: Matrix): Pair<Matrix, Double> {
    val n = M.numRows
    val minValue = 2.220446049250314e-16
    
    // Ensure minimum values
    val result = Matrix(Array(n) { i ->
        DoubleArray(n) { j -> maxOf(M.get(i, j), minValue) }
    })
    
    val X = Matrix.eye(n)
    val Y = Matrix.eye(n)
    
    var maxRowError = Inf
    var maxColError = Inf
    val tolerance = 0.001
    
    while (maxRowError > tolerance || maxColError > tolerance) {
        // Normalize columns
        val colSums = DoubleArray(n) { j ->
            (0 until n).sumOf { i -> result.get(i, j) }
        }
        
        for (j in 0 until n) {
            if (colSums[j] > 0) {
                for (i in 0 until n) {
                    result.set(i, j, result.get(i, j) / colSums[j])
                    Y.set(i, j, Y.get(i, j) / colSums[j])
                }
            }
        }
        
        // Normalize rows
        val rowSums = DoubleArray(n) { i ->
            (0 until n).sumOf { j -> result.get(i, j) }
        }
        
        for (i in 0 until n) {
            if (rowSums[i] > 0) {
                for (j in 0 until n) {
                    result.set(i, j, result.get(i, j) / rowSums[i])
                    X.set(i, j, X.get(i, j) / rowSums[i])
                }
            }
        }
        
        // Check convergence
        val newColSums = DoubleArray(n) { j ->
            (0 until n).sumOf { i -> result.get(i, j) }
        }
        val newRowSums = DoubleArray(n) { i ->
            (0 until n).sumOf { j -> result.get(i, j) }
        }
        
        maxColError = newColSums.maxOf { kotlin.math.abs(it - 1.0) }
        maxRowError = newRowSums.maxOf { kotlin.math.abs(it - 1.0) }
    }
    
    // Compute rescaling factor
    var rescalingFactor = 1.0
    for (i in 0 until n) {
        rescalingFactor *= X.get(i, i) * Y.get(i, i)
    }
    
    return Pair(result, rescalingFactor)
}

/**
 * Queueing network model without think time.
 * 
 * This class models a closed queueing network where jobs circulate between
 * queues without any external think time. The network is characterized by:
 * - Number of queues and job classes
 * - Population of each class
 * - Mean service demands at each queue for each class
 * 
 * The permanent of matrices constructed from service demands gives the
 * unnormalized probability of network states.
 * 
 * @param numberOfQueues Number of service stations in the network
 * @param numberOfClasses Number of job classes
 * @param numberPerClass Population of each job class
 * @param meanServiceDemand Matrix of mean service demands [queue][class]
 * @param progress Whether to show progress during computations (default: true)
 */
class NetworkNoThink(
    private val numberOfQueues: Int,
    private val numberOfClasses: Int,
    private val numberPerClass: IntArray,
    private val meanServiceDemand: Array<DoubleArray>,
    private val progress: Boolean = true
) {
    
    /**
     * Compute the unnormalized probability of a joint state.
     * 
     * @param state Matrix representing the state [queue][class] = number of jobs
     * @return Unnormalized probability of the state
     */
    fun joint(state: Array<IntArray>): Double {
        var termMSD = 1.0
        var invFactorial = 1.0
        var factorialTerm = 1.0
        
        // Compute MSD^(n_m) term
        for (i in 0 until numberOfQueues) {
            for (j in 0 until numberOfClasses) {
                termMSD *= meanServiceDemand[i][j].pow(state[i][j])
            }
        }
        
        // Compute factorial terms
        for (i in 0 until numberOfQueues) {
            for (j in 0 until numberOfClasses) {
                invFactorial *= factorial(state[i][j])
            }
            
            val queueTotal = state[i].sum()
            factorialTerm *= factorial(queueTotal)
        }
        
        return termMSD * factorialTerm / invFactorial
    }
    
    /**
     * Compute marginal probability of a state using specified solver.
     * 
     * @param solver The permanent solver to use
     * @param state Vector representing marginal state [queue] = total jobs
     * @param preprocessing Whether to apply doubly stochastic preprocessing
     * @return Triple of (probability, computation time, memory usage)
     */
    fun marginal(solver: PermSolver, state: IntArray, preprocessing: Boolean = false): Triple<Double, Long, Long> {
        // Construct matrix from service demands repeated according to state
        val matrixData = Array(state.sum()) { DoubleArray(state.sum()) }
        
        var rowIndex = 0
        for (i in 0 until numberOfQueues) {
            repeat(state[i]) {
                var colIndex = 0
                for (j in 0 until numberOfQueues) {
                    repeat(state[j]) {
                        matrixData[rowIndex][colIndex] = meanServiceDemand[i][j]
                        colIndex++
                    }
                }
                rowIndex++
            }
        }
        
        val matrix = Matrix(matrixData)
        
        // Apply preprocessing if requested
        val (processedMatrix, rescalingFactor) = if (preprocessing) {
            preprocessingDS(matrix)
        } else {
            Pair(matrix, 1.0)
        }
        
        // Create solver instance with processed matrix
        val solverInstance = when (solver) {
            is BethePermanent -> BethePermanent(processedMatrix, solve = true)
            is NaivePermanent -> NaivePermanent(processedMatrix, solve = true)
            is RyzerPermanent -> RyzerPermanent(processedMatrix, solve = true)
            is AdaPartSampler -> AdaPartSampler(processedMatrix, solve = true)
            is HuberLawSampler -> HuberLawSampler(processedMatrix, solve = true)
            else -> throw IllegalArgumentException("Unsupported solver type")
        }
        
        // Calculate factorial normalization
        var factorialNormalization = 1.0
        for (jobs in state) {
            factorialNormalization *= factorial(jobs)
        }
        
        val probability = solverInstance.value * rescalingFactor / factorialNormalization
        
        return Triple(probability, solverInstance.time, solverInstance.memory)
    }
    
    /**
     * Generate dictionary of all marginal states with their probabilities.
     * 
     * @param solver The permanent solver to use
     * @param preprocessing Whether to apply doubly stochastic preprocessing
     * @return Map of state vectors to (probability, time, memory) triples
     */
    fun generateMarginal(solver: PermSolver, preprocessing: Boolean = false): Map<IntArray, Triple<Double, Long, Long>> {
        val results = mutableMapOf<IntArray, Triple<Double, Long, Long>>()
        
        // Generate all possible marginal states
        val states = generateAllMarginalStates()
        
        for (state in states) {
            val result = marginal(solver, state, preprocessing)
            results[state] = result
            
            if (progress) {
                println("Processed state: ${state.contentToString()}")
            }
        }
        
        return results
    }
    
    /**
     * Generate all possible marginal states (total jobs per queue).
     */
    private fun generateAllMarginalStates(): List<IntArray> {
        val totalJobs = numberPerClass.sum()
        val states = mutableListOf<IntArray>()
        
        // Generate all ways to distribute totalJobs among numberOfQueues queues
        generateStatesRecursive(IntArray(numberOfQueues), 0, totalJobs, states)
        
        return states
    }
    
    /**
     * Recursively generate all state combinations.
     */
    private fun generateStatesRecursive(
        currentState: IntArray,
        queueIndex: Int,
        remainingJobs: Int,
        results: MutableList<IntArray>
    ) {
        if (queueIndex == numberOfQueues - 1) {
            currentState[queueIndex] = remainingJobs
            results.add(currentState.copyOf())
            return
        }
        
        for (jobs in 0..remainingJobs) {
            currentState[queueIndex] = jobs
            generateStatesRecursive(currentState, queueIndex + 1, remainingJobs - jobs, results)
        }
    }
    
    /**
     * Simple factorial computation.
     */
    private fun factorial(n: Int): Double {
        if (n <= 1) return 1.0
        var result = 1.0
        for (i in 2..n) {
            result *= i
        }
        return result
    }
}

/**
 * Queueing network model with think time.
 * 
 * This class extends the basic queueing network model to include think time,
 * representing time spent by jobs outside the queueing network (e.g., user
 * think time in interactive systems).
 * 
 * @param numberOfQueues Number of service stations in the network
 * @param numberOfClasses Number of job classes
 * @param numberPerClass Population of each job class
 * @param thinkTime Think time for each job class
 * @param meanServiceDemand Matrix of mean service demands [queue][class]
 * @param progress Whether to show progress during computations (default: true)
 */
class NetworkThink(
    private val numberOfQueues: Int,
    private val numberOfClasses: Int,
    private val numberPerClass: IntArray,
    private val thinkTime: DoubleArray,
    private val meanServiceDemand: Array<DoubleArray>,
    private val progress: Boolean = true
) {
    
    /**
     * Compute the unnormalized probability of a joint state including think time effects.
     * 
     * @param state Matrix representing the state [queue][class] = number of jobs
     * @return Unnormalized probability of the state
     */
    fun joint(state: Array<IntArray>): Double {
        var termMSD = 1.0
        var invFactorial = 1.0
        var factorialTerm = 1.0
        var thinkTimeTerm = 1.0
        
        // Compute MSD^(n_m) term
        for (i in 0 until numberOfQueues) {
            for (j in 0 until numberOfClasses) {
                termMSD *= meanServiceDemand[i][j].pow(state[i][j])
            }
        }
        
        // Compute factorial terms
        for (i in 0 until numberOfQueues) {
            for (j in 0 until numberOfClasses) {
                invFactorial *= factorial(state[i][j])
            }
            
            val queueTotal = state[i].sum()
            factorialTerm *= factorial(queueTotal)
        }
        
        // Compute think time term
        for (j in 0 until numberOfClasses) {
            val totalInClass = state.sumOf { it[j] }
            val thinkingJobs = numberPerClass[j] - totalInClass
            thinkTimeTerm *= thinkTime[j].pow(thinkingJobs)
        }
        
        return termMSD * factorialTerm * thinkTimeTerm / invFactorial
    }
    
    /**
     * Compute marginal probability including think time effects.
     * 
     * @param solver The permanent solver to use
     * @param state Vector representing marginal state [queue] = total jobs
     * @param preprocessing Whether to apply doubly stochastic preprocessing
     * @return Triple of (probability, computation time, memory usage)
     */
    fun marginal(solver: PermSolver, state: IntArray, preprocessing: Boolean = false): Triple<Double, Long, Long> {
        // Extended matrix includes think time stations
        val totalStations = numberOfQueues + numberOfClasses
        val totalJobs = state.sum()
        
        val matrixData = Array(totalJobs) { DoubleArray(totalJobs) }
        
        var rowIndex = 0
        
        // Service stations
        for (i in 0 until numberOfQueues) {
            repeat(state[i]) {
                var colIndex = 0
                
                // Demands at service stations
                for (j in 0 until numberOfQueues) {
                    repeat(state[j]) {
                        matrixData[rowIndex][colIndex] = meanServiceDemand[i][j]
                        colIndex++
                    }
                }
                
                // Think time "demands" (jobs in think stations)
                for (j in 0 until numberOfClasses) {
                    val thinkingJobs = numberPerClass[j] - state.sum() // Simplified - assumes single class
                    repeat(thinkingJobs) {
                        matrixData[rowIndex][colIndex] = thinkTime[j]
                        colIndex++
                    }
                }
                
                rowIndex++
            }
        }
        
        val matrix = Matrix(matrixData)
        
        // Apply preprocessing if requested
        val (processedMatrix, rescalingFactor) = if (preprocessing) {
            preprocessingDS(matrix)
        } else {
            Pair(matrix, 1.0)
        }
        
        // Create solver instance
        val solverInstance = when (solver) {
            is BethePermanent -> BethePermanent(processedMatrix, solve = true)
            is NaivePermanent -> NaivePermanent(processedMatrix, solve = true)
            is RyzerPermanent -> RyzerPermanent(processedMatrix, solve = true)
            is AdaPartSampler -> AdaPartSampler(processedMatrix, solve = true)
            is HuberLawSampler -> HuberLawSampler(processedMatrix, solve = true)
            else -> throw IllegalArgumentException("Unsupported solver type")
        }
        
        // Calculate factorial normalization
        var factorialNormalization = 1.0
        for (jobs in state) {
            factorialNormalization *= factorial(jobs)
        }
        
        val probability = solverInstance.value * rescalingFactor / factorialNormalization
        
        return Triple(probability, solverInstance.time, solverInstance.memory)
    }
    
    /**
     * Generate dictionary of all marginal states with their probabilities.
     * 
     * @param solver The permanent solver to use
     * @param preprocessing Whether to apply doubly stochastic preprocessing
     * @return Map of state vectors to (probability, time, memory) triples
     */
    fun generateMarginal(solver: PermSolver, preprocessing: Boolean = false): Map<IntArray, Triple<Double, Long, Long>> {
        val results = mutableMapOf<IntArray, Triple<Double, Long, Long>>()
        
        // Generate all possible marginal states
        val states = generateAllMarginalStates()
        
        for (state in states) {
            val result = marginal(solver, state, preprocessing)
            results[state] = result
            
            if (progress) {
                println("Processed state: ${state.contentToString()}")
            }
        }
        
        return results
    }
    
    /**
     * Generate all possible marginal states.
     */
    private fun generateAllMarginalStates(): List<IntArray> {
        val totalJobs = numberPerClass.sum()
        val states = mutableListOf<IntArray>()
        
        generateStatesRecursive(IntArray(numberOfQueues), 0, totalJobs, states)
        
        return states
    }
    
    /**
     * Recursively generate all state combinations.
     */
    private fun generateStatesRecursive(
        currentState: IntArray,
        queueIndex: Int,
        remainingJobs: Int,
        results: MutableList<IntArray>
    ) {
        if (queueIndex == numberOfQueues - 1) {
            currentState[queueIndex] = remainingJobs
            results.add(currentState.copyOf())
            return
        }
        
        for (jobs in 0..remainingJobs) {
            currentState[queueIndex] = jobs
            generateStatesRecursive(currentState, queueIndex + 1, remainingJobs - jobs, results)
        }
    }
    
    /**
     * Simple factorial computation.
     */
    private fun factorial(n: Int): Double {
        if (n <= 1) return 1.0
        var result = 1.0
        for (i in 2..n) {
            result *= i
        }
        return result
    }
}