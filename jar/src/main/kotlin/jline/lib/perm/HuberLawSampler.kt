package jline.lib.perm

import jline.GlobalConstants.Inf
import jline.util.matrix.Matrix
import kotlin.math.*
import kotlin.random.Random

/**
 * Sampling method to approximate the permanent using the Huber-Law bound.
 * 
 * The Huber-Law algorithm uses acceptance-rejection sampling with a carefully constructed
 * bound to estimate the permanent. The method rescales the matrix to be doubly stochastic
 * and then uses the Huber-Law bound for efficient sampling.
 * 
 * Key features:
 * - Provides theoretical guarantees on approximation quality
 * - Three sampling modes: classic (with acceptance guarantee), time-limited, and sample-limited
 * - Automatic matrix rescaling to doubly stochastic form
 * 
 * @param matrix The matrix for which to compute the permanent approximation
 * @param delta Approximation parameter: result will be (1+delta) times the true value with high probability (default: 0.1)
 * @param alpha2 Accuracy of doubly stochastic rescaling (default: 0.000001)
 * @param epsilon Confidence parameter: approximation holds with probability (1-epsilon) (default: 0.1)
 * @param mode Sampling mode: "classic", "time", or "sample" (default: "classic")
 * @param numberOfSamples Maximum number of samples for sample mode (default: 1000)
 * @param maximumTime Maximum time in seconds for time mode (default: 30)
 * @param solve Whether to automatically run solve() after construction (default: false)
 */
class HuberLawSampler(
    matrix: Matrix,
    private val delta: Double = 0.1,
    private val alpha2: Double = 0.000001,
    private val epsilon: Double = 0.1,
    private val mode: String = "classic",
    private val numberOfSamples: Int = 1000,
    private val maximumTime: Long = 30000, // milliseconds
    solve: Boolean = false
) : PermSolver(matrix) {
    
    private val random = Random.Default
    private lateinit var C: Matrix // Rescaled doubly stochastic matrix
    private var rescalingConstant: Double = 1.0
    
    // Tracking variables for analysis
    var sampleAccepted: List<Int> = emptyList()
        private set
    var sampleTime: List<Long> = emptyList()
        private set
    var permStep: List<Double> = emptyList()
        private set
    
    init {
        if (solve) {
            solve()
        }
    }
    
    override fun compute() {
        value = when (mode) {
            "classic" -> samplingClassic()
            "time" -> samplingTime()
            "sample" -> samplingSample()
            else -> samplingClassic()
        }
    }
    
    /**
     * Classic Huber-Law sampling with theoretical guarantees.
     */
    private fun samplingClassic(): Double {
        rescale()
        
        // Number of accepted samples needed for theoretical guarantee
        val K = (14 * delta.pow(-2) * ln(2 / epsilon)).toInt()
        
        val startTime = System.currentTimeMillis()
        val acceptedList = mutableListOf<Int>()
        val timeList = mutableListOf<Long>()
        
        var acceptedCount = 0
        
        while (acceptedCount < K) {
            val sigma = sample()
            val isAccepted = if (sigma.any { it == n }) 0 else 1
            
            acceptedList.add(isAccepted)
            timeList.add(System.currentTimeMillis() - startTime)
            
            if (isAccepted == 1) {
                acceptedCount++
            }
        }
        
        sampleAccepted = acceptedList
        sampleTime = timeList
        computePermStep()
        
        return sampleAccepted.sum().toDouble() / sampleAccepted.size * rescalingConstant
    }
    
    /**
     * Time-constrained Huber-Law sampling.
     */
    private fun samplingTime(): Double {
        rescale()
        
        val startTime = System.currentTimeMillis()
        val acceptedList = mutableListOf<Int>()
        val timeList = mutableListOf<Long>()
        
        while (System.currentTimeMillis() - startTime < maximumTime) {
            val sigma = sample()
            val isAccepted = if (sigma.any { it == n }) 0 else 1
            
            acceptedList.add(isAccepted)
            timeList.add(System.currentTimeMillis() - startTime)
        }
        
        sampleAccepted = acceptedList
        sampleTime = timeList
        computePermStep()
        
        return if (sampleAccepted.isNotEmpty()) {
            sampleAccepted.sum().toDouble() / sampleAccepted.size * rescalingConstant
        } else 0.0
    }
    
    /**
     * Sample-count-constrained Huber-Law sampling.
     */
    private fun samplingSample(): Double {
        rescale()
        
        val startTime = System.currentTimeMillis()
        val acceptedList = mutableListOf<Int>()
        val timeList = mutableListOf<Long>()
        
        while (acceptedList.size < numberOfSamples) {
            val sigma = sample()
            val isAccepted = if (sigma.any { it == n }) 0 else 1
            
            acceptedList.add(isAccepted)
            timeList.add(System.currentTimeMillis() - startTime)
        }
        
        sampleAccepted = acceptedList
        sampleTime = timeList
        computePermStep()
        
        return sampleAccepted.sum().toDouble() / sampleAccepted.size * rescalingConstant
    }
    
    /**
     * Generate one sample using Huber-Law method.
     * 
     * @return Array representing a permutation, or array of n's if rejected
     */
    private fun sample(): IntArray {
        val M = C.copy()
        val sigma = IntArray(n)
        
        for (j in 0 until n) {
            // Compute upper bound
            val rowSums = DoubleArray(n) { i ->
                (0 until n).sumOf { k -> M.get(i, k) }
            }
            val ub = rowSums.map { h(it) }.reduce { acc, v -> acc * v } / exp(n.toDouble())
            
            // Compute pseudo probabilities
            val p = precomputing(M, j)
            
            // Normalize by upper bound
            val normalizedP = p.map { it / ub }
            
            // Add slack term
            val probSum = normalizedP.sum()
            val prob = normalizedP.toMutableList()
            prob.add(1.0 - probSum)
            
            // Correct negative probabilities
            if (prob.last() < 0) {
                val posSum = prob.dropLast(1).sum()
                for (i in 0 until prob.size - 1) {
                    prob[i] = prob[i] / posSum
                }
                prob[prob.size - 1] = 0.0
            }
            
            // Random choice
            val randVal = random.nextDouble()
            var cumSum = 0.0
            var selectedI = n // Default to slack (rejection)
            
            for (i in prob.indices) {
                cumSum += prob[i]
                if (randVal <= cumSum) {
                    selectedI = i
                    break
                }
            }
            
            // If slack selected, return rejection
            if (selectedI == n) {
                return IntArray(n) { n }
            }
            
            // Add value to permutation
            sigma[j] = selectedI
            
            // Modify current matrix (zero out row selectedI and column j, except intersection)
            val newMatrix = Matrix.zeros(n, n)
            for (row in 0 until n) {
                for (col in 0 until n) {
                    if ((row == selectedI && col == j) || (row != selectedI && col != j)) {
                        newMatrix.set(row, col, M.get(row, col))
                    }
                }
            }
            
            // Update M for next iteration
            for (row in 0 until n) {
                for (col in 0 until n) {
                    M.set(row, col, newMatrix.get(row, col))
                }
            }
        }
        
        return sigma
    }
    
    /**
     * Compute part of the Huber-Law bound.
     */
    private fun h(r: Double): Double {
        return if (r >= 1.0) {
            r + 0.5 * ln(maxOf(r, 1.0)) + E - 1
        } else {
            1 + (E - 1) * r
        }
    }
    
    /**
     * Compute pseudo probabilities for column j.
     */
    private fun precomputing(M: Matrix, j: Int): DoubleArray {
        val c = DoubleArray(n) { i -> M.get(i, j) }
        val r = DoubleArray(n) { i ->
            (0 until n).sumOf { k -> M.get(i, k) } - c[i]
        }
        
        val hr = r.map { h(it) }
        val hrProduct = hr.reduce { acc, v -> acc * v }
        
        return DoubleArray(n) { i ->
            hrProduct / hr[i] * c[i] / exp(n - 1.0)
        }
    }
    
    /**
     * Rescale the matrix to doubly stochastic form.
     */
    private fun rescale() {
        // Solve maximum assignment problem for scaling constant
        val logMatrix = Matrix(Array(n) { i ->
            DoubleArray(n) { j -> ln(maxOf(matrix.get(i, j), 1e-10)) }
        })
        
        // Simple greedy assignment for approximation
        val assignment = hungarianAssignment(logMatrix)
        var alpha3 = 1.0
        for (i in assignment.indices) {
            alpha3 *= matrix.get(i, assignment[i])
        }
        
        // Scale by maximum element
        val maxElement = (0 until n).maxOf { i ->
            (0 until n).maxOf { j -> matrix.get(i, j) }
        }
        
        val MScaled = Matrix(Array(n) { i ->
            DoubleArray(n) { j -> matrix.get(i, j) / maxElement }
        })
        
        // Minimum threshold
        val alpha1 = alpha3 * delta / 3 / factorial(n)
        
        // Replace small values
        for (i in 0 until n) {
            for (j in 0 until n) {
                if (MScaled.get(i, j) < alpha1) {
                    MScaled.set(i, j, alpha1)
                }
            }
        }
        
        // Make doubly stochastic using Sinkhorn algorithm
        val (doublyStochastic, X, Y) = makeDoublyStochastic(MScaled)
        
        // Normalize by maximum of each row
        val Z = Matrix.eye(n)
        for (i in 0 until n) {
            val maxVal = (0 until n).maxOf { j -> doublyStochastic.get(i, j) }
            Z.set(i, i, 1.0 / maxVal)
        }
        
        C = Z.mult(doublyStochastic)
        
        // Compute rescaling constant
        val rowSums = DoubleArray(n) { i ->
            (0 until n).sumOf { j -> C.get(i, j) }
        }
        val hProduct = rowSums.map { h(it) / E }.reduce { acc, v -> acc * v }
        
        var diagonalProduct = 1.0
        for (i in 0 until n) {
            diagonalProduct *= X.get(i, i) * Y.get(i, i) * Z.get(i, i)
        }
        
        rescalingConstant = hProduct / diagonalProduct * maxElement.pow(n)
    }
    
    /**
     * Simple factorial computation.
     */
    private fun factorial(n: Int): Double {
        var result = 1.0
        for (i in 1..n) {
            result *= i
        }
        return result
    }
    
    /**
     * Simple greedy Hungarian assignment approximation.
     */
    private fun hungarianAssignment(costMatrix: Matrix): IntArray {
        val assignment = IntArray(n) { -1 }
        val usedCols = BooleanArray(n)
        
        for (i in 0 until n) {
            var bestCol = -1
            var bestCost = Inf
            
            for (j in 0 until n) {
                if (!usedCols[j] && -costMatrix.get(i, j) < bestCost) {
                    bestCost = -costMatrix.get(i, j)
                    bestCol = j
                }
            }
            
            if (bestCol != -1) {
                assignment[i] = bestCol
                usedCols[bestCol] = true
            }
        }
        
        return assignment
    }
    
    /**
     * Make a matrix doubly stochastic using Sinkhorn algorithm.
     */
    private fun makeDoublyStochastic(M: Matrix): Triple<Matrix, Matrix, Matrix> {
        val result = M.copy()
        val X = Matrix.eye(n)
        val Y = Matrix.eye(n)
        
        var maxRowError = Inf
        var maxColError = Inf
        
        while (maxRowError > alpha2 || maxColError > alpha2) {
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
            
            maxColError = newColSums.maxOf { abs(it - 1.0) }
            maxRowError = newRowSums.maxOf { abs(it - 1.0) }
        }
        
        return Triple(result, X, Y)
    }
    
    /**
     * Compute permanent approximation at each step for analysis.
     */
    private fun computePermStep() {
        val steps = mutableListOf<Double>()
        var cumAccepted = 0
        
        for (i in sampleAccepted.indices) {
            cumAccepted += sampleAccepted[i]
            steps.add(rescalingConstant * cumAccepted / (i + 1))
        }
        
        permStep = steps
    }
}