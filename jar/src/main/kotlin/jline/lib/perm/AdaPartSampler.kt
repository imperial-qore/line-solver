package jline.lib.perm

import jline.util.matrix.Matrix
import kotlin.math.abs
import kotlin.math.pow
import kotlin.random.Random

/**
 * Sampling method to approximate the permanent using Adaptive Partitioning (AdaPart).
 * 
 * The AdaPart algorithm uses importance sampling with adaptive partitioning to estimate
 * the permanent. It iteratively refines the sampling space by partitioning based on
 * upper bounds, leading to more efficient sampling than naive approaches.
 * 
 * The algorithm provides three sampling modes:
 * - Classic: Samples until reaching maximum accepted samples
 * - Time: Samples for a maximum time duration
 * - Sample: Samples for a maximum total number of samples
 * 
 * @param matrix The matrix for which to compute the permanent approximation
 * @param maximumAcceptedSamples Maximum number of accepted samples for classic mode (default: 100)
 * @param maximumTime Maximum time in seconds for time mode (default: 30)
 * @param maximumSamples Maximum total samples for sample mode (default: 450)
 * @param mode Sampling mode: "classic", "time", or "sample" (default: "classic")
 * @param solve Whether to automatically run solve() after construction (default: false)
 */
class AdaPartSampler(
    matrix: Matrix,
    private val maximumAcceptedSamples: Int = 100,
    private val maximumTime: Long = 30000, // milliseconds
    private val maximumSamples: Int = 450,
    private val mode: String = "classic",
    solve: Boolean = false
) : PermSolver(matrix) {
    
    private val random = Random.Default
    
    // Tracking variables for analysis
    var sampleAccepted: List<Int> = emptyList()
    var sampleTime: List<Long> = emptyList()
    var permStep: List<Double> = emptyList()
    
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
     * Classic AdaPart sampling - samples until reaching maximum accepted samples.
     */
    private fun samplingClassic(): Double {
        val zUB = soulesBound(matrix)
        
        var accepted = 0
        var total = 0
        
        val startTime = System.currentTimeMillis()
        val acceptedList = mutableListOf<Int>()
        val timeList = mutableListOf<Long>()
        
        while (accepted < maximumAcceptedSamples) {
            val sample = sample()
            
            accepted += sample
            total++
            
            acceptedList.add(sample)
            timeList.add(System.currentTimeMillis() - startTime)
        }
        
        sampleAccepted = acceptedList
        sampleTime = timeList
        computePermStep(zUB)
        
        return zUB * accepted / total
    }
    
    /**
     * Time-constrained AdaPart sampling.
     */
    private fun samplingTime(): Double {
        val zUB = soulesBound(matrix)
        
        var accepted = 0
        var total = 0
        
        val startTime = System.currentTimeMillis()
        val acceptedList = mutableListOf<Int>()
        val timeList = mutableListOf<Long>()
        
        while (System.currentTimeMillis() - startTime < maximumTime) {
            val sample = sample()
            
            accepted += sample
            total++
            
            acceptedList.add(sample)
            timeList.add(System.currentTimeMillis() - startTime)
        }
        
        sampleAccepted = acceptedList
        sampleTime = timeList
        computePermStep(zUB)
        
        return if (total > 0) zUB * accepted / total else 0.0
    }
    
    /**
     * Sample-count-constrained AdaPart sampling.
     */
    private fun samplingSample(): Double {
        val zUB = soulesBound(matrix)
        
        var accepted = 0
        var total = 0
        
        val startTime = System.currentTimeMillis()
        val acceptedList = mutableListOf<Int>()
        val timeList = mutableListOf<Long>()
        
        while (acceptedList.size < maximumSamples) {
            val sample = sample()
            
            accepted += sample
            total++
            
            acceptedList.add(sample)
            timeList.add(System.currentTimeMillis() - startTime)
        }
        
        sampleAccepted = acceptedList
        sampleTime = timeList
        computePermStep(zUB)
        
        return zUB * accepted / total
    }
    
    /**
     * Generate one sample using adaptive partitioning.
     * 
     * @return 1 if sample is accepted, 0 if rejected
     */
    private fun sample(): Int {
        // Store current restriction space
        var S = mutableSetOf(IntArray(n) { n }.toList())
        
        val startTime = System.currentTimeMillis()
        
        // Until all elements of the permutation are defined
        while (S.any { it.contains(n) } && 
               (mode != "time" || System.currentTimeMillis() - startTime < maximumTime)) {
            
            // Get the current partial permutation
            val sInit = S.first()
            
            // Get the state matrix
            val sMatrix = modifyMatrix(matrix, sInit)
            
            // Compute the bound
            var ub = soulesBound(sMatrix)
            val zubS = ub
            
            // Partitioning loop
            var init = true
            while ((ub >= zubS || init) && 
                   (mode != "time" || System.currentTimeMillis() - startTime < maximumTime)) {
                
                init = false
                
                // Get a random partial permutation
                val sSub = S.random(random)
                S.remove(sSub)
                
                // Get the current state matrix
                val subMatrix = modifyMatrix(matrix, sSub)
                
                // Select column and modify bound
                val (newUb, j) = selectColumn(subMatrix, zubS, ub)
                ub = newUb
                
                // Modify the restriction space
                for (i in 0 until n) {
                    val sAdd = sSub.toMutableList()
                    if (i !in sSub) {
                        sAdd[j] = i
                        S.add(sAdd)
                    }
                }
            }
            
            // Compute probabilities and select random
            val c = computeProbabilities(S, zubS)
            
            // If slack selected, reject
            if (c == S.size) {
                return 0
            }
            
            // Otherwise, initialize restriction space with selected permutation
            S = subset(S, c)
        }
        
        return 1
    }
    
    /**
     * Select the column to fix during partitioning.
     */
    private fun selectColumn(sMatrix: Matrix, zubS: Double, ub: Double): Pair<Double, Int> {
        val ubi = DoubleArray(n)
        
        for (i in 0 until n) {
            for (j in 0 until n) {
                val L = IntArray(n) { n }
                L[i] = j
                val fB = modifyMatrix(sMatrix, L.toList())
                ubi[i] += soulesBound(fB)
            }
        }
        
        val j = ubi.indices.minByOrNull { ubi[it] } ?: 0
        val newUb = ub - zubS + ubi[j]
        
        return Pair(newUb, j)
    }
    
    /**
     * Compute probabilities and randomly select a permutation.
     */
    private fun computeProbabilities(S: Set<List<Int>>, zubS: Double): Int {
        val p = S.map { soulesBound(modifyMatrix(matrix, it)) }.toMutableList()
        
        // Normalize
        val sum = p.sum()
        if (sum > 0) {
            for (i in p.indices) {
                p[i] /= zubS
            }
        }
        
        // Add slack
        val slack = 1.0 - p.sum()
        p.add(slack)
        
        // Correct for rounding errors
        val totalProb = p.sum()
        if (totalProb > 0) {
            for (i in p.indices) {
                p[i] = abs(p[i]) / totalProb
            }
        }
        
        // Randomly select one
        val rand = random.nextDouble()
        var cumSum = 0.0
        for (i in p.indices) {
            cumSum += p[i]
            if (rand <= cumSum) {
                return i
            }
        }
        
        return p.size - 1
    }
    
    /**
     * Generate new subset from current set and selected index.
     */
    private fun subset(S: Set<List<Int>>, c: Int): MutableSet<List<Int>> {
        val sInter = S.elementAt(c).toMutableList()
        
        // If only one missing item in the permutation, fill it
        val missingCount = sInter.count { it == n }
        if (missingCount == 1) {
            val existingValues = sInter.filter { it != n }.toSet()
            val missingValue = (0 until n).first { it !in existingValues }
            val missingIndex = sInter.indexOf(n)
            sInter[missingIndex] = missingValue
        }
        
        return mutableSetOf(sInter)
    }
    
    /**
     * Modify a matrix using a partial permutation.
     */
    private fun modifyMatrix(M: Matrix, t: List<Int>): Matrix {
        val result = Matrix.zeros(n, n)
        
        // Create mask for fixed parts
        val mask = Array(n) { BooleanArray(n) }
        
        for (j in t.indices) {
            if (t[j] != n) {
                mask[t[j]][j] = true
            }
        }
        
        // Add non-fixed parts
        for (i in 0 until n) {
            if (i !in t) {
                for (j in t.indices) {
                    if (t[j] == n) {
                        mask[i][j] = true
                    }
                }
            }
        }
        
        // Apply mask to matrix
        for (i in 0 until n) {
            for (j in 0 until n) {
                if (mask[i][j]) {
                    result.set(i, j, M.get(i, j))
                }
            }
        }
        
        return result
    }
    
    /**
     * Compute Soules bound for permanent estimation.
     */
    private fun soulesBound(M: Matrix): Double {
        // Compute gamma values
        val gamma = DoubleArray(n + 1)
        gamma[0] = 0.0
        
        for (k in 1..n) {
            var factorial = 1.0
            for (i in 1..k) {
                factorial *= i
            }
            gamma[k] = factorial.pow(1.0 / k)
        }
        
        // Compute delta
        val delta = DoubleArray(n)
        for (i in 0 until n) {
            delta[i] = gamma[n - i] - gamma[n - i - 1]
        }
        
        // Sort matrix columns and compute bound
        val mSum = DoubleArray(n)
        for (j in 0 until n) {
            val column = DoubleArray(n) { i -> M.get(i, j) }
            column.sort()
            
            for (i in 0 until n) {
                mSum[j] += column[i] * delta[i]
            }
        }
        
        var product = 1.0
        for (sum in mSum) {
            product *= sum
        }
        
        return product
    }
    
    /**
     * Compute permanent approximation at each step for analysis.
     */
    private fun computePermStep(zUB: Double) {
        val steps = mutableListOf<Double>()
        var cumAccepted = 0
        
        for (i in sampleAccepted.indices) {
            cumAccepted += sampleAccepted[i]
            steps.add(zUB * cumAccepted / (i + 1))
        }
        
        permStep = steps
    }
}