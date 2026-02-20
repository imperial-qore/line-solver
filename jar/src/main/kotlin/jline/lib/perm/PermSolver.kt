package jline.lib.perm

import jline.util.matrix.Matrix

/**
 * Abstract base class for permanent computation solvers.
 * 
 * This abstract class provides the common interface for all permanent computation algorithms.
 * All solvers inherit from this class and implement the compute() method for their specific algorithm.
 * The base class provides the solve() method that measures time and memory usage during computation.
 * 
 * @property matrix The matrix for which to compute the permanent
 * @property n The dimension of the matrix
 * @property value The computed permanent value (set by compute() method)
 * @property time The computation time in milliseconds
 * @property memory The memory usage in bytes
 */
abstract class PermSolver(val matrix: Matrix) {
    
    val n: Int = matrix.getNumRows()
    var value: Double = 0.0
    var time: Long = 0L
    var memory: Long = 0L
    
    /**
     * Compute the permanent or approximation for the given matrix.
     * This method must be implemented by all concrete solver classes.
     * The result should be stored in the value property.
     */
    abstract fun compute()
    
    /**
     * Measure time and memory usage of the compute method.
     * This method wraps the compute() call with performance monitoring.
     */
    fun solve() {
        // Get initial memory usage
        val runtime = Runtime.getRuntime()
        runtime.gc() // Suggest garbage collection before measurement
        val memoryBefore = runtime.totalMemory() - runtime.freeMemory()
        
        // Measure computation time
        val startTime = System.currentTimeMillis()
        
        compute()
        
        val endTime = System.currentTimeMillis()
        this.time = endTime - startTime
        
        // Measure memory usage after computation
        val memoryAfter = runtime.totalMemory() - runtime.freeMemory()
        this.memory = memoryAfter - memoryBefore
    }
}