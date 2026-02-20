package jline.lib.perm

import jline.util.matrix.Matrix

/**
 * Example usage of the permanent computation algorithms.
 * 
 * This example demonstrates:
 * 1. Creating a test matrix
 * 2. Computing permanent using different algorithms
 * 3. Comparing results and performance
 */
fun main() {
    // Create a small test matrix
    val data = arrayOf(
        doubleArrayOf(1.0, 2.0, 3.0),
        doubleArrayOf(4.0, 5.0, 6.0),
        doubleArrayOf(7.0, 8.0, 9.0)
    )
    val matrix = Matrix(data)
    
    println("Test Matrix:")
    println("1.0  2.0  3.0")
    println("4.0  5.0  6.0")
    println("7.0  8.0  9.0")
    println()
    
    // Test each solver
    val solvers = listOf(
        "Naive" to NaivePermanent(matrix),
        "Ryzer (Gray Code)" to RyzerPermanent(matrix, mode = "graycode"),
        "Ryzer (Naive)" to RyzerPermanent(matrix, mode = "naive"),
        "Bethe" to BethePermanent(matrix),
        "AdaPart" to AdaPartSampler(matrix, maximumAcceptedSamples = 1000),
        "Huber-Law" to HuberLawSampler(matrix)
    )
    
    println("Permanent Computation Results:")
    println("=" * 50)
    
    for ((name, solver) in solvers) {
        solver.solve()
        println("$name Solver:")
        println("  Permanent value: ${solver.value}")
        println("  Computation time: ${solver.time} ms")
        println("  Memory usage: ${solver.memory} bytes")
        println()
    }
    
    // Example with queueing network
    println("\nQueueing Network Example:")
    println("=" * 50)
    
    val network = NetworkNoThink(
        numberOfQueues = 2,
        numberOfClasses = 1,
        numberPerClass = intArrayOf(3),
        meanServiceDemand = arrayOf(
            doubleArrayOf(1.0),
            doubleArrayOf(2.0)
        )
    )
    
    val state = intArrayOf(1, 2)  // 1 job at queue 1, 2 jobs at queue 2
    val betheSolver = BethePermanent(matrix) // Dummy matrix, will be replaced
    
    val (prob, time, memory) = network.marginal(betheSolver, state, preprocessing = true)
    println("Marginal probability for state [1, 2]: $prob")
    println("Time: $time ms, Memory: $memory bytes")
}

/**
 * String multiplication extension for creating separators.
 */
private operator fun String.times(n: Int): String = this.repeat(n)