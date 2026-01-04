package jline.lib.mom

import jline.lib.mom.solver.LinearSolver
import jline.lib.mom.solver.MomSolver
import org.apache.commons.math3.linear.MatrixUtils

fun main() {
    println("Testing MOM Solver...")
    
    // Test 1: Simple 1 station, 2 class network
    println("\n=== Test 1: Single Station, Two Classes ===")
    val L1 = MatrixUtils.createRealMatrix(arrayOf(
        doubleArrayOf(1.0, 2.0)  // Service rates for class 1 and 2
    ))
    val N1 = intArrayOf(5, 3)  // 5 customers of class 1, 3 of class 2
    val Z1 = doubleArrayOf(0.0, 0.0)  // No think time
    
    try {
        val solver1 = MomSolver()
        val result1 = solver1.solve(L1, N1, Z1)
        
        println("MomSolver Results:")
        println("Throughput: ${result1.X.getRow(0).contentToString()}")
        println("Queue lengths: ${result1.Q.getRow(0).contentToString()}")
        println("Normalizing constants size: ${result1.G.size}")
    } catch (e: Exception) {
        println("MomSolver failed: ${e.message}")
        e.printStackTrace()
    }
    
    // Test 2: Two stations, two classes with LinearSolver
    println("\n=== Test 2: Two Stations, Two Classes (LinearSolver) ===")
    val L2 = MatrixUtils.createRealMatrix(arrayOf(
        doubleArrayOf(10.0, 5.0),
        doubleArrayOf(8.0, 12.0)
    ))
    val N2 = intArrayOf(10, 10)
    val Z2 = doubleArrayOf(1.0, 2.0)
    
    try {
        val solver2 = LinearSolver()
        val result2 = solver2.solve(L2, N2, Z2)
        
        println("LinearSolver Results:")
        for (i in 0 until result2.X.rowDimension) {
            println("Station $i throughput: ${result2.X.getRow(i).contentToString()}")
        }
        for (i in 0 until result2.Q.rowDimension) {
            println("Station $i queue length: ${result2.Q.getRow(i).contentToString()}")
        }
    } catch (e: Exception) {
        println("LinearSolver failed: ${e.message}")
        e.printStackTrace()
    }
    
    println("\n=== Tests Complete ===")
}