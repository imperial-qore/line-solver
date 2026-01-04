package jline.lib.mom

import jline.util.Maths
import jline.util.matrix.Matrix

fun main() {
    println("=== MOM Solver Integration Test ===\n")
    
    // Test 1: Using jline.util.Maths functions
    println("Test 1: Binomial coefficient (10 choose 3)")
    val binCoeff = Maths.binomialCoeff(10, 3)
    println("Result: $binCoeff")
    
    // Test 2: Using multichoose
    println("\nTest 2: Multichoose(3, 2)")
    val multiChoose = Maths.multichoose(3.0, 2.0)
    println("Results (${multiChoose.getNumRows()} combinations):")
    for (i in 0 until multiChoose.getNumRows()) {
        print("  [")
        for (j in 0 until multiChoose.getNumCols()) {
            print("${multiChoose.get(i, j).toInt()}")
            if (j < multiChoose.getNumCols() - 1) print(", ")
        }
        println("]")
    }
    
    // Test 3: Using Matrix operations
    println("\nTest 3: Matrix operations")
    val matrix1 = Matrix(2, 2)
    matrix1.set(0, 0, 1.0)
    matrix1.set(0, 1, 2.0)
    matrix1.set(1, 0, 3.0)
    matrix1.set(1, 1, 4.0)
    
    val matrix2 = Matrix(2, 2)
    matrix2.set(0, 0, 5.0)
    matrix2.set(0, 1, 6.0)
    matrix2.set(1, 0, 7.0)
    matrix2.set(1, 1, 8.0)
    
    val result = matrix1.mult(matrix2)
    println("Matrix multiplication result:")
    for (i in 0 until result.getNumRows()) {
        print("  [")
        for (j in 0 until result.getNumCols()) {
            print("${result.get(i, j)}")
            if (j < result.getNumCols() - 1) print(", ")
        }
        println("]")
    }
    
    // Test 4: MOM-specific utilities
    println("\nTest 4: MOM utilities")
    val combinations = arrayOf(
        intArrayOf(2, 0, 0),
        intArrayOf(1, 1, 0),
        intArrayOf(0, 2, 0),
        intArrayOf(1, 0, 1),
        intArrayOf(0, 1, 1),
        intArrayOf(0, 0, 2)
    )
    
    println("Original combinations:")
    combinations.forEach { println("  ${it.contentToString()}") }
    
    val sorted = jline.lib.mom.util.MomUtils.sortByNnzPos(combinations)
    println("\nSorted by non-zeros:")
    sorted.forEach { println("  ${it.contentToString()}") }
    
    println("\n=== All tests completed successfully ===")
}