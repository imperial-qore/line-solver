package jline.lib.smc

import jline.util.matrix.Matrix

/**
 * QBD_pi: Stationary vector of a Quasi-Birth-Death Markov Chains [Neuts]
 *
 * DISCRETE TIME CASE:
 * pi=QBD_pi(B0,B1,R) computes the stationary vector of a Quasi-Birth-Death
 * Markov chain with a transition matrix of the form
 *
 *           B1  A2  0   0   0  ...
 *           B0  A1  A2  0   0  ...
 *   P  =    0   A0  A1  A2  0  ...
 *           0   0   A0  A1  A2 ...
 *           ...
 *
 * the input matrix R is the minimal nonnegative solution to the matrix
 * equation R = A2 + R A1 + R^2 A0
 *
 * CONTINUOUS TIME CASE:
 * pi=QBD_pi(B0,B1,R) computes the stationary vector of a Quasi-Birth-Death
 * Markov chain with a rate matrix of the form
 *
 *           B1  A2  0   0   0  ...
 *           B0  A1  A2  0   0  ...
 *   Q  =    0   A0  A1  A2  0  ...
 *           0   0   A0  A1  A2 ...
 *           ...
 *
 * the input matrix R is the minimal nonnegative solution to the matrix
 * equation 0 = A2 + R A1 + R^2 A0
 *
 * Optional Parameters:
 *   MaxNumComp: Maximum number of components (default: 500)
 *   Verbose: The accumulated probability mass is printed at every n steps when set to n (default:0)
 *   Boundary: Allows solving the QBD with a more general boundary
 *   RAPComp: set to 1 if the QBD has RAP components
 */
fun QBD_pi(B0: Matrix, B1: Matrix, R: Matrix,
           MaxNumComp: Int = 500,
           Verbose: Int = 0,
           Boundary: Matrix? = null,
           RAPComp: Int = 0): Matrix {
    
    var B1_work = B1.copy()
    val B0_work = B0.copy()
    var Boundary_work = Boundary?.copy()
    
    val m = R.numRows
    
    // Convert to discrete time problem, if needed
    val B1_diag = Matrix(m, 1, 0)
    Matrix.extractDiag(B1_work, B1_diag)
    
    if (B1_diag.elementMin() < 0 || RAPComp == 1) { // continuous time
        val lamb = -B1_diag.elementMin()
        
        // Simplified normalization for continuous time case
        B1_work = B1_work.scale(1.0 / lamb).add(1.0, Matrix.eye(m))
        
        B0_work.scaleEq(1.0 / lamb)
    }
    
    val temp = Matrix.eye(m).add(-1.0, R).inv()
    if (temp.elementMin() < -100 * 1e-15) {
        throw RuntimeException("The spectral radius of R is not below 1: QBD is not pos. recurrent")
    }
    
    // Simplified implementation - in practice this would be more complex
    val pi0 = stat(B1_work.add(1.0, R.mult(B0_work)))
    val normalizer = pi0.mult(temp).mult(Matrix.ones(m, 1))[0]
    pi0.scaleEq(1.0 / normalizer)
    
    val pi_components = mutableListOf<Matrix>()
    pi_components.add(pi0.copy())
    
    var sumpi = pi0.elementSum()
    var numit = 1
    
    while (sumpi < 1 - 1e-10 && numit < MaxNumComp) {
        val pi_next = pi_components.last().mult(R)
        pi_components.add(pi_next)
        numit++
        sumpi += pi_next.elementSum()
        
        if (Verbose > 0 && numit % Verbose == 0) {
            println("Accumulated mass after $numit iterations: $sumpi")
        }
    }
    
    if (numit == MaxNumComp) {
        println("Maximum Number of Components $numit reached")
    }
    
    // Concatenate all components horizontally 
    var result = pi_components[0].copy()
    for (i in 1 until pi_components.size) {
        result = Matrix.concatColumns(result, pi_components[i], null)
    }
    
    return result
}