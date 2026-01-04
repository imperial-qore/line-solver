/**
 * @file Q_RAP_RAP_1 - RAP/RAP/1 Queue Analyzer
 *
 * Computes queue length distribution for a RAP/RAP/1/FCFS queue.
 * RAP = Rational Arrival Process, a generalization of MAP.
 *
 * Based on the Q-MAM library by Benny Van Houdt.
 *
 * @since LINE 3.1.0
 */
package jline.lib.qmam

import jline.lib.smc.QBD_CR
import jline.lib.smc.QBD_FI
import jline.lib.smc.QBD_pi
import jline.lib.smc.stat
import jline.util.matrix.Matrix

/**
 * Result of RAP/RAP/1 queue analysis
 */
data class RAPRAP1Result(
    val queueLength: Matrix   // Queue length distribution
)

/**
 * Options for RAP/RAP/1 queue analysis
 */
data class RAPRAP1Options(
    val mode: String = "CR",
    val maxNumComp: Int = 50,
    val verbose: Int = 0
)

/**
 * Computes queue length distribution for a RAP/RAP/1/FCFS queue.
 *
 * @param C0 RAP arrival process matrix C0 (mA x mA)
 * @param C1 RAP arrival process matrix C1 (mA x mA)
 * @param D0 RAP service process matrix D0 (mS x mS)
 * @param D1 RAP service process matrix D1 (mS x mS)
 * @param options Solver options
 * @return RAPRAP1Result containing queue length distribution
 */
fun qRapRap1(
    C0: Matrix,
    C1: Matrix,
    D0: Matrix,
    D1: Matrix,
    options: RAPRAP1Options = RAPRAP1Options()
): RAPRAP1Result {
    val mA = C0.numRows
    val mS = D0.numRows
    val mtot = mA * mS

    // Validate dimensions
    require(C0.numCols == mA && C1.numRows == mA && C1.numCols == mA) {
        "Arrival process matrices must be mA x mA"
    }
    require(D0.numCols == mS && D1.numRows == mS && D1.numCols == mS) {
        "Service process matrices must be mS x mS"
    }

    // Test the load of the queue
    val eyeA = Matrix.eye(mA)
    val eyeS = Matrix.eye(mS)

    val sumA = C0.add(C1).add(eyeA)
    val sumS = D0.add(D1).add(eyeS)

    val piA = stat(sumA)
    val piS = stat(sumS)

    val lambda = piA.mult(C1).elementSum()
    val mu = piS.mult(D1).elementSum()

    val load = lambda / mu
    require(load < 1) { "The load $load of the system exceeds one" }

    // Compute QBD blocks A0, A1, A2
    val A0 = eyeA.kron(D1)
    val A1 = C0.kron(eyeS).add(eyeA.kron(D0))
    val A2 = C1.kron(eyeS)

    val B0 = A0.copy()
    val B1 = C0.kron(eyeS)

    // Compute G and R using appropriate solver with RAPComp=1
    val qbdResult = if (options.mode.contains("FI")) {
        QBD_FI(A0, A1, A2, null, if (options.verbose > 0) 1 else null, null, null, 1)
    } else {
        QBD_CR(A0, A1, A2, null, if (options.verbose > 0) 1 else null, null, 1)
    }

    val R = qbdResult["R"]!!

    // Compute stationary distribution
    val pi = QBD_pi(B0, B1, R, options.maxNumComp, options.verbose, null, 1)

    // Compute queue length distribution
    val numLevels = pi.numCols / mtot
    val ql = Matrix(1, numLevels)
    for (i in 0 until numLevels) {
        var sum = 0.0
        for (j in 0 until mtot) {
            sum += pi[0, i * mtot + j]
        }
        ql[0, i] = sum
    }

    return RAPRAP1Result(ql)
}
