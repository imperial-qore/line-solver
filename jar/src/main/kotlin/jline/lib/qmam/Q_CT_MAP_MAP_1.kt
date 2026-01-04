/**
 * @file Q_CT_MAP_MAP_1 - Continuous-Time MAP/MAP/1 Queue Analyzer
 *
 * Computes queue length, sojourn time and waiting time distribution for a
 * continuous-time MAP/MAP/1/FCFS queue.
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
 * Result of MAP/MAP/1 queue analysis
 */
data class MAPMAP1Result(
    val queueLength: Matrix,        // Queue length distribution
    val sojAlpha: Matrix?,          // Sojourn time PH alpha vector
    val waitAlpha: Matrix?,         // Waiting time PH alpha vector
    val Smat: Matrix?               // Service time matrix
)

/**
 * Options for MAP/MAP/1 queue analysis
 */
data class MAPMAP1Options(
    val mode: String = "SylvesCR",
    val maxNumComp: Int = 1000,
    val verbose: Int = 0
)

/**
 * Computes queue length and time distributions for a MAP/MAP/1/FCFS queue.
 *
 * @param C0 MAP arrival process matrix D0 (ma x ma)
 * @param C1 MAP arrival process matrix D1 (ma x ma)
 * @param D0 MAP service process matrix D0 (ms x ms)
 * @param D1 MAP service process matrix D1 (ms x ms)
 * @param options Solver options
 * @return MAPMAP1Result containing queue length and time distributions
 */
fun qCtMapMap1(
    C0: Matrix,
    C1: Matrix,
    D0: Matrix,
    D1: Matrix,
    options: MAPMAP1Options = MAPMAP1Options()
): MAPMAP1Result {
    val ma = C0.numRows
    val ms = D0.numRows
    val mtot = ma * ms

    // Validate dimensions
    require(C0.numCols == ma && C1.numRows == ma && C1.numCols == ma) {
        "Arrival process matrices must be ma x ma"
    }
    require(D0.numCols == ms && D1.numRows == ms && D1.numCols == ms) {
        "Service process matrices must be ms x ms"
    }

    // Test the load of the queue
    val invC0 = C0.scale(-1.0).inv()
    val piA = stat(C1.mult(invC0))
    val lambda = piA.mult(C1).elementSum()

    val invD0 = D0.scale(-1.0).inv()
    val piS = stat(D1.mult(invD0))
    val mu = piS.mult(D1).elementSum()

    val load = lambda / mu
    require(load < 1) { "The load $load of the system exceeds one" }

    // Compute classic QBD blocks A0, A1, A2
    val eyeMa = Matrix.eye(ma)
    val eyeMs = Matrix.eye(ms)

    val Am1 = eyeMa.kron(D1)
    val A0 = eyeMa.kron(D0).add(C0.kron(eyeMs))
    val A1 = C1.kron(eyeMs)
    val B0 = C0.kron(eyeMs)

    // Compute G and R using appropriate solver
    val qbdResult = if (options.mode.contains("FI")) {
        QBD_FI(Am1, A0, A1, null, if (options.verbose > 0) 1 else null, null, null, null)
    } else {
        QBD_CR(Am1, A0, A1, null, if (options.verbose > 0) 1 else null, null, null)
    }

    val R = qbdResult["R"]!!

    // Compute stationary distribution
    val stv = QBD_pi(Am1, B0, R, options.maxNumComp, options.verbose)

    // Queue length distribution
    val numLevels = stv.numCols / mtot
    val ql = Matrix(1, numLevels)
    for (i in 0 until numLevels) {
        var sum = 0.0
        for (j in 0 until mtot) {
            sum += stv[0, i * mtot + j]
        }
        ql[0, i] = sum
    }

    // Compute Sojourn and Waiting time PH representation
    val LM = C1.kron(D1)
    val eyeD0 = eyeMa.kron(D0)
    val C0eye = C0.kron(eyeMs)

    // Compute T iteratively
    var Told = Matrix.zeros(mtot, mtot)
    var Tnew = eyeD0.copy()
    var L: Matrix

    if (options.mode.contains("Direct")) {
        // Direct method
        val eyeC0eye = Matrix.eye(mtot).kron(C0eye)
        while (Matrix.infNorm(Told.sub(Tnew)) > 1e-10) {
            Told = Tnew.copy()

            // Solve: vec(L) = -vec(I) * (kron(T',I) + kron(I,C0eye))^{-1}
            val kronTI = Tnew.transpose().kron(Matrix.eye(mtot))
            val system = kronTI.add(eyeC0eye)
            val negEyeVec = Matrix(mtot * mtot, 1)
            for (i in 0 until mtot) {
                negEyeVec[i * mtot + i, 0] = -1.0
            }
            val Lvec = Matrix(mtot * mtot, 1)
            Matrix.solve(system, negEyeVec, Lvec)

            L = Matrix(mtot, mtot)
            for (i in 0 until mtot) {
                for (j in 0 until mtot) {
                    L[i, j] = Lvec[j * mtot + i, 0]
                }
            }

            Tnew = eyeD0.add(L.mult(LM))
        }
        L = Matrix.zeros(mtot, mtot)  // Placeholder for L
    } else {
        // Sylvester method using Schur decomposition
        val schur = schurDecomposition(C0)
        val U = schur.first
        val Tr = schur.second
        val Ukron = U.kron(eyeMs)
        val Trkron = Tr.kron(eyeMs)

        while (Matrix.infNorm(Told.sub(Tnew)) > 1e-10) {
            Told = Tnew.copy()
            L = qSylvest(Ukron, Trkron, Tnew)
            Tnew = eyeD0.add(L.mult(LM))
        }
    }

    // Recompute L for final time distributions
    val schur = schurDecomposition(C0)
    val U = schur.first
    val Tr = schur.second
    L = qSylvest(U.kron(eyeMs), Tr.kron(eyeMs), Tnew)

    // Compute Smat
    val thetaTot = piA.mult(C1).kron(piS).scale(1.0 / (mu * load))

    // Find non-zero entries
    val nonzIndices = mutableListOf<Int>()
    for (i in 0 until mtot) {
        if (thetaTot[0, i] > 0) {
            nonzIndices.add(i)
        }
    }

    if (nonzIndices.isEmpty()) {
        return MAPMAP1Result(ql, null, null, null)
    }

    val nz = nonzIndices.size
    val thetaTotRed = Matrix(1, nz)
    for (i in 0 until nz) {
        thetaTotRed[0, i] = thetaTot[0, nonzIndices[i]]
    }

    val TnewReduced = Matrix(nz, nz)
    for (i in 0 until nz) {
        for (j in 0 until nz) {
            TnewReduced[i, j] = Tnew[nonzIndices[i], nonzIndices[j]]
        }
    }

    val diagThetaRed = Matrix.diag(*thetaTotRed.getRow(0).toArray1D())
    val diagThetaRedInv = diagThetaRed.inv()
    val Smat = diagThetaRedInv.mult(TnewReduced.transpose()).mult(diagThetaRed)

    // Alpha vector of PH representation of Sojourn time
    val D1sumCol = D1.sumRows()
    val sojAlphaFull = Matrix.diag(*thetaTot.getRow(0).toArray1D()).mult(Matrix.ones(ma, 1).kron(D1sumCol)).scale(load / lambda)
    val sojAlpha = Matrix(1, nz)
    for (i in 0 until nz) {
        sojAlpha[0, i] = sojAlphaFull[nonzIndices[i], 0]
    }

    // Alpha vector of PH representation of Waiting time
    val C1sumCol = C1.sumRows()
    val waitAlphaFull = Matrix.diag(*thetaTot.getRow(0).toArray1D()).mult(L).mult(C1sumCol.kron(D1sumCol)).scale(load / lambda)
    val waitAlpha = Matrix(1, nz)
    for (i in 0 until nz) {
        waitAlpha[0, i] = waitAlphaFull[nonzIndices[i], 0]
    }

    return MAPMAP1Result(ql, sojAlpha.transpose(), waitAlpha.transpose(), Smat)
}
