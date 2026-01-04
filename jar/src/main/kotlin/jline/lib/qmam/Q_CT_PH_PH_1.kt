/**
 * @file Q_CT_PH_PH_1 - Continuous-Time PH/PH/1 Queue Analyzer
 *
 * Computes queue length and waiting time distribution for a continuous-time
 * PH/PH/1/FCFS queue.
 *
 * Based on the Q-MAM library by Benny Van Houdt.
 *
 * @since LINE 3.1.0
 */
package jline.lib.qmam

import jline.io.line_warning
import jline.lib.smc.stat
import jline.util.matrix.Matrix

/**
 * Result of PH/PH/1 queue analysis
 */
data class PHPH1Result(
    val queueLength: Matrix,        // Queue length distribution
    val waitAlpha: Matrix?,         // Waiting time PH alpha vector
    val waitT: Matrix?              // Waiting time PH matrix
)

/**
 * Options for PH/PH/1 queue analysis
 */
data class PHPH1Options(
    val maxNumComp: Int = 1000,
    val verbose: Int = 0
)

/**
 * Computes queue length and waiting time distribution for a PH/PH/1/FCFS queue.
 *
 * @param alpha PH arrival process initial vector (1 x ma)
 * @param T PH arrival process rate matrix (ma x ma)
 * @param beta PH service time initial vector (1 x ms)
 * @param S PH service time rate matrix (ms x ms)
 * @param options Solver options
 * @return PHPH1Result containing queue length and optionally waiting time distribution
 */
fun qCtPhPh1(
    alpha: Matrix,
    T: Matrix,
    beta: Matrix,
    S: Matrix,
    options: PHPH1Options = PHPH1Options()
): PHPH1Result {
    // Validate dimensions
    val ma = alpha.numCols
    val ms = beta.numCols

    require(T.numRows == ma && T.numCols == ma) { "T matrix dimensions must match alpha" }
    require(S.numRows == ms && S.numCols == ms) { "S matrix dimensions must match beta" }

    // Arrival process
    val t = T.mult(Matrix.ones(ma, 1)).scale(-1.0)
    val invNegT = T.scale(-1.0).inv()
    val avgT = alpha.mult(invNegT).mult(Matrix.ones(ma, 1))[0, 0]

    // Service process
    val s = S.mult(Matrix.ones(ms, 1)).scale(-1.0)
    val invNegS = S.scale(-1.0).inv()
    val avgS = beta.mult(invNegS).mult(Matrix.ones(ms, 1))[0, 0]

    val mtot = ms * ma
    val rho = avgS / avgT

    require(rho < 1) { "The load $rho of the system exceeds one" }

    // Compute classic QBD blocks A0, A1, A2
    val eyeMa = Matrix.eye(ma)
    val eyeMs = Matrix.eye(ms)

    val A0 = t.mult(alpha).kron(eyeMs)
    val A1 = T.kron(eyeMs).add(eyeMa.kron(S))

    // Compute QBD blocks in Latouche-Ramaswami approach
    val invmA1 = A1.scale(-1.0).inv()
    val tEye = t.kron(eyeMs)
    val alphaEye = alpha.kron(eyeMs)
    val eyeBeta = eyeMa.kron(beta)
    val eyeS = eyeMa.kron(s)

    val A0pp = alphaEye.mult(invmA1).mult(tEye)
    val A0mp = eyeBeta.mult(invmA1).mult(tEye)
    val A2pm = alphaEye.mult(invmA1).mult(eyeS)
    val A2mm = eyeBeta.mult(invmA1).mult(eyeS)

    // Construct block matrices
    val A0n = Matrix(ms + ma, ms + ma)
    val A2n = Matrix(ms + ma, ms + ma)

    // A0n = [A0pp zeros(ms,ma); A0mp zeros(ma)]
    for (i in 0 until ms) {
        for (j in 0 until ms) {
            A0n[i, j] = A0pp[i, j]
        }
    }
    for (i in 0 until ma) {
        for (j in 0 until ms) {
            A0n[ms + i, j] = A0mp[i, j]
        }
    }

    // A2n = [zeros(ms) A2pm; zeros(ma,ms) A2mm]
    for (i in 0 until ms) {
        for (j in 0 until ma) {
            A2n[i, ms + j] = A2pm[i, j]
        }
    }
    for (i in 0 until ma) {
        for (j in 0 until ma) {
            A2n[ms + i, ms + j] = A2mm[i, j]
        }
    }

    // Compute matrix Gamma: NE corner of matrix G using cyclic reduction
    var itB0 = A0n.copy()
    var itB2 = A2n.copy()
    val Gamma = Matrix(ms, ma)

    // Initialize Gamma from itB2
    for (i in 0 until ms) {
        for (j in 0 until ma) {
            Gamma[i, j] = itB2[i, ms + j]
        }
    }

    var itT = itB0.copy()
    var check = 1.0
    var numit = 1
    val eyeMaMsTotal = Matrix.eye(ma + ms)

    while (check > 1e-13) {
        val itA1 = itB0.mult(itB2).add(itB2.mult(itB0))
        val invFactor = eyeMaMsTotal.sub(itA1).inv()
        itB0 = invFactor.mult(itB0.mult(itB0))
        itB2 = invFactor.mult(itB2.mult(itB2))

        val tmp = itT.mult(itB2)
        for (i in 0 until ms) {
            for (j in 0 until ma) {
                Gamma[i, j] = Gamma[i, j] + tmp[i, ms + j]
            }
        }
        itT = itT.mult(itB0)

        check = Matrix.ones(1, ms).mult(Gamma).mult(Matrix.ones(ma, 1)).sub(Matrix.singleton(1.0)).infinityNorm()
        numit++

    }

    val Gm = Matrix.eye(ma).sub(A0mp.mult(Gamma)).inv().mult(A2mm)
    val RGam = A0pp.mult(Matrix.eye(ms).sub(Gamma.mult(A0mp)).inv())

    // Compute queue length distribution
    val Gstar = invmA1.mult(eyeS.mult(eyeBeta))
        .add(invmA1.mult(tEye).mult(Gamma).mult(Gm).mult(eyeBeta))
    val Rstar = A0.mult(A1.add(A0.mult(Gstar)).scale(-1.0).inv())

    // Compute pi_0
    val betaGamma = beta.mult(Gamma)
    val betaGammaInvNegT = betaGamma.mult(invNegT)
    val normFactor = betaGammaInvNegT.mult(Matrix.ones(ma, 1))[0, 0]
    val pi0Unnorm = betaGammaInvNegT.scale((1 - rho) / normFactor)
    val pi0 = pi0Unnorm.kron(beta)

    // Compute pi_1, pi_2, ...
    val piLevels = mutableListOf<Matrix>()
    piLevels.add(pi0)

    var sumpi = pi0.elementSum()
    numit = 1

    while (sumpi < 1 - 1e-10 && numit < 1 + options.maxNumComp) {
        val piNext = piLevels.last().mult(Rstar)
        piLevels.add(piNext)
        numit++
        sumpi += piNext.elementSum()

        if (options.verbose > 0 && numit % options.verbose == 0) {
            println("Accumulated mass after $numit iterations: $sumpi")
        }
    }

    // Compute queue length distribution
    val ql = Matrix(1, piLevels.size)
    for (i in piLevels.indices) {
        ql[0, i] = piLevels[i].elementSum()
    }

    if (numit == 1 + options.maxNumComp) {
        line_warning("Q_CT_PH_PH_1", "Maximum Number of Components %d reached", numit - 1)
    }

    // Compute waiting time PH representation
    val sigtilde = beta.mult(invNegS).scale(1.0 / beta.mult(invNegS).mult(Matrix.ones(ms, 1))[0, 0])
    val Delta = Matrix.diag(*sigtilde.getRow(0).toArray1D())

    // wait_T = inv(Delta) * (S + RGam * s * beta)' * Delta
    val innerMatrix = S.add(RGam.mult(s).mult(beta))
    val DeltaInv = Delta.inv()
    val waitT = DeltaInv.mult(innerMatrix.transpose()).mult(Delta)

    // Compute wait_alpha
    val sigrho = sigtilde.scale(rho)
    val theta = s.transpose().mult(Delta).scale(-1.0 / beta.mult(invNegS).mult(Matrix.ones(ms, 1))[0, 0])
    val D = DeltaInv.mult(RGam.transpose()).mult(Delta)
    val thetaD = theta.mult(D)
    val thetaDOnes = thetaD.mult(Matrix.ones(ms, 1))[0, 0]

    val prob0 = beta.mult(Matrix.eye(ms).sub(RGam).inv()).mult(Matrix.ones(ms, 1))[0, 0]
    val waitAlpha = thetaD.scale((1 - 1.0 / prob0) / thetaDOnes)

    return PHPH1Result(ql, waitAlpha, waitT)
}
