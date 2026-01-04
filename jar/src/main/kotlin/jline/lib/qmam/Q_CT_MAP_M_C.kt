/**
 * @file Q_CT_MAP_M_C - Continuous-Time MAP/M/c Queue Analyzer
 *
 * Computes queue length and waiting time distribution for a
 * continuous-time MAP/M/c/FCFS queue.
 *
 * Based on the Q-MAM library by Benny Van Houdt.
 *
 * @since LINE 3.1.0
 */
package jline.lib.qmam

import jline.lib.smc.QBD_CR
import jline.lib.smc.QBD_FI
import jline.lib.smc.stat
import jline.util.matrix.Matrix

/**
 * Result of MAP/M/c queue analysis
 */
data class MAPMcResult(
    val queueLength: Matrix,        // Queue length distribution
    val waitAlpha: Matrix?,         // Waiting time PH alpha vector
    val Smat: Matrix?               // Service time matrix
)

/**
 * Options for MAP/M/c queue analysis
 */
data class MAPMcOptions(
    val mode: String = "SylvesCR",
    val maxNumComp: Int = 1000,
    val verbose: Int = 0
)

/**
 * Computes queue length and waiting time distribution for a MAP/M/c/FCFS queue.
 *
 * @param D0 MAP arrival process matrix D0 (m x m)
 * @param D1 MAP arrival process matrix D1 (m x m)
 * @param mu Exponential service rate
 * @param c Number of servers
 * @param options Solver options
 * @return MAPMcResult containing queue length and waiting time distribution
 */
fun qCtMapMC(
    D0: Matrix,
    D1: Matrix,
    mu: Double,
    c: Int,
    options: MAPMcOptions = MAPMcOptions()
): MAPMcResult {
    val m = D0.numRows

    // Validate inputs
    require(D0.numCols == m && D1.numRows == m && D1.numCols == m) {
        "D0 and D1 must be m x m matrices"
    }
    require(mu > 1e-14) { "Service rate mu must be strictly positive" }
    require(c >= 1) { "Number of servers c must be at least 1" }

    // Test the load of the queue
    val invD0 = D0.scale(-1.0).inv()
    val piA = stat(D1.mult(invD0))
    val lambda = piA.mult(D1).elementSum()

    val load = lambda / (mu * c)
    require(load < 1) { "The load $load of the system exceeds one" }

    // Compute QBD blocks
    val eyeM = Matrix.eye(m)
    val A0 = eyeM.scale(c * mu)
    val A1 = D0.sub(eyeM.scale(c * mu))
    val A2 = D1.copy()

    // Compute G and R using appropriate solver
    val qbdResult = if (options.mode.contains("FI")) {
        QBD_FI(A0, A1, A2, null, if (options.verbose > 0) 1 else null, null, null, null)
    } else {
        QBD_CR(A0, A1, A2, null, if (options.verbose > 0) 1 else null, null, null)
    }

    val R = qbdResult["R"]!!

    // Gaver-Jacobs-Latouche Level-Dependent QBD approach
    val piGJL = Matrix(1, c * m)

    if (c > 1) {
        // Compute invC matrices
        val invC = mutableListOf<Matrix>()
        invC.add(D0.scale(-1.0).inv())

        for (i in 2 until c) {
            val innerMatrix = D0.sub(eyeM.scale((i - 1) * mu))
                .add(invC[i - 2].mult(D1).scale((i - 1) * mu))
                .scale(-1.0)
            invC.add(innerMatrix.inv())
        }

        // Compute pi at level c-1
        val boundaryMatrix = D0.sub(eyeM.scale((c - 1) * mu))
            .add(R.mult(A0))
            .add(invC[c - 2].mult(D1).scale((c - 1) * mu))
            .add(eyeM)
        val piCm1 = stat(boundaryMatrix)

        // Copy to piGJL at position (c-1)*m
        for (j in 0 until m) {
            piGJL[0, (c - 1) * m + j] = piCm1[0, j]
        }

        // Compute lower levels
        for (i in c - 1 downTo 1) {
            val piNext = Matrix(1, m)
            for (j in 0 until m) {
                piNext[0, j] = piGJL[0, i * m + j]
            }
            val piPrev = piNext.mult(invC[i - 1]).scale(i * mu)
            for (j in 0 until m) {
                piGJL[0, (i - 1) * m + j] = piPrev[0, j]
            }
        }
    } else {
        // c = 1 case
        val boundaryMatrix = D0.add(R.mult(A0)).add(eyeM)
        val pi0 = stat(boundaryMatrix)
        for (j in 0 until m) {
            piGJL[0, j] = pi0[0, j]
        }
    }

    // Normalize
    val ImR = eyeM.sub(R)
    val ImRinv = ImR.inv()

    var K = 0.0
    for (i in 0 until (c - 1) * m) {
        K += piGJL[0, i]
    }

    val piCm1 = Matrix(1, m)
    for (j in 0 until m) {
        piCm1[0, j] = piGJL[0, (c - 1) * m + j]
    }
    K += piCm1.mult(ImRinv).mult(Matrix.ones(m, 1))[0, 0]

    piGJL.scaleEq(1.0 / K)

    // Compute higher levels
    var piC1 = Matrix(1, m)
    for (j in 0 until m) {
        piC1[0, j] = piGJL[0, (c - 1) * m + j]
    }

    val piLevels = mutableListOf<Matrix>()
    piLevels.add(piC1.copy())

    var sumpi = piGJL.elementSum()
    var numit = 1

    while (sumpi < 1 - 1e-10 && numit < 1 + options.maxNumComp - c) {
        val piNext = piLevels.last().mult(R)
        piLevels.add(piNext)
        numit++
        sumpi += piNext.elementSum()

        if (options.verbose > 0 && numit % options.verbose == 0) {
            println("Accumulated mass after $numit iterations: $sumpi")
        }
    }

    // Queue length distribution
    val ql = Matrix(1, c - 1 + piLevels.size)

    // First c-1 levels from piGJL
    for (i in 0 until c - 1) {
        var levelSum = 0.0
        for (j in 0 until m) {
            levelSum += piGJL[0, i * m + j]
        }
        ql[0, i] = levelSum
    }

    // Remaining levels from piLevels
    for (i in piLevels.indices) {
        ql[0, c - 1 + i] = piLevels[i].elementSum()
    }

    // Waiting time distribution
    // Probability of zero waiting
    var probZero = 0.0
    for (i in 0 until c * m) {
        probZero += piGJL[0, i]
    }
    val D1sumCol = D1.sumRows()

    var numerator = 0.0
    for (i in 0 until c) {
        val piLevel = Matrix(1, m)
        for (j in 0 until m) {
            piLevel[0, j] = piGJL[0, i * m + j]
        }
        numerator += piLevel.mult(D1sumCol)[0, 0]
    }

    // Concatenate all piT
    val piTsize = (c - 1) * m + piLevels.size * m
    val piT = Matrix(1, piTsize)
    for (i in 0 until (c - 1) * m) {
        piT[0, i] = piGJL[0, i]
    }
    for (i in piLevels.indices) {
        for (j in 0 until m) {
            piT[0, (c - 1) * m + i * m + j] = piLevels[i][0, j]
        }
    }

    var denominator = 0.0
    val numPiTLevels = piTsize / m
    for (i in 0 until numPiTLevels) {
        val piLevel = Matrix(1, m)
        for (j in 0 until m) {
            piLevel[0, j] = piT[0, i * m + j]
        }
        denominator += piLevel.mult(D1sumCol)[0, 0]
    }

    probZero = numerator / denominator

    // Compute alpha vector
    val piCm1ForAlpha = Matrix(1, m)
    for (j in 0 until m) {
        piCm1ForAlpha[0, j] = piGJL[0, (c - 1) * m + j]
    }
    val temp = piCm1ForAlpha.mult(ImRinv).mult(D1)
    val alphaVec = temp.scale(1.0 / temp.elementSum())

    // Compute T matrix iteratively
    var Told = Matrix.zeros(m, m)
    var Tnew = A0.scale(-1.0).copy()
    var L: Matrix
    var TIterCount = 0

    val maxTIterations = 1000
    if (options.mode.contains("Direct")) {
        val eyeD0eye = eyeM.kron(D0)
        var tNorm = Matrix.infNorm(Told.sub(Tnew))
        while (tNorm > 1e-10 && TIterCount < maxTIterations) {
            Told = Tnew.copy()

            val kronTI = Tnew.transpose().kron(eyeM)
            val system = kronTI.add(eyeD0eye)
            val negEyeVec = Matrix(m * m, 1)
            for (i in 0 until m) {
                negEyeVec[i * m + i, 0] = -1.0
            }
            val Lvec = Matrix(m * m, 1)
            Matrix.solve(system, negEyeVec, Lvec)

            L = Matrix(m, m)
            for (i in 0 until m) {
                for (j in 0 until m) {
                    L[i, j] = Lvec[j * m + i, 0]
                }
            }

            Tnew = A0.scale(-1.0).add(L.mult(D1).scale(mu * c))
            TIterCount++
            tNorm = Matrix.infNorm(Told.sub(Tnew))
        }
    } else {
        val schur = schurDecomposition(D0)
        val U = schur.first
        val Tr = schur.second
        var tNorm = Matrix.infNorm(Told.sub(Tnew))
        while (tNorm > 1e-10 && TIterCount < maxTIterations) {
            Told = Tnew.copy()
            L = qSylvest(U, Tr, Tnew)
            Tnew = A0.scale(-1.0).add(L.mult(D1).scale(mu * c))
            TIterCount++
            tNorm = Matrix.infNorm(Told.sub(Tnew))
        }
        if (TIterCount >= maxTIterations && tNorm > 1e-10) {
            // Algorithm did not converge, use direct method as fallback
            // For scalar case (m=1), compute directly
            if (m == 1) {
                // Solve: T + L*D1*(mu*c) = -A0 where L satisfies T*L + D0*L = -I
                // For m=1: L = -1/(T + D0), so T = -A0 + (-1/(T+D0))*D1*(mu*c)
                // Fixed point: T = -A0 - D1*(mu*c)/(T + D0)
                // Let x = T[0,0], a = -A0[0,0], b = D1[0,0]*(mu*c), d = D0[0,0]
                // x = a - b/(x + d)  =>  x(x+d) = a(x+d) - b  =>  x^2 + dx = ax + ad - b
                // x^2 + (d-a)x - (ad-b) = 0
                val a = -A0[0, 0]  // = c*mu
                val b = D1[0, 0] * (mu * c)
                val d = D0[0, 0]
                val p = d - a
                val q = -(a * d - b)
                // x = (-p +/- sqrt(p^2 + 4q)) / 2
                val discriminant = p * p + 4 * q
                if (discriminant >= 0) {
                    val sqrtD = kotlin.math.sqrt(discriminant)
                    // Take the negative root (T should be negative for a valid generator)
                    val x1 = (-p + sqrtD) / 2
                    val x2 = (-p - sqrtD) / 2
                    Tnew[0, 0] = if (x2 < 0) x2 else x1
                }
            }
        }
    }

    // Compute Smat for non-zero entries
    val nonzIndices = mutableListOf<Int>()
    for (i in 0 until m) {
        if (alphaVec[0, i] > 0) {
            nonzIndices.add(i)
        }
    }

    if (nonzIndices.isEmpty()) {
        return MAPMcResult(ql, null, null)
    }

    val nz = nonzIndices.size
    val thetaTotRed = Matrix(1, nz)
    for (i in 0 until nz) {
        thetaTotRed[0, i] = alphaVec[0, nonzIndices[i]]
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

    // Waiting time alpha
    val rhoVec = Tnew.add(A0).sumRows()
    var alphaRhoSum = 0.0
    for (i in 0 until m) {
        alphaRhoSum += alphaVec[0, i] * rhoVec[i, 0]
    }

    val waitAlphaFull = Matrix(1, m)
    for (i in 0 until m) {
        waitAlphaFull[0, i] = (1 - probZero) * alphaVec[0, i] * rhoVec[i, 0] / alphaRhoSum
    }

    val waitAlpha = Matrix(1, nz)
    for (i in 0 until nz) {
        waitAlpha[0, i] = waitAlphaFull[0, nonzIndices[i]]
    }

    return MAPMcResult(ql, waitAlpha, Smat)
}
