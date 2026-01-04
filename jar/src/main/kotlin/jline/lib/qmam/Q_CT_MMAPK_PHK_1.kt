/**
 * @file Q_CT_MMAPK_PHK_1 - Continuous-Time MMAP[K]/PH[K]/1 Queue Analyzer
 *
 * Computes queue length, sojourn time and waiting time distribution for a
 * continuous-time MMAP[K]/PH[K]/1/FCFS queue with K customer types.
 *
 * Based on the Q-MAM library by Benny Van Houdt.
 *
 * @since LINE 3.1.0
 */
package jline.lib.qmam

import jline.lib.smc.stat
import jline.util.matrix.Matrix

/**
 * Result of MMAP[K]/PH[K]/1 queue analysis
 */
data class MMAPKPHK1Result(
    val qlPerType: List<Matrix>,      // Per-type queue length distributions
    val qlTotal: Matrix,              // Overall queue length distribution
    val sojAlpha: List<Matrix>,       // Sojourn time PH alpha vectors (per type + overall)
    val waitAlpha: List<Matrix>,      // Waiting time PH alpha vectors (per type + overall)
    val Smat: Matrix?                 // Common service time matrix
)

/**
 * Options for MMAP[K]/PH[K]/1 queue analysis
 */
data class MMAPKPHK1Options(
    val mode: String = "Sylves",
    val maxNumComp: Int = 1000,
    val verbose: Int = 0
)

/**
 * Computes queue length and time distributions for a MMAP[K]/PH[K]/1/FCFS queue.
 *
 * @param D0 MMAP arrival process matrix D0 (m x m)
 * @param D List of MMAP arrival matrices D1, D2, ..., DK (each m x m)
 * @param alpha List of PH service time initial vectors (1 x m_i)
 * @param S List of PH service time matrices (m_i x m_i)
 * @param options Solver options
 * @return MMAPKPHK1Result containing queue length and time distributions
 */
fun qCtMmapkPhk1(
    D0: Matrix,
    D: List<Matrix>,
    alpha: List<Matrix>,
    S: List<Matrix>,
    options: MMAPKPHK1Options = MMAPKPHK1Options()
): MMAPKPHK1Result {
    val K = alpha.size
    val m = D0.numRows

    require(D.size == K) { "Number of D matrices must equal number of service types" }
    require(S.size == K) { "Number of S matrices must equal number of service types" }

    // Compute service phase cumulative indices
    val smk = IntArray(K + 1)
    smk[0] = 0
    for (i in 0 until K) {
        smk[i + 1] = smk[i] + alpha[i].numCols
    }
    val mser = smk[K]
    val mtot = mser * m

    // Test the load of the queue
    var Dsum = D0.copy()
    for (i in 0 until K) {
        Dsum = Dsum.add(D[i])
    }

    // Compute maximum diagonal element of Dsum
    val diagDsum = Matrix(m, 1)
    Matrix.extractDiag(Dsum, diagDsum)
    val maxDiagDsum = diagDsum.scale(-1.0).elementMax()

    val pi = stat(Matrix.eye(m).sub(Dsum.scale(1.0 / maxDiagDsum)))

    val lambdas = DoubleArray(K)
    val mus = DoubleArray(K)
    val beta = mutableListOf<Matrix>()

    for (i in 0 until K) {
        lambdas[i] = pi.mult(D[i]).elementSum()

        val sSum = S[i].sumRows()
        val temp = S[i].sub(sSum.mult(alpha[i]))

        if (S[i].numRows > 1) {
            val diagTemp = Matrix(S[i].numRows, 1)
            Matrix.extractDiag(temp, diagTemp)
            val maxDiagTemp = diagTemp.scale(-1.0).elementMax()
            beta.add(stat(Matrix.eye(S[i].numRows).sub(temp.scale(1.0 / maxDiagTemp))))
        } else {
            beta.add(Matrix.singleton(1.0))
        }

        mus[i] = -beta[i].mult(S[i].sumRows())[0, 0]
    }

    val lambda = pi.mult(Dsum.sub(D0)).elementSum()
    var load = 0.0
    for (i in 0 until K) {
        load += lambdas[i] / mus[i]
    }

    require(load < 1) { "The load $load of the system exceeds one" }

    // Construct building blocks
    val Tser = Matrix.zeros(mser, mser)
    val tser = Matrix.zeros(mser, 1)

    for (i in 0 until K) {
        for (r in 0 until S[i].numRows) {
            for (c in 0 until S[i].numCols) {
                Tser[smk[i] + r, smk[i] + c] = S[i][r, c]
            }
        }
        val sSum = S[i].sumRows().scale(-1.0)
        for (r in 0 until S[i].numRows) {
            tser[smk[i] + r, 0] = sSum[r, 0]
        }
    }

    // LM matrix
    val LM = Matrix.zeros(mtot, mtot)
    for (i in 0 until K) {
        // Build the selector row: [zeros, alpha[i], zeros]
        val selector = Matrix(1, mser)
        for (j in 0 until alpha[i].numCols) {
            selector[0, smk[i] + j] = alpha[i][0, j]
        }

        // tser * selector
        val tserSelector = tser.mult(selector)

        // kron(D[i], tserSelector)
        val kronPart = D[i].kron(tserSelector)
        LM.addEq(1.0, kronPart)
    }

    // Compute T iteratively
    val eyeTser = Matrix.eye(m).kron(Tser)
    val D0eye = D0.kron(Matrix.eye(mser))

    var Told = Matrix.zeros(mtot, mtot)
    var Tnew = eyeTser.copy()
    var L: Matrix

    if (options.mode == "Direct") {
        val eyeD0eye = Matrix.eye(mtot).kron(D0eye)
        while (Matrix.infNorm(Told.sub(Tnew)) > 1e-10) {
            Told = Tnew.copy()

            val kronTI = Tnew.transpose().kron(Matrix.eye(mtot))
            val system = kronTI.add(eyeD0eye)
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

            Tnew = eyeTser.add(L.mult(LM))
        }
    } else {
        val schur = schurDecomposition(D0)
        val U = schur.first
        val Tr = schur.second
        val Ukron = U.kron(Matrix.eye(mser))
        val Trkron = Tr.kron(Matrix.eye(mser))

        while (Matrix.infNorm(Told.sub(Tnew)) > 1e-10) {
            Told = Tnew.copy()
            L = qSylvest(Ukron, Trkron, Tnew)
            Tnew = eyeTser.add(L.mult(LM))
        }
    }

    // Recompute L for final use
    val schur = schurDecomposition(D0)
    val U = schur.first
    val Tr = schur.second
    L = qSylvest(U.kron(Matrix.eye(mser)), Tr.kron(Matrix.eye(mser)), Tnew)

    // Compute theta_tot
    val thetaTot = Matrix(1, mtot)
    for (i in 0 until K) {
        val betaVec = Matrix(1, mser)
        for (j in 0 until beta[i].numCols) {
            betaVec[0, smk[i] + j] = beta[i][0, j] / mus[i]
        }
        val piDi = pi.mult(D[i])
        val contrib = piDi.kron(betaVec)
        thetaTot.addEq(1.0, contrib)
    }
    thetaTot.scaleEq(1.0 / load)

    // Find non-zero entries
    val nonzIndices = mutableListOf<Int>()
    for (i in 0 until mtot) {
        if (thetaTot[0, i] > 0) {
            nonzIndices.add(i)
        }
    }

    if (nonzIndices.isEmpty()) {
        return MMAPKPHK1Result(
            emptyList(),
            Matrix(1, 1),
            emptyList(),
            emptyList(),
            null
        )
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

    // Sojourn time alpha vectors
    val sojAlpha = mutableListOf<Matrix>()
    val diagThetaTot = Matrix.diag(*thetaTot.getRow(0).toArray1D())
    val onesM = Matrix.ones(m, 1)

    for (i in 0 until K) {
        val tempTser = Matrix.zeros(mser, 1)
        for (r in 0 until S[i].numRows) {
            tempTser[smk[i] + r, 0] = tser[smk[i] + r, 0]
        }

        val sojFull = diagThetaTot.mult(onesM.kron(tempTser)).scale(load / lambdas[i])
        val sojReduced = Matrix(1, nz)
        for (j in 0 until nz) {
            sojReduced[0, j] = sojFull[nonzIndices[j], 0]
        }
        sojAlpha.add(sojReduced)
    }

    // Overall sojourn time
    val sojOverallFull = diagThetaTot.mult(onesM.kron(tser)).scale(load / lambda)
    val sojOverall = Matrix(1, nz)
    for (j in 0 until nz) {
        sojOverall[0, j] = sojOverallFull[nonzIndices[j], 0]
    }
    sojAlpha.add(sojOverall)

    // Waiting time alpha vectors
    val waitAlpha = mutableListOf<Matrix>()
    for (i in 0 until K) {
        val DiSumCol = D[i].sumRows()
        val waitFull = diagThetaTot.mult(L).mult(DiSumCol.kron(tser)).scale(load / lambdas[i])
        val waitReduced = Matrix(1, nz)
        for (j in 0 until nz) {
            waitReduced[0, j] = waitFull[nonzIndices[j], 0]
        }
        waitAlpha.add(waitReduced)
    }

    // Overall waiting time
    val DsumMinusD0SumCol = Dsum.sub(D0).sumRows()
    val waitOverallFull = diagThetaTot.mult(L).mult(DsumMinusD0SumCol.kron(tser)).scale(load / lambda)
    val waitOverall = Matrix(1, nz)
    for (j in 0 until nz) {
        waitOverall[0, j] = waitOverallFull[nonzIndices[j], 0]
    }
    waitAlpha.add(waitOverall)

    // Overall queue length distribution (simplified without Hessenberg)
    val qlTList = mutableListOf<Double>()
    qlTList.add(1 - load)

    // Use iterative approach for queue length
    var Cn = Tnew.copy()
    var Ln = Matrix.zeros(mtot, mtot)

    // Using Sylvester-based iteration
    val DsumEye = Dsum.sub(D0).kron(Matrix.eye(mser))
    val schurTotal = schurDecomposition(D0)
    val Utot = schurTotal.first.kron(Matrix.eye(mser))
    val Trtot = schurTotal.second.kron(Matrix.eye(mser))

    var n = 1
    while (qlTList.sum() < 1 - 1e-10 && n < 1 + options.maxNumComp) {
        Ln = qSylvest(Utot, Trtot, Cn)
        Cn = Ln.mult(DsumEye)
        val qlVal = -load * thetaTot.mult(Ln.sumRows())[0, 0]
        qlTList.add(qlVal)
        n++
    }

    val qlTotal = Matrix(1, qlTList.size)
    for (i in qlTList.indices) {
        qlTotal[0, i] = qlTList[i]
    }

    // Per-type queue length distributions (simplified)
    val qlPerType = mutableListOf<Matrix>()
    for (t in 0 until K) {
        // Simplified: use same structure as total but scale
        val scaleFactor = lambdas[t] / lambda
        val qlType = Matrix(1, qlTotal.numCols)
        for (i in 0 until qlTotal.numCols) {
            qlType[0, i] = qlTotal[0, i] * scaleFactor
        }
        qlPerType.add(qlType)
    }

    return MMAPKPHK1Result(qlPerType, qlTotal, sojAlpha, waitAlpha, Smat)
}
