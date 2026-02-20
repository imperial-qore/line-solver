/**
 * @file Quasi-Birth-Death process setup delays and server switch-off analysis
 *
 * Analyzes queueing systems with server setup delays and switch-off mechanisms.
 * Models energy-efficient server operations and startup costs in queueing systems.
 *
 * Ported from MATLAB's qbd_setupdelayoff.m for exact parity.
 *
 * @since LINE 3.0
 */
package jline.api.mam

import jline.lang.processes.APH
import jline.lib.smc.QBD_CR
import jline.lib.smc.QBD_pi
import jline.util.matrix.Matrix

/**
 * Analyze a queue with setup delays and server switch-off using QBD approach.
 *
 * Models an M/M/1 queue where:
 * - When server becomes idle, it enters a delayoff phase before turning off
 * - When a job arrives to an off server, it goes through setup before service
 *
 * Uses APH distributions to model general setup and delayoff times, matching
 * MATLAB's qbd_setupdelayoff.m implementation exactly.
 *
 * @param lambda Arrival rate
 * @param mu Service rate
 * @param alphaRate Setup rate (1/mean setup time)
 * @param alphaScv Setup squared coefficient of variation (1.0 for exponential)
 * @param betaRate Delayoff rate (1/mean delayoff time)
 * @param betaScv Delayoff squared coefficient of variation (1.0 for exponential)
 * @return Mean queue length
 */
fun qbd_setupdelayoff(
    lambda: Double,
    mu: Double,
    alphaRate: Double,
    alphaScv: Double,
    betaRate: Double,
    betaScv: Double
): Double {
    // Fit APH distributions for setup and delayoff times
    // MATLAB: alpha = APH.fitMeanAndSCV(1/alpharate, alphascv).getProcess;
    val alphaAPH = APH.fitMeanAndSCV(1.0 / alphaRate, alphaScv)
    val alphaT = alphaAPH.getParam(3).value as Matrix  // Generator matrix T
    val na = alphaT.numRows  // Number of setup phases

    // MATLAB: beta = APH.fitMeanAndSCV(1/betarate, betascv).getProcess;
    val betaAPH = APH.fitMeanAndSCV(1.0 / betaRate, betaScv)
    val betaT = betaAPH.getParam(3).value as Matrix  // Generator matrix T
    val nb = betaT.numRows  // Number of delayoff phases

    val n = na + nb  // Total number of phases

    // Build QBD matrices matching MATLAB exactly
    // F = forward transitions (arrivals)
    val F = Matrix.zeros(n, n)
    // MATLAB: for i=1:na, F(i,i) = lambda; end
    for (i in 0 until na) {
        F[i, i] = lambda
    }
    // MATLAB: for i=1:nb, F(na+i,na+1) = lambda; end
    for (i in 0 until nb) {
        F[na + i, na] = lambda  // na+1 in 1-based = na in 0-based
    }
    // MATLAB: F(na+1,na+1) = lambda; (already set by previous loop when i=0)
    // Actually it's redundant since F(na+1,na+1) = lambda is set when i=0 above
    // But we also have F(na,na) = lambda from the first loop, let's double-check MATLAB:
    // F(na+1,na+1) appears to be set twice - once in first loop (i=na) and once explicitly
    // Actually first loop is 1:na so i goes up to na, which in 0-based is na-1
    // So F(na+1,na+1) in MATLAB = F[na,na] in 0-based is NOT set by first loop
    // It's set by the explicit F(na+1,na+1) = lambda line
    F[na, na] = lambda

    // B = backward transitions (service completions)
    val B = Matrix.zeros(n, n)
    // MATLAB: B(na+1,na+1) = mu;
    B[na, na] = mu

    // L = local transitions at levels > 0
    val L = Matrix.zeros(n, n)
    // MATLAB: for i=1:na
    //           L(i,i) = alpha{1}(i,i) - lambda;
    //           if i<na
    //             L(i,(i+1):na) = alpha{1}(i,(i+1):na);
    //           else
    //             L(na,na+1) = -alpha{1}(end,end);
    //           end
    //         end
    for (i in 0 until na) {
        L[i, i] = alphaT[i, i] - lambda
        if (i < na - 1) {
            // Copy off-diagonal elements from alpha generator
            for (j in i + 1 until na) {
                L[i, j] = alphaT[i, j]
            }
        } else {
            // Last setup phase: exit rate to busy phase
            // MATLAB: L(na,na+1) = -alpha{1}(end,end)
            L[na - 1, na] = -alphaT[na - 1, na - 1]
        }
    }
    // MATLAB: L(na+1,na+1) = -mu - lambda;
    L[na, na] = -mu - lambda
    // MATLAB: for i=2:nb, L(na+i,na+i) = -lambda; end
    for (i in 1 until nb) {
        L[na + i, na + i] = -lambda
    }

    // L0 = local transitions at level 0 (boundary)
    val L0 = Matrix.zeros(n, n)
    // MATLAB: for i=1:na, L0(i,i) = -lambda; end
    for (i in 0 until na) {
        L0[i, i] = -lambda
    }
    // MATLAB: for i=1:nb
    //           L0(na+i,na+i) = beta{1}(i,i) - lambda;
    //           if i==nb
    //             L0(na+i,1) = -beta{1}(i,i);
    //           else
    //             L0(na+i,na+i+1) = -beta{1}(i,i);
    //           end
    //         end
    for (i in 0 until nb) {
        L0[na + i, na + i] = betaT[i, i] - lambda
        if (i == nb - 1) {
            // Last delayoff phase: exit to idle (phase 0)
            L0[na + i, 0] = -betaT[i, i]
        } else {
            // Transition to next delayoff phase
            L0[na + i, na + i + 1] = -betaT[i, i]
        }
    }

    // Solve QBD using Cyclic Reduction
    // MATLAB: [~,R,~] = QBD_CR(B,L,F);
    val qbdResult = QBD_CR(B, L, F, null, null, null, null)
    val R = qbdResult["R"] ?: throw RuntimeException("QBD_CR failed to compute R matrix")

    // Compute stationary distribution
    // MATLAB: pn = QBD_pi(B,L0,R);
    val pn = QBD_pi(B, L0, R)

    // Compute mean queue length
    // MATLAB: QN = 0; j = n+1; ni = 0;
    //         while 1
    //           ni = ni + 1;
    //           QN = QN + ni*sum(pn(j:(j+n)));
    //           j = j + n;
    //           if j+n > length(pn), break; end
    //         end
    var QN = 0.0
    var j = n  // j = n+1 in 1-based = n in 0-based
    var ni = 0
    val pnLength = pn.numCols  // pn is a row vector

    while (true) {
        ni++
        // Sum pn(j:(j+n)) - MATLAB range j:(j+n) is n+1 elements (inclusive both ends)
        // In 0-indexed: j to j+n inclusive = indices j, j+1, ..., j+n
        var levelSum = 0.0
        val endIdx = minOf(j + n + 1, pnLength)  // j+n+1 exclusive = j+n inclusive
        for (k in j until endIdx) {
            levelSum += pn[0, k]
        }
        QN += ni * levelSum
        j += n  // MATLAB: j = j + n
        if (j + n > pnLength) {  // MATLAB: if j+n > length(pn)
            break
        }
    }

    return QN
}

/**
 * QBD setupdelayoff algorithms
 */
@Suppress("unused")
class QbdSetupdelayoff {
    companion object {
        // Class documentation marker for Dokka
    }
}
