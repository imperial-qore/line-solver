package jline.api.sn

import jline.GlobalConstants.Inf
import jline.io.Ret
import jline.lang.NetworkStruct
import jline.lang.nodes.Sink
import jline.util.Utils
import jline.util.matrix.Matrix

/**
 * Calculate class-based performance metrics for a queueing network based on performance measures of its chains
 *
 * @param sn      - NetworkStruct object for the queueing network model
 * @param Lchain  - service demands per chain
 * @param ST      - mean service times per class
 * @param STchain - mean service times per chain
 * @param Vchain  - mean visits per chain
 * @param alpha   - class aggregation coefficients
 * @param Qchain  - mean queue-lengths per chain
 * @param Uchain  - mean utilization per chain
 * @param Rchain  - mean response time per chain
 * @param Tchain  - mean throughput per chain
 * @param Cchain  - mean system response time per chain
 * @param Xchain  - mean system throughput per chain
 * @return chain performance metrics
 */

fun snDeaggregateChainResults(sn: NetworkStruct,
                              Lchain: Matrix,
                              ST: Matrix?,
                              STchain: Matrix,
                              Vchain: Matrix,
                              alpha: Matrix,
                              Qchain: Matrix?,
                              Uchain: Matrix?,
                              Rchain: Matrix,
                              Tchain: Matrix,
                              Cchain: Matrix?,
                              Xchain: Matrix): Ret.snDeaggregateChainResults {
    var ST = ST
    if (ST == null || ST.isEmpty) {
        ST = Matrix(0, 0)
        sn.rates.divide(1.0, ST, false)
        ST.removeNaN()
    }

    if (Cchain != null && !Cchain.isEmpty) throw RuntimeException("Cchain input to snDeaggregateChainResults not yet supported")

    val M = sn.nstations
    val K = sn.nclasses
    val X = Matrix(1, K)
    val U = Matrix(M, K)
    val Q = Matrix(M, K)
    val T = Matrix(M, K)
    val R = Matrix(M, K)
    val C = Matrix(1, K)

    var idxSink = 0
    for (nodeIter in sn.nodes) {
        if (nodeIter is Sink) idxSink = nodeIter.getNodeIndex()
    }

    var Vsinktmp = Matrix(sn.nodevisits[0]!!.numRows, sn.nodevisits[0]!!.numCols)
    for (i in 0..<sn.nodevisits.size) {
        Vsinktmp = Vsinktmp.add(1.0, sn.nodevisits[i])
    }
    val Vsink = Matrix.extractRows(Vsinktmp, idxSink, idxSink + 1, null)
    for (c in 0..<sn.nchains) {
        val inchain_c = sn.inchain[c]
        var sum = 0.0
        for (idx in 0..<inchain_c!!.numCols) sum += sn.njobs[inchain_c[idx].toInt()]
        for (idx in 0..<inchain_c.numCols) {
            val k = inchain_c[0, idx].toInt()
            if (Utils.isInf(sum)) X[0, k] = Xchain[0, c] * Vsink[0, k]
            else X[0, k] = Xchain[0, c] * alpha[sn.refstat[k, 0].toInt(), k]
            for (i in 0..<M) {
                if (Uchain == null || Uchain.isEmpty) {
                    if (Utils.isInf(sn.nservers[i, 0])) U[i, k] =
                        ST[i, k] * (Xchain[0, c] * Vchain[i, c] / Vchain[sn.refstat[k, 0].toInt(), c]) * alpha[i, k]
                    else U[i, k] =
                        ST[i, k] * (Xchain[0, c] * Vchain[i, c] / Vchain[sn.refstat[k, 0].toInt(), c]) * alpha[i, k] / sn.nservers[i, 0]
                } else {
                    if (Utils.isInf(sn.nservers[i, 0])) U[i, k] =
                        ST[i, k] * (Xchain[0, c] * Vchain[i, c] / Vchain[sn.refstat[k, 0].toInt(), c]) * alpha[i, k]
                    else U[i, k] = Uchain[i, c] * alpha[i, k]
                }

                if (Lchain[i, c] > 0) {
                    if (Qchain != null && !Qchain.isEmpty) Q[i, k] = Qchain[i, c] * alpha[i, k]
                    else Q[i, k] =
                        Rchain[i, c] * ST[i, k] / STchain[i, c] * Xchain[0, c] * Vchain[i, c] / Vchain[sn.refstat[k, 0].toInt(), c] * alpha[i, k]
                    T[i, k] = Tchain[i, c] * alpha[i, k]
                    R[i, k] = Q[i, k] / T[i, k]
                } else {
                    T.remove(i, k)
                    R.remove(i, k)
                    Q.remove(i, k)
                }
            }
            C[0, k] = sn.njobs[0, k] / X[0, k]
        }
    }

    Q.absEq()
    R.absEq()
    X.absEq()
    U.absEq()
    T.absEq()
    C.absEq()
    Q.removeNaN()
    Q.apply(Inf, 0.0, "equal")
    R.removeNaN()
    R.apply(Inf, 0.0, "equal")
    X.removeNaN()
    X.apply(Inf, 0.0, "equal")
    U.removeNaN()
    U.apply(Inf, 0.0, "equal")
    T.removeNaN()
    T.apply(Inf, 0.0, "equal")
    C.removeNaN()
    C.apply(Inf, 0.0, "equal")

    return Ret.snDeaggregateChainResults(Q, U, R, T, C, X)
}
/**
 * Stochastic network DeaggregateChainResults algorithms
 */
@Suppress("unused")
class SndeaggregatechainresultsAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}