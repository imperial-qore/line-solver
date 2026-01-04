/**
 * @file Product-form queueing network parameter extraction
 * 
 * Extracts essential parameters (service demands, populations, visit ratios) from
 * a queueing network model for product-form analysis algorithms like MVA, convolution,
 * and mean value analysis. Transforms network topology into algorithmic parameters.
 *
 * @since LINE 3.0
 */
package jline.api.sn

import jline.io.Ret
import jline.lang.NetworkStruct
import jline.lang.constant.NodeType
import jline.util.Maths
import jline.util.Utils
import jline.util.matrix.Matrix
import kotlin.math.ceil

/**
 * Calculate the parameters at class level for a queueing network model
 *
 * @param sn NetworkStruct object for the queueing network model.
 * @return queueing network parameters
 */
fun snGetProductFormParams(sn: NetworkStruct): Ret.snGetProductFormParams {
    val R = sn.nclasses
    val N = sn.njobs
    val queueIndices = ArrayList<Int>()
    val delayIndices = ArrayList<Int>()
    var sourceIndex = -1
    for (i in sn.nodetype.indices) {
        if (sn.nodetype[i] == NodeType.Queue) {
            queueIndices.add(i)
        } else if (sn.nodetype[i] == NodeType.Delay) {
            delayIndices.add(i)
        } else if (sn.nodetype[i] == NodeType.Source) {
            sourceIndex = i
        }
    }
    val Mq = queueIndices.size
    val Mz = delayIndices.size

    val lambda = Matrix(1, R)
    for (r in 0..<R) {
        if (Utils.isInf(N[r])) {
            lambda[0, r] = sn.rates[sn.nodeToStation[sourceIndex].toInt(), r]
        }
    }
    val S = Matrix(queueIndices.size, 1)
    var Smax = Double.MIN_VALUE
    for (i in queueIndices.indices) {
        S[i, 0] = sn.nservers[sn.nodeToStation[queueIndices[i]].toInt()]
        if (java.lang.Double.isFinite(S[i, 0]) && S[i, 0] > Smax) Smax = S[i, 0]
    }
    val D = Matrix(Mq, R)
    var Nct = 0.0
    for (i in 0..<N.numRows) {
        for (j in 0..<N.numCols) {
            if (java.lang.Double.isFinite(N[i, j])) Nct += N[i, j]
        }
    }
    val mu = Matrix.ones(Mq, (ceil(Nct) + Smax).toInt())
    for (ist in 0..<Mq) {
        for (r in 0..<R) {
            var c = 0
            while (c < sn.chains.numRows) {
                if (sn.chains[c, r] != 0.0) break
                c++
            }
            if (sn.refclass[c] > 0) {
                D[ist, r] =
                    sn.visits[c]!![sn.nodeToStateful[queueIndices[ist]].toInt(), r] / sn.rates[sn.nodeToStation[queueIndices[ist]].toInt(), r] / sn.visits[c]!![sn.stationToStateful[sn.refstat[r].toInt()].toInt(), sn.refclass[c].toInt()]
            } else {
                D[ist, r] =
                    sn.visits[c]!![sn.nodeToStateful[queueIndices[ist]].toInt(), r] / sn.rates[sn.nodeToStation[queueIndices[ist]].toInt(), r]
            }
        }
        for (j in 0..<mu.numCols) {
            mu[ist, j] = Maths.min((j + 1).toDouble(), sn.nservers[sn.nodeToStation[queueIndices[ist]].toInt()])
        }
    }
    val Z = Matrix(Maths.max(1.0, Mz.toDouble()).toInt(), R)
    for (ist in 0..<Mz) {
        for (r in 0..<R) {
            var c = 0
            while (c < sn.chains.numRows) {
                if (sn.chains[c, r] != 0.0) break
                c++
            }
            if (sn.refclass[c] > 0) {
                Z[ist, r] =
                    sn.visits[c]!![sn.nodeToStateful[delayIndices[ist]].toInt(), r] / sn.rates[sn.nodeToStation[delayIndices[ist]].toInt(), r] / sn.visits[c]!![sn.stationToStateful[sn.refstat[r].toInt()].toInt(), sn.refclass[c].toInt()]
            } else {
                Z[ist, r] =
                    sn.visits[c]!![sn.nodeToStateful[delayIndices[ist]].toInt(), r] / sn.rates[sn.nodeToStation[delayIndices[ist]].toInt(), r]
            }
        }
    }
    var V: Matrix? = null
    for (i in sn.visits.keys) {
        val visitsValue = sn.visits[i]
        V = if (V == null) {
            visitsValue
        } else {
            V.add(1.0, visitsValue)
        }
    }
    D.apply(Double.NaN, 0.0, "equal")
    Z.apply(Double.NaN, 0.0, "equal")
    val ret = Ret.snGetProductFormParams()
    ret.lambda = lambda
    ret.D = D
    ret.N = N
    ret.Z = Z
    ret.mu = mu
    ret.S = S
    ret.V = V
    return ret
}
/**
 * Stochastic network GetProductFormParams algorithms
 */
@Suppress("unused")
class SngetproductformparamsAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}