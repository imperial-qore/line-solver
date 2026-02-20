package jline.api.sn

import jline.io.Ret
import jline.lang.NetworkStruct
import jline.lang.constant.NodeType
import jline.util.matrix.Matrix

/**
 * Calculate the parameters at class and chain level for a queueing network model
 *
 * @param sn - NetworkStruct object for the queueing network model.
 * @return queueing network parameters
 */
fun snGetProductFormChainParams(sn: NetworkStruct): Ret.snGetProductFormParams {
    val ret1 = snGetProductFormParams(sn)
    val lambda = ret1.lambda
    val mu = ret1.mu
    val queueIndex = ArrayList<Int>()
    val delayIndex = ArrayList<Int>()
    val ignoreIndex = HashSet<Int>()
    for (i in sn.nodetype.indices) {
        when (sn.nodetype[i]) {
            NodeType.Queue -> queueIndex.add(i)
            NodeType.Delay -> delayIndex.add(i)
            NodeType.Source, NodeType.Join -> ignoreIndex.add(sn.nodeToStation[i].toInt())
            else -> {}
        }
    }
    val ret2 = snGetDemandsChain(sn)
    val Dchain = ret2.Dchain
    val Vchain = ret2.Vchain.copy()
    val Nchain = ret2.Nchain
    val lambda_chains = Matrix(1, sn.nchains)

    val D_chains = Matrix(queueIndex.size, sn.nchains)
    val Z_chains = Matrix(delayIndex.size, sn.nchains)

    for (c in 0..<sn.nchains) {
        var lambdaSum = 0.0
        val inc = sn.inchain[c]
        for (i in 0..<inc!!.numCols) {
            if (!java.lang.Double.isNaN(lambda[inc[i].toInt()])) {
                lambdaSum += lambda[inc[i].toInt()]
            }
        }
        lambda_chains[0, c] = lambdaSum
        for (i in queueIndex.indices) {
            D_chains[i, c] = Dchain[sn.nodeToStation[queueIndex[i]].toInt(), c]
        }
        for (i in delayIndex.indices) {
            Z_chains[i, c] = Dchain[sn.nodeToStation[delayIndex[i]].toInt(), c]
        }
    }
    val S = Matrix(queueIndex.size, 1)
    for (i in queueIndex.indices) {
        S[i, 0] = sn.nservers[sn.nodeToStation[queueIndex[i]].toInt()]
    }
    val ret = Ret.snGetProductFormParams()
    ret.S = S
    ret.lambda = lambda_chains
    ret.N = Nchain
    Vchain.removeRows(ignoreIndex)
    ret.D = D_chains
    ret.Z = Z_chains
    ret.mu = mu
    ret.V = Vchain
    if (ret.Z.isEmpty) {
        ret.Z = Matrix(ret.N.numRows, ret.N.numCols)
    }
    return ret
}
/**
 * Stochastic network GetProductFormChainParams algorithms
 */
@Suppress("unused")
class SngetproductformchainparamsAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}