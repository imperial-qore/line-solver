package jline.api.sn

import jline.lang.NetworkStruct
import jline.lang.constant.NodeType
import jline.solvers.AvgHandle
import jline.util.matrix.Matrix

/**
 * Calculates the average arrival rates at each station from the network throughputs.
 *
 * @param sn the NetworkStruct object for the queueing network model
 * @param TN the matrix of throughputs
 * @param TH throughput handles for the solver
 * @return a matrix of average arrival rates
 */
fun snGetArvRFromTput(sn: NetworkStruct, TN: Matrix?, TH: AvgHandle?): Matrix {
    val M = sn.nstations
    val R = sn.nclasses
    var AN = Matrix(M, R)

    if (TH != null && TN != null) {
        for (ist in 0..<M) {
            if (sn.nodetype[sn.stationToNode[ist].toInt()] == NodeType.Source) continue  // same as AN(ist,:)=0.

            for (k in 0..<R) {
                var a = 0.0
                for (jst in 0..<M) {
                    for (r in 0..<R) {
                        a += TN[jst, r] * sn.rt[(jst * R) + r, (ist * R) + k]
                    }
                }
                AN[ist, k] = a
            }
        }
    } else {
        AN = Matrix(0, 0)
    }

    if (sn.fj.any()) {
        val ANn = snGetNodeArvRFromTput(sn, TN!!, TH, AN)
        for (ist in 0..<M) {
            for (r in 0..<R) {
                AN[ist, r] = ANn[sn.stationToNode[ist].toInt(), r]
            }
        }
    }

    return AN
}
/**
 * Stochastic network GetArvRFromTput algorithms
 */
@Suppress("unused")
class SngetarvrfromtputAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}