package jline.api.sn

import jline.lang.NetworkStruct
import jline.lang.constant.NodeType
import jline.lang.nodeparam.CacheNodeParam
import jline.solvers.AvgHandle
import jline.util.matrix.Matrix

fun snGetNodeArvRFromTput(sn: NetworkStruct, TN: Matrix, TH: AvgHandle?, AN: Matrix?): Matrix {
    val I = sn.nnodes
    val C = sn.nchains
    val M = sn.nstations
    val R = sn.nclasses
    val ANn = Matrix(I, R)

    // First, copy station arrival rates to station nodes (MATLAB lines 48-51)
    // This preserves the AN values computed in snGetArvRFromTput (including fork scaling)
    if (AN != null && AN.numRows > 0) {
        for (ist in 0..<M) {
            val ind = sn.stationToNode[ist].toInt()
            if (ind >= 0 && ind < I) {
                for (r in 0..<R) {
                    ANn[ind, r] = AN[ist, r]
                }
            }
        }
    }

    // Process non-station nodes using nodevisits formula (MATLAB lines 52-81)
    for (ind in 0..<I) {
        // Skip if this is a station node - keep the AN values copied above
        val nodeToStation = sn.nodeToStation[ind].toInt()
        if (nodeToStation >= 0) continue

        for (c in 0..<C) {
            val inchain = sn.inchain[c]
            val refstat = sn.refstat[c].toInt()
            for (r in 0..<inchain!!.length()) {
                val rIdx = inchain[r].toInt()
                if (sn.nodetype[ind] != NodeType.Source) {
                    if (sn.nodetype[ind] == NodeType.Cache) {
                        val np = sn.nodeparam[sn.nodes[ind]] as CacheNodeParam?
                        val isHitClass = !np!!.hitclass.findNumber(rIdx.toDouble()).isEmpty
                        val isMissClass = !np.missclass.findNumber(rIdx.toDouble()).isEmpty
                        if (!(isHitClass || isMissClass)) {
                            // ANn(ind, r) =  (sn.nodevisits{c}(ind,r) / sum(sn.visits{c}(sn.stationToStateful(refstat),inchain))) * sum(TN(refstat,inchain));
                            val num = sn.nodevisits[c]!![ind, rIdx]
                            var den = 0.0
                            for (r2 in 0..<inchain.length()) {
                                val rprime = inchain[r2].toInt()
                                den = den + sn.visits[c]!![sn.stationToStateful[refstat].toInt(), rprime]
                            }
                            var coeff = 0.0
                            for (r2 in 0..<inchain.length()) {
                                val rprime = inchain[r2].toInt()
                                coeff = coeff + TN[refstat, rprime]
                            }
                            if (den > 0) {
                                ANn[ind, rIdx] = num / den * coeff
                            } else {
                                ANn[ind, rIdx] = 0
                            }
                        }
                    } else {
                        // Non-station, non-Cache nodes (e.g., Sink, ClassSwitch)
                        // ANn(ind, r) =  (sn.nodevisits{c}(ind,r) / sum(sn.visits{c}(sn.stationToStateful(refstat),inchain))) * sum(TN(refstat,inchain));
                        val num = sn.nodevisits[c]!![ind, rIdx]
                        var den = 0.0
                        for (r2 in 0..<inchain.length()) {
                            val rprime = inchain[r2].toInt()
                            den = den + sn.visits[c]!![sn.stationToStateful[refstat].toInt(), rprime]
                        }
                        var coeff = 0.0
                        for (r2 in 0..<inchain.length()) {
                            val rprime = inchain[r2].toInt()
                            coeff = coeff + TN[refstat, rprime]
                        }
                        if (den > 0) {
                            ANn[ind, rIdx] = num / den * coeff
                        } else {
                            ANn[ind, rIdx] = 0
                        }
                    }
                }
            }
        }
    }
    return ANn
}
/**
 * Stochastic network GetNodeArvRFromTput algorithms
 */
@Suppress("unused")
class SngetnodearvrfromtputAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}
