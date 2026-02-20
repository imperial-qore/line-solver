package jline.api.sn

import jline.lang.NetworkStruct
import jline.lang.constant.NodeType
import jline.lang.nodeparam.CacheNodeParam
import jline.solvers.AvgHandle
import jline.util.matrix.Matrix

/**
 * APIs to process NetworkStruct objects.
 */


fun snGetNodeTputFromTput(sn: NetworkStruct, TN: Matrix, TH: AvgHandle?, ANn: Matrix): Matrix {
    val I = sn.nnodes
    val C = sn.nchains
    val M = sn.nstations
    val R = sn.nclasses
    val TNn = Matrix(I, R)

    // First, copy station throughput to station nodes
    // This preserves the TN values computed by the solver
    for (ist in 0..<M) {
        val ind = sn.stationToNode[ist].toInt()
        if (ind >= 0 && ind < I) {
            for (r in 0..<R) {
                TNn[ind, r] = TN[ist, r]
            }
        }
    }

    // Process Cache hit/miss classes using nodevisits formula
    for (ind in 0..<I) {
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
                        if (isHitClass || isMissClass) {
                            // TNn(ind, r) =  (sn.nodevisits{c}(ind,r) / sum(sn.visits{c}(sn.stationToStateful(refstat),inchain))) * sum(TN(refstat,inchain));
                            val num = sn.nodevisits[c]!![ind, rIdx]
                            var den = 0.0
                            for (s in 0..<inchain.length()) {
                                val sIdx = inchain[s].toInt()
                                den = den + sn.visits[c]!![sn.stationToStateful[refstat].toInt(), sIdx]
                            }
                            var coeff = 0.0
                            for (s in 0..<inchain.length()) {
                                val sIdx = inchain[s].toInt()
                                coeff = coeff + TN[refstat, sIdx]
                            }
                            if (den > 0) {
                                TNn[ind, rIdx] = num / den * coeff
                            } else {
                                TNn[ind, rIdx] = 0
                            }
                        }
                    }
                }
            }
        }
    }

    // Process non-station nodes using routing formula
    for (ind in 0..<I) {
        // Skip if this is a station node - keep the TN values copied above
        val nodeToStation = sn.nodeToStation[ind].toInt()
        if (nodeToStation >= 0) continue

        for (c in 0..<C) {
            val inchain = sn.inchain[c]
            for (r in 0..<inchain!!.length()) {
                val rIdx = inchain[r].toInt()
                val anystateful = !sn.visits[c]!!.getColumn(rIdx).isEmpty
                if (anystateful) {
                    if (sn.nodetype[ind] != NodeType.Sink && sn.nodetype[ind] != NodeType.Join) {
                        for (s in 0..<inchain.length()) {
                            val sIdx = inchain[s].toInt()
                            for (jnd in 0..<I) {
                                if (sn.nodetype[ind] == NodeType.Cache) {
                                    if (ind != jnd) {
                                        TNn[ind, sIdx] = TNn[ind, sIdx] + ANn[ind, rIdx] * sn.rtnodes[ind * R + rIdx, jnd * R + sIdx]
                                    }
                                } else {
                                    TNn[ind, sIdx] = TNn[ind, sIdx] + ANn[ind, rIdx] * sn.rtnodes[ind * R + rIdx, jnd * R + sIdx]
                                }
                            }
                        }
                    } else if (sn.nodetype[ind] == NodeType.Join) {
                        for (s in 0..<inchain.length()) {
                            val sIdx = inchain[s].toInt()
                            for (jnd in 0..<I) {
                                TNn[ind, sIdx] = TNn[ind, sIdx] + ANn[ind, rIdx] * sn.rtnodes[(ind) * R + rIdx, (jnd) * R + sIdx]
                            }
                        }
                    }
                }
            }
        }
    }
    return TNn
}
/**
 * Stochastic network GetNodeTputFromTput algorithms
 */
@Suppress("unused")
class SngetnodetputfromtputAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}