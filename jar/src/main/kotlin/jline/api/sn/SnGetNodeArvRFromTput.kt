package jline.api.sn

import jline.lang.NetworkStruct
import jline.lang.constant.NodeType
import jline.lang.nodeparam.CacheNodeParam
import jline.solvers.AvgHandle
import jline.util.matrix.Matrix

fun snGetNodeArvRFromTput(sn: NetworkStruct, TN: Matrix, TH: AvgHandle?, AN: Matrix?): Matrix {
    val I = sn.nnodes
    val C = sn.nchains
    sn.nstations
    val R = sn.nclasses
    val ANn = Matrix(I, R)

    // Debug: Log input parameters
    val debug = false
    if (debug) {
        println("DEBUG snGetNodeArvRFromTput: I=$I, C=$C, R=$R")
        println("DEBUG TN (station throughputs):")
        for (i in 0..<TN.numRows) {
            for (j in 0..<TN.numCols) {
                print("  TN[$i,$j]=${TN[i,j]}")
            }
            println()
        }
        if (AN != null && AN.numRows > 0) {
            println("DEBUG AN (station arrival rates):")
            for (i in 0..<AN.numRows) {
                for (j in 0..<AN.numCols) {
                    print("  AN[$i,$j]=${AN[i,j]}")
                }
                println()
            }
        }
        println("DEBUG nodevisits[0] (chain 0):")
        if (sn.nodevisits.containsKey(0)) {
            val nv = sn.nodevisits[0]!!
            for (i in 0..<nv.numRows) {
                for (j in 0..<nv.numCols) {
                    print("  nv[$i,$j]=${nv[i,j]}")
                }
                println()
            }
        }
    }

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
                        if (debug && ind == 0 && rIdx == 0) {
                            println("DEBUG node[$ind] class[$rIdx]: num=$num, den=$den, coeff=$coeff, ANn=${ANn[ind, rIdx]}")
                        }
                    }
                }
            }
        }
    }
    if (debug) {
        println("DEBUG ANn (result):")
        for (i in 0..<ANn.numRows) {
            for (j in 0..<ANn.numCols) {
                print("  ANn[$i,$j]=${ANn[i,j]}")
            }
            println()
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
