package jline.api.sn

import jline.lang.NetworkStruct
import jline.lang.constant.NodeType
import jline.lang.nodeparam.CacheNodeParam
import jline.solvers.AvgHandle
import jline.util.matrix.Matrix

/**
 * Calculates the average arrival rates at each station from the network throughputs.
 *
 * This function computes arrival rates using the routing matrix (sn.rt) which is indexed
 * by stateful nodes, not stations. For Cache nodes (stateful but not stations), throughputs
 * are computed from hit/miss probabilities.
 *
 * @param sn the NetworkStruct object for the queueing network model
 * @param TN the matrix of throughputs (indexed by stations)
 * @param TH throughput handles for the solver
 * @return a matrix of average arrival rates (indexed by stations)
 */
fun snGetArvRFromTput(sn: NetworkStruct, TN: Matrix?, TH: AvgHandle?): Matrix {
    val M = sn.nstations
    val R = sn.nclasses
    var AN = Matrix(M, R)

    if (TH != null && TN != null) {
        // Build list of stateful node indices
        val statefulNodes = mutableListOf<Int>()
        for (i in 0..<sn.isstateful.length()) {
            if (sn.isstateful[i] == 1.0) {
                statefulNodes.add(i)
            }
        }
        val nStateful = statefulNodes.size

        // Build throughput vector for all stateful nodes
        // Stations have TN directly; non-station stateful nodes need throughput computed
        val TN_stateful = Matrix(nStateful, R)

        // First pass: set throughputs for stations
        for (sf in 0..<nStateful) {
            val ind = statefulNodes[sf]
            val ist = sn.nodeToStation[ind].toInt()
            if (ist >= 0) {
                // This stateful node is a station - use station throughput
                for (r in 0..<R) {
                    TN_stateful[sf, r] = TN[ist, r]
                }
            }
        }

        // Second pass: compute throughputs for non-station stateful nodes (e.g., Cache)
        for (sf in 0..<nStateful) {
            val ind = statefulNodes[sf]
            val ist = sn.nodeToStation[ind].toInt()
            if (ist < 0) {
                // This is a non-station stateful node
                if (sn.nodetype[ind] == NodeType.Cache) {
                    // For Cache nodes, compute outgoing throughputs using hit/miss probabilities
                    // Use sn.stateful[sf] to get the correct StatefulNode key for nodeparam
                    val statefulNode = sn.stateful[sf]
                    val cacheParam = sn.nodeparam[statefulNode] as? CacheNodeParam
                    if (cacheParam != null) {
                        val hitclass = cacheParam.hitclass
                        val missclass = cacheParam.missclass
                        val actualHitProb = cacheParam.actualhitprob
                        val actualMissProb = cacheParam.actualmissprob

                        if (actualHitProb != null && actualMissProb != null && !actualHitProb.isEmpty) {
                            // Compute throughput for each chain
                            for (c in 0..<sn.nchains) {
                                val inchain = sn.inchain[c] ?: continue
                                val refstat = sn.refstat[c].toInt()

                                // Get total throughput from reference station for this chain
                                var totalTput = 0.0
                                for (idx in 0..<inchain.length()) {
                                    val classIdx = inchain[idx]?.toInt() ?: continue
                                    if (classIdx in 0..<R) {
                                        totalTput += TN[refstat, classIdx]
                                    }
                                }

                                // Set throughput for hit/miss classes based on hit/miss probabilities
                                // JAR hitclass/missclass are 0-indexed (-1 = invalid)
                                // MATLAB hitclass/missclass are 1-indexed (0 = invalid)
                                for (k in 0..<R) {
                                    if (hitclass.length() > k && missclass.length() > k) {
                                        val hc = hitclass[k].toInt()
                                        val mc = missclass[k].toInt()
                                        if (hc >= 0 && hc < R && actualHitProb.length() > k) {
                                            val hitProbVal = actualHitProb[k]
                                            if (!hitProbVal.isNaN()) {
                                                TN_stateful[sf, hc] = totalTput * hitProbVal
                                            }
                                        }
                                        if (mc >= 0 && mc < R && actualMissProb.length() > k) {
                                            val missProbVal = actualMissProb[k]
                                            if (!missProbVal.isNaN()) {
                                                TN_stateful[sf, mc] = totalTput * missProbVal
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        // Propagate throughputs to non-station, non-Cache stateful nodes
        // (e.g., Router) using already-known station throughputs and the rt matrix
        for (sf in 0..<nStateful) {
            val ind = statefulNodes[sf]
            val ist = sn.nodeToStation[ind].toInt()
            if (ist < 0 && sn.nodetype[ind] != NodeType.Cache) {
                // Compute throughput from upstream stateful nodes
                for (sf_src in 0..<nStateful) {
                    for (k in 0..<R) {
                        for (r in 0..<R) {
                            TN_stateful[sf, k] = TN_stateful[sf, k] + TN_stateful[sf_src, r] * sn.rt[sf_src * R + r, sf * R + k]
                        }
                    }
                }
            }
        }

        // Compute arrival rates using stateful node throughputs and rt matrix
        for (ist in 0..<M) {
            val ind_ist = sn.stationToNode[ist].toInt()
            if (sn.nodetype[ind_ist] == NodeType.Source) continue  // AN(ist,:)=0

            // Find the position of this station in the stateful node list
            val sf_ist = statefulNodes.indexOf(ind_ist)
            if (sf_ist < 0) continue

            for (k in 0..<R) {
                var a = 0.0
                for (sf_jst in 0..<nStateful) {
                    for (r in 0..<R) {
                        a += TN_stateful[sf_jst, r] * sn.rt[sf_jst * R + r, sf_ist * R + k]
                    }
                }
                AN[ist, k] = a
            }
        }
    } else {
        AN = Matrix(0, 0)
    }

    if (sn.fj.any()) {
        // Apply node-level arrival rates (for Join nodes, etc.)
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