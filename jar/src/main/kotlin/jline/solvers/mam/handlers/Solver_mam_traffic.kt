package jline.solvers.mam.handlers

import jline.api.mam.mmap_lambda
import jline.api.mam.mmap_normalize
import jline.api.mam.mmap_super
import jline.api.mc.dtmc_stochcomp
import jline.api.npfqn.npfqn_traffic_merge
import jline.api.npfqn.npfqn_traffic_split_cs
import jline.lang.NetworkStruct
import jline.GlobalConstants
import jline.lang.constant.NodeType
import jline.solvers.SolverOptions
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

fun solver_mam_traffic(sn: NetworkStruct,
                       DEP: MutableMap<Int?, MutableMap<Int?, MatrixCell?>?>,
                       config: SolverOptions.Config): MutableMap<Int?, MatrixCell?> {
    val I = sn.nnodes
    sn.nchains
    val R = sn.nclasses
    val non_cs_classes = Matrix(1, I * R, I * R)
    val isNCS = Matrix(1, I, I)
    val nodeToNCS = Matrix(1, I, I)
    var end = 0
    for (ind in 0..<I) {
        if (sn.nodetype.get(ind) != NodeType.ClassSwitch) {
            for (i in 0..<R) {
                non_cs_classes.set(i + end, (ind * R + i).toDouble())
            }
            end = end + R
            isNCS.set(ind, 1.0)
            nodeToNCS.set(ind, isNCS.elementSum())
        } else {
            isNCS.set(ind, 0.0)
        }
    }

    val non_cs_classes_list: MutableList<Int?> = ArrayList<Int?>()
    for (i in 0..<non_cs_classes.length()) {
        non_cs_classes_list.add(non_cs_classes.get(i).toInt())
    }

    val rtncs = dtmc_stochcomp(sn.rtnodes, non_cs_classes_list)
    var Inc = I
    for (i in sn.nodetype.indices) {
        if (sn.nodetype.get(i) == NodeType.ClassSwitch) {
            Inc = Inc - 1
        }
    }

    for (ist in 0..<DEP.size) {
        for (r in 0..<DEP.get(ist)!!.size) {
            if (DEP.get(ist)!!.get(r)!!.isEmpty || DEP.get(ist)!!.get(r)!!.get(0).hasNaN()) {
                DEP.get(ist)!!.get(r)!!.set(0, Matrix(1, 1, 0))
                DEP.get(ist)!!.get(r)!!.set(1, Matrix(1, 1, 0))
                DEP.get(ist)!!.get(r)!!.set(2, Matrix(1, 1, 0))
            } else {
                DEP.get(ist)!!.get(r)!!.set(2, DEP.get(ist)!!.get(r)!!.get(1))
            }
        }
    }


    val DEP_new: MutableMap<Int?, MutableMap<Int?, MatrixCell?>?> = HashMap<Int?, MutableMap<Int?, MatrixCell?>?>()
    val LINKS: MutableMap<Int?, MutableMap<Int?, MatrixCell?>?> = HashMap<Int?, MutableMap<Int?, MatrixCell?>?>()
    val ARV: MutableMap<Int?, MatrixCell?> = HashMap<Int?, MatrixCell?>()
    for (ind in 0..<I) {
        if (isNCS.get(ind) == 1.0) {
            val inc = nodeToNCS.get(ind).toInt()
            if (sn.nodetype.get(ind) == NodeType.Source || sn.nodetype.get(ind) == NodeType.Delay || sn.nodetype.get(ind) == NodeType.Queue) {
                val ist = sn.nodeToStation.get(ind).toInt()
                DEP_new.put(inc, HashMap<Int?, MatrixCell?>())
                if (R > 1) {
                    // Collect all class departure processes to superpose
                    var superposedMMAP: MatrixCell? = DEP.get(ist)!!.get(0)!!
                    for (r in 1..<R) {
                        superposedMMAP = mmap_super(superposedMMAP!!, DEP.get(ist)!!.get(r)!!)
                    }
                    DEP_new.get(inc)!!.put(0, superposedMMAP)
                } else {
                    DEP_new.get(inc)!!.put(0, DEP.get(ist)!!.get(0))
                }
                val Psplit = Matrix(R, Inc * R, R * Inc * R)
                for (r in 0..<R) {
                    for (jnd in 0..<I) {
                        if (isNCS.get(jnd) == 1.0) {
                            val jnc = nodeToNCS.get(jnd).toInt()
                            for (s in 0..<R) {
                                Psplit.set(r, (jnc - 1) * R + s, rtncs.get((inc - 1) * R + r, (jnc - 1) * R + s))
                            }
                        }
                    }
                }
                @Suppress("UNCHECKED_CAST")
                val Fsplit: MutableMap<Int?, MatrixCell?> =
                    npfqn_traffic_split_cs(DEP_new.get(inc)!!.get(0)!!, Psplit) as MutableMap<Int?, MatrixCell?>
                for (jnc in 0..<Inc) {
                    LINKS.put(inc, HashMap<Int?, MatrixCell?>())
                    LINKS.get(inc)!!.put(jnc, Fsplit.get(jnc))
                    LINKS.get(inc)!!.put(jnc, mmap_normalize(LINKS.get(inc)!!.get(jnc)!!))
                }
            }
        }
    }

    for (ind in 0..<I) {
        val FLOWS: MutableMap<Int?, MatrixCell?> = HashMap<Int?, MatrixCell?>()
        if (isNCS.get(ind) == 1.0 && sn.nodetype.get(ind) != NodeType.Source) {
            val inc = nodeToNCS.get(ind).toInt()
            for (jnd in 0..<Inc) {
                if (!LINKS.get(jnd)!!.get(inc)!!.isEmpty && mmap_lambda(LINKS.get(jnd)!!
                        .get(inc)!!).elementSum() > GlobalConstants.FineTol) {
                    FLOWS.put(FLOWS.size, LINKS.get(jnd)!!.get(inc))
                }
            }
            if (FLOWS.size > 1) {
                ARV.put(ind, npfqn_traffic_merge(FLOWS, config.merge, config.compress))
            } else if (FLOWS.size == 1) {
                ARV.put(ind, FLOWS.get(0))
            } else {
                // Use the last valid jnd from the previous loop
                var lastJnd = 0
                for (j in 0..<Inc) {
                    if (LINKS.get(j) != null && LINKS.get(j)!!.get(inc) != null) {
                        lastJnd = j
                    }
                }
                ARV.put(ind, LINKS.get(lastJnd)!!.get(0))
            }
        } else {
            ARV.put(ind, MatrixCell())
        }
    }

    return ARV
}