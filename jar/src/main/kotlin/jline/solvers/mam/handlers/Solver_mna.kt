package jline.solvers.mam.handlers

import jline.api.mam.*
import jline.lang.NetworkStruct
import jline.GlobalConstants
import jline.lang.constant.NodeType
import jline.lang.constant.SchedStrategy
import jline.lang.processes.APH
import jline.lib.butools.MMAPPH1FCFS
import jline.solvers.SolverOptions
import jline.solvers.SolverResult
import jline.solvers.mam.MAMResult
import jline.util.Utils
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import kotlin.math.pow
import kotlin.math.sqrt

fun solver_mna(sn: NetworkStruct, options: SolverOptions): SolverResult {
    val config = options.config
    config.space_max = 1

    val K = sn.nclasses
    val rt = sn.rt.copy()
    val S = sn.rates.elementPow(-1.0)
    val scv = sn.scv.copy()
    scv.removeNaN()
    val PH = sn.proc
    val I = sn.nnodes
    val M = sn.nstations
    val C = sn.nchains
    val V = Matrix.cellsum(sn.visits)
    val Q = Matrix(M, K, M * K)
    val pie: MutableMap<Int?, MatrixCell?> = HashMap<Int?, MatrixCell?>()
    val DO: MutableMap<Int?, MatrixCell?> = HashMap<Int?, MatrixCell?>()


    val U = Matrix(M, K, M * K)
    val R = Matrix(M, K, M * K)
    val T = Matrix(M, K, M * K)
    val X = Matrix(1, K, K)
    for (ist in 0..<M) {
        if (sn.sched.get(sn.getStations().get(ist)) == SchedStrategy.FCFS || sn.sched.get(sn.getStations().get(ist)) == SchedStrategy.HOL || sn.sched.get(
                sn.getStations().get(ist)) == SchedStrategy.PS) {
            pie.put(ist, MatrixCell())
            DO.put(ist, MatrixCell())
            for (k in 0..<K) {
                PH.get(sn.getStations().get(ist))!!.put(sn.jobclasses.get(k),
                    map_scale(PH.get(sn.getStations().get(ist))!!.get(sn.jobclasses.get(k))!!.get(0),
                        PH.get(sn.getStations().get(ist))!!.get(sn.jobclasses.get(k))!!.get(1),
                        S.get(ist, k) / sn.nservers.get(ist)))
                pie.get(ist)!!.set(k,
                    map_pie(PH.get(sn.getStations().get(ist))!!.get(sn.jobclasses.get(k))!!.get(0),
                        PH.get(sn.getStations().get(ist))!!.get(sn.jobclasses.get(k))!!.get(1)))
                DO.get(ist)!!.set(k, PH.get(sn.getStations().get(ist))!!.get(sn.jobclasses.get(k))!!.get(0))
            }
        }
    }

    val lambda = Matrix(1, C, C)

    var it = 0

    val a1 = Matrix(M, K, M * K)
    val a2 = Matrix(M, K, M * K)
    var a1_1 = a1.elementIncrease(Double.Companion.POSITIVE_INFINITY)

    var a2_1 = a2.elementIncrease(Double.Companion.POSITIVE_INFINITY)

    val d2 = Matrix(M, 1, M)
    val f2 = Matrix(M * K, M * K, (M * K).toDouble().pow(2.0).toInt())
    for (i in 0..<M) {
        for (j in 0..<M) {
            if (sn.nodetype.get(sn.stationToNode.get(j).toInt()) != NodeType.Source) {
                for (r in 0..<K) {
                    for (s in 0..<K) {
                        if (rt.get(i * K + r, j * K + s) > 0) {
                            f2.set(i * K + r, j * K + s, 1)
                        }
                    }
                }
            }
        }
    }
    val lambdas_inchain = MatrixCell()
    val scvs_inchain = MatrixCell()
    val d2c = Matrix(1, C, C)
    var last_source_idx = 0
    for (c in 0..<C) {
        val inchain: Matrix = sn.inchain.get(c)!!
        val sourceIdx = sn.refstat.get(inchain.get(0).toInt()).toInt()
        last_source_idx = sourceIdx
        lambdas_inchain.set(c, Matrix(1, inchain.length(), inchain.length()))
        for (i in 0..<inchain.length()) {
            lambdas_inchain.get(c).set(0, i, sn.rates.get(sourceIdx, inchain.get(i).toInt()))
        }
        scvs_inchain.set(c, Matrix(1, inchain.length(), inchain.length()))
        for (i in 0..<inchain.length()) {
            scvs_inchain.get(c).set(0, i, scv.get(sourceIdx, inchain.get(i).toInt()))
        }
        val lambdas_inchain_C = lambdas_inchain.get(c).copy()
        lambdas_inchain_C.removeInfinite()
        lambda.set(c, lambdas_inchain_C.elementSum())
        d2c.set(c, qna_superpos(lambdas_inchain.get(c), scvs_inchain.get(c)))
        var openChain = false
        for (i in 0..<inchain.length()) {
            if (Utils.isInf(sn.njobs.get(inchain.get(i).toInt()))) {
                openChain = true
            }
        }

        for (i in 0..<inchain.length()) {
            T.set(sourceIdx, inchain.get(i).toInt(), lambdas_inchain.get(c).get(i))
        }
    }
    d2.set(last_source_idx,
        Matrix.extractRows(d2c, last_source_idx, last_source_idx + 1, null).mult(lambda.transpose())
            .get(0) / lambda.elementSum())
    var a1_diff = a1.add(-1.0, a1_1)
    a1_diff.absEq()
    var a2_diff = a2.add(-1.0, a2_1)
    a2_diff.absEq()

    while ((a1_diff.elementMax() > options.iter_tol || a2_diff.elementMax() > options.iter_tol) && it <= options.iter_max) {
        it = it + 1
        a1_1 = a1.copy()
        a2_1 = a2.copy()

        if (it == 1) {
            for (c in 0..<C) {
                val inchain: Matrix = sn.inchain.get(c)!!
                for (m in 0..<M) {
                    for (i in 0..<inchain.length()) {
                        T.set(m, inchain.get(i).toInt(), V.get(m, inchain.get(i).toInt()) * lambda.get(c))
                    }
                }
            }
        }
        for (i in 0..<M) {
            for (j in 0..<K) {
                a1.set(i, j, 0)
                a2.set(i, j, 0)
            }
            val lambda_i = T.sumRows(i)
            for (j in 0..<M) {
                for (r in 0..<K) {
                    for (s in 0..<K) {
                        a1.set(i, r, a1.get(i, r) + T.get(j, s) * rt.get(j * K + s, i * K + r))
                        a2.set(i,
                            r,
                            a2.get(i, r) + (1 / lambda_i) * f2.get(j * K + s, i * K + r) * T.get(j,
                                s) * rt.get(j * K + s, i * K + r))
                    }
                }
            }
        }

        for (ind in 0..<I) {
            if (sn.isstation.get(ind) == 1.0) {
                val ist = sn.nodeToStation.get(ind).toInt()
                if (sn.nodetype.get(ind) == NodeType.Fork) {
                    // Fork nodes: fan-out operation with zero queueing metrics like Join nodes
                    for (k in 0..<K) {
                        Q.set(ist, k, 0.0)
                        U.set(ist, k, 0.0)
                        T.set(ist, k, a1.get(ist, k))
                    }
                    d2.set(ist, 1.0)
                } else if (sn.nodetype.get(ind) != NodeType.Join) {
                    if (sn.sched.get(sn.getStations().get(ist)) == SchedStrategy.INF) {
                        for (i in 0..<M) {
                            for (r in 0..<K) {
                                for (s in 0..<K) {
                                    d2.set(ist, s, a2.get(ist, s))
                                }
                            }
                        }
                        for (c in 0..<C) {
                            val inchain: Matrix = sn.inchain.get(c)!!
                            for (k1 in 0..<inchain.length()) {
                                val k = inchain.get(k1).toInt()
                                T.set(ist, k, a1.get(ist, k))
                                U.set(ist, k, S.get(ist, k) * T.get(ist, k))
                                Q.set(ist, k, T.get(ist, k) * S.get(ist, k) * V.get(ist, k))
                                R.set(ist, k, Q.get(ist, k) / T.get(ist, k))
                            }
                        }
                    } else if (sn.sched.get(sn.getStations().get(ist)) == SchedStrategy.FCFS) {
                        val mu_ist = Matrix(1, K, K)
                        for (i in 0..<K) {
                            mu_ist.set(i, sn.rates.get(ist, i))
                        }
                        mu_ist.removeNaN()
                        val rho_ist_class = Matrix(1, K, K)
                        for (i in 0..<K) {
                            rho_ist_class.set(i, a1.get(ist, i) / (GlobalConstants.FineTol + sn.rates.get(ist, i)))
                        }
                        rho_ist_class.removeNaN()
                        val lambda_ist = a1.sumRows(ist)
                        val mi = sn.nservers.get(ist).toInt()
                        val rho_ist = rho_ist_class.elementSum() / mi
                        var c2 = 0.0
                        if (rho_ist < 1 - options.tol) {
                            for (k in 0..<K) {
                                val mubar = lambda_ist / rho_ist
                                c2 = -1.0
                                for (r in 0..<K) {
                                    if (mu_ist.get(r) > 0) {
                                        c2 = c2 + a1.get(ist,
                                            r) / lambda_ist * (mubar / mi / mu_ist.get(r)).pow(2.0) * (scv.get(ist,
                                            r) + 1)
                                    }
                                }
                            }
                            d2.set(ist,
                                1 + rho_ist.pow(2.0) * (c2 - 1) / sqrt(mi.toDouble()) + (1 - rho_ist.pow(2.0)) * (a2.sumRows(
                                    ist) - 1))
                        } else {
                            for (k in 0..<K) {
                                Q.set(ist, k, sn.njobs.get(k))
                            }
                            d2.set(ist, 1.0)
                        }
                        for (k in 0..<K) {
                            T.set(ist, k, a1.get(ist, k))
                            U.set(ist, k, T.get(ist, k) * S.get(ist, k) / sn.nservers.get(ist))
                        }
                    }
                }
            }
        }

        for (i in 0..<M) {
            for (j in 0..<M) {
                if (sn.nodetype.get(sn.stationToNode.get(j).toInt()) != NodeType.Source) {
                    for (r in 0..<K) {
                        for (s in 0..<K) {
                            if (rt.get(i * K + r, j * K + s) > 0) {
                                f2.set(i * K + r, j * K + s, 1 + rt.get(i * K + r, j * K + s) * (d2.get(i) - 1))
                            }
                        }
                    }
                }
            }
        }
        a1_diff = a1.add(-1.0, a1_1)
        a1_diff.absEq()
        a2_diff = a2.add(-1.0, a2_1)
        a2_diff.absEq()
    }

    for (ind in 0..<I) {
        if (sn.isstation.get(ind) == 1.0) {
            val ist = sn.nodeToStation.get(ind).toInt()
            if (sn.sched.get(sn.getStations().get(ist)) == SchedStrategy.FCFS) {
                val mu_ist = Matrix(1, K, K)
                for (i in 0..<K) {
                    mu_ist.set(i, sn.rates.get(ist, i))
                }
                mu_ist.removeNaN()
                val rho_ist_class = Matrix(1, K, K)
                for (i in 0..<K) {
                    rho_ist_class.set(i, a1.get(ist, i) / (GlobalConstants.FineTol + sn.rates.get(ist, i)))
                }
                rho_ist_class.removeNaN()
                a1.sumRows(ist)
                val mi = sn.nservers.get(ist).toInt()
                val rho_ist = rho_ist_class.elementSum() / mi
                if (rho_ist < 1 - options.tol) {
                    var arri_class_anx = MatrixCell()
                    val arri_class = MatrixCell()
                    var arri_node: MatrixCell? = MatrixCell()
                    for (k in 0..<K) {
                        if (a1.get(ist, k) == 0.0) {
                            arri_class_anx = map_exponential(Double.Companion.POSITIVE_INFINITY)
                        } else {
                            arri_class_anx = APH.fitMeanAndSCV(1 / a1.get(ist, k), a2.get(ist, k)).getProcess()
                        }
                        arri_class.set(0, arri_class_anx.get(0))
                        arri_class.set(1, arri_class_anx.get(1))
                        arri_class.set(2, arri_class_anx.get(1))
                        if (k == 0) {
                            arri_node!!.set(0, arri_class.get(0))
                            arri_node.set(1, arri_class.get(1))
                            arri_node.set(2, arri_class.get(2))
                        } else {
                            arri_node = mmap_super(arri_node!!, arri_class)
                        }
                    }
                    var Qret: MutableMap<Int?, Matrix?> = HashMap<Int?, Matrix?>()
                    val result = MMAPPH1FCFS(mmap_shorten(arri_node!!),
                        pie.get(ist)!!.toMap(),
                        DO.get(ist)!!.toMap(),
                        1,
                        null,
                        null,
                        null,
                        false,
                        false,
                        null,
                        null).get("ncMoms")
                    Qret = result?.let { map -> 
                        map.mapKeys { it.key as Int? }.mapValues { it.value as Matrix? }.toMutableMap()
                    } ?: mutableMapOf()
                    for (i in 0..<Qret.size) {
                        Q.set(ist, i, Qret.get(i)!!.get(0))
                    }
                } else {
                    for (k in 0..<K) {
                        Q.set(ist, k, sn.njobs.get(k))
                    }
                }
                for (k in 0..<K) {
                    R.set(ist, k, Q.get(ist, k) / T.get(ist, k))
                }
            }
        }
    }

    val result = MAMResult()
    val CN = R.sumCols()
    Q.absEq()
    Q.removeNaN()
    U.removeNaN()
    R.removeNaN()
    CN.removeNaN()
    X.removeNaN()
    result.XN = X
    result.QN = Q
    result.CN = CN
    result.UN = U
    result.RN = R
    result.TN = T
    result.iter = it

    return result
}