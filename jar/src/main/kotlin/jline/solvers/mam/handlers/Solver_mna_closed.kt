package jline.solvers.mam.handlers

import jline.api.mam.*
import jline.lang.NetworkStruct
import jline.GlobalConstants
import jline.lang.constant.NodeType
import jline.lang.constant.SchedStrategy
import jline.lang.processes.APH
import jline.lib.butools.MMAPPH1FCFS
import jline.api.mam.mmap_super_safe
import jline.solvers.SolverOptions
import jline.solvers.mam.MAMResult
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import org.apache.commons.math3.util.FastMath
import kotlin.math.abs
import kotlin.math.max
import kotlin.math.min
import kotlin.math.pow

fun solver_mna_closed(sn: NetworkStruct, options: SolverOptions): MAMResult {
    val config = options.config
    config.space_max = 16

    val K = sn.nclasses
    val rt = sn.rt.copy()
    val S = sn.rates.elementPow(-1.0)
    val scv = sn.scv.copy()
    scv.removeNaN()
    val PH = sn.proc
    val I = sn.nnodes
    val M = sn.nstations
    val C = sn.nchains
    val N = sn.njobs.copy()
    val V = Matrix.cellsum(sn.visits)

    val pie: MutableMap<Int?, MatrixCell?> = HashMap<Int?, MatrixCell?>()
    val D0: MutableMap<Int?, MatrixCell?> = HashMap<Int?, MatrixCell?>()
    var Q = Matrix(M, K, M * K)
    var U = Matrix(M, K, M * K)
    var R = Matrix(M, K, M * K)
    var T = Matrix(M, K, M * K)
    var X = Matrix(1, K, K)


    for (ist in 0..<M) {
        val sched = sn.sched.get(sn.stations.get(ist))
        if (sched == SchedStrategy.FCFS || sched == SchedStrategy.INF || sched == SchedStrategy.PS) {
            pie.put(ist, MatrixCell())
            D0.put(ist, MatrixCell())
            for (k in 0..<K) {
                val phProc = PH.get(sn.stations.get(ist))!!.get(sn.jobclasses.get(k))!!
                pie.get(ist)!!.set(k, map_pie(phProc.get(0), phProc.get(1)))
                val d0Val = phProc.get(0)
                if (d0Val.hasNaN()) {
                    D0.get(ist)!!.set(k, Matrix.singleton(-GlobalConstants.Immediate))
                    pie.get(ist)!!.set(k, Matrix.singleton(1.0))
                    PH.get(sn.stations.get(ist))!!.put(sn.jobclasses.get(k),
                        map_exponential(GlobalConstants.Immediate))
                } else {
                    D0.get(ist)!!.set(k, d0Val)
                }
            }
        }
    }

    var lambda = Matrix(1, C, C)
    val QNc = sn.njobs.copy()
    var it_out = 0
    val lambda_lb = Matrix(1, K, K)
    val lambda_ub = Matrix(1, K, K)
    for (k in 0..<K) {
        for (i in 0..<I) {
            if (sn.nservers.get(i) != Double.Companion.POSITIVE_INFINITY) {
                if (lambda_ub.get(k) == 0.0) {
                    lambda_ub.set(k, sn.rates.get(i, k))
                } else {
                    lambda_ub.set(k, min(lambda_ub.get(k), sn.rates.get(i, k)))
                }
            }
        }
    }


    HashMap<Int?, Matrix?>()
    HashMap<Int?, Matrix?>()
    Matrix(1, C, C)


    var QN = Matrix(1, K, K)
    var QN_diff = QN.add(-1.0, QNc)
    QN_diff.absEq()

    while (QN_diff.elementMax() > options.iter_tol && it_out < options.iter_max) {
        it_out++
        if (it_out == 1) {
            lambda = lambda_ub.copy()
        } else {
            for (k in 0..<K) {
                if (QN.get(k) < QNc.get(k)) {
                    lambda_lb.set(k, lambda.get(k))
                } else {
                    lambda_ub.set(k, lambda.get(k))
                }
                lambda.set(k, (lambda_ub.get(k) + lambda_lb.get(k)) / 2)
            }
        }
        var it = 0
        Q = Matrix(M, K, M * K)
        U = Matrix(M, K, M * K)
        R = Matrix(M, K, M * K)
        T = Matrix(M, K, M * K)
        X = Matrix(1, K, K)
        val a1 = Matrix(M, K, M * K)
        val a2 = Matrix(M, K, M * K)
        var a1_1 = a1.elementIncrease(Double.Companion.POSITIVE_INFINITY)

        var a2_1 = a2.elementIncrease(Double.Companion.POSITIVE_INFINITY)
        var a1_diff = a1.add(-1.0, a1_1)
        a1_diff.absEq()
        var a2_diff = a2.add(-1.0, a2_1)
        a2_diff.absEq()
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
        while ((a1_diff.elementMax() > options.iter_tol || a2_diff.elementMax() > options.iter_tol) && it <= options.iter_max) {
            it = it + 1

            // Chain-based queue normalization
            for (c in 0..<C) {
                val inchain = sn.inchain.get(c)!!
                var chainSum = 0.0
                for (k in 0..<inchain.length()) {
                    chainSum += Q.sumCols().get(inchain.get(k).toInt())
                }
                if (chainSum > 0) {
                    for (k in 0..<inchain.length()) {
                        val classIdx = inchain.get(k).toInt()
                        for (m in 0..<M) {
                            Q.set(m, classIdx, sn.njobs.get(c) * Q.get(m, classIdx) / chainSum)
                        }
                    }
                }
            }
            
            // Handle self-looping classes
            for (k in 0..<K) {
                if (sn.isslc.get(k) == 1.0) {
                    for (m in 0..<M) {
                        Q.set(m, k, 0.0)
                    }
                    Q.set(sn.refstat.get(k).toInt(), k, sn.njobs.get(k))
                }
            }
            
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
                    if (sn.nodetype.get(ind) == NodeType.Join) {
                        // no-op
                    } else {
                        if (sn.sched.get(sn.stations.get(ist)) == SchedStrategy.INF) {
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
                        } else if (sn.sched.get(sn.stations.get(ist)) == SchedStrategy.PS) {
                            for (c in 0..<C) {
                                val inchain = sn.inchain.get(c)!!
                                for (k1 in 0..<inchain.length()) {
                                    val k = inchain.get(k1).toInt()
                                    T.set(ist, k, lambda.get(c) * V.get(ist, k))
                                    U.set(ist, k, S.get(ist, k) * T.get(ist, k))
                                }
                                val Uden = min(1 - GlobalConstants.FineTol, U.sumRows(ist))
                                for (k1 in 0..<inchain.length()) {
                                    val k = inchain.get(k1).toInt()
                                    val Nc = sn.njobs.get(c)
                                    Q.set(ist, k, (U.get(ist, k) - U.get(ist, k).pow(Nc + 1)) / (1 - Uden))
                                    R.set(ist, k, Q.get(ist, k) / T.get(ist, k))
                                }
                            }
                        } else if (sn.sched.get(sn.stations.get(ist)) == SchedStrategy.FCFS) {
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
                                            c2 =
                                                c2 + a1.get(ist, r) / lambda_ist * FastMath.pow(mubar / mi / mu_ist.get(
                                                    r), 2) * (scv.get(ist, r) + 1)
                                        }
                                    }
                                }
                                d2.set(ist,
                                    1 + FastMath.pow(rho_ist,
                                        2) * (c2 - 1) / FastMath.sqrt(mi.toDouble()) + (1 - FastMath.pow(rho_ist,
                                        2)) * (a2.sumRows(ist) - 1))
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
                } else {
                    // not a station
                    if (sn.nodetype.get(ind) == NodeType.Fork) {
                        throw IllegalArgumentException("Fork nodes not supported yet by MNA solver")
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
                if (sn.sched.get(sn.stations.get(ist)) == SchedStrategy.FCFS) {
                    val rho_ist_class = Matrix(1, K, K)
                    for (i in 0..<K) {
                        rho_ist_class.set(i, a1.get(ist, i) / (GlobalConstants.FineTol + sn.rates.get(ist, i)))
                    }
                    rho_ist_class.removeNaN()

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
                                val mmapMap: Map<Int?, MatrixCell> = mapOf(0 to arri_node!!, 1 to arri_class)
                                arri_node = mmap_super_safe(mmapMap, config.space_max, "default")
                            }
                        }
                        val finite_N = N.copy()
                        finite_N.removeInfinite()
                        val maxLevel = finite_N.elementMax() + 1
                        val D = mmap_shorten(arri_node!!)
                        var pdistr: MutableMap<Int?, Matrix?> = HashMap<Int?, Matrix?>()
                        val Qret: MutableMap<Int?, Matrix?> = HashMap<Int?, Matrix?>()
                        if (map_lambda(D.get(0), D.get(1)) < GlobalConstants.FineTol) {
                            for (k in 0..<K) {
                                val pdistrK = Matrix(1, 2, 2)
                                pdistrK.set(0, 1 - GlobalConstants.FineTol)
                                pdistrK.set(1, GlobalConstants.FineTol)
                                pdistr.put(k, pdistrK)
                                Qret.put(k, Matrix.singleton(GlobalConstants.FineTol / sn.rates.get(ist)))
                            }
                        } else {
                            val result = MMAPPH1FCFS(D,
                                pie.get(ist)!!.toMap(),
                                D0.get(ist)!!.toMap(),
                                null,
                                maxLevel.toInt(),
                                null,
                                null,
                                false,
                                false,
                                null,
                                null).get("ncDistr")
                            pdistr = result?.let { map -> 
                                map.mapKeys { it.key as Int? }.mapValues { it.value as Matrix? }.toMutableMap()
                            } ?: mutableMapOf()
                            for (k in 0..<K) {
                                pdistr.put(k,
                                    Matrix.extractRows(pdistr.get(k)!!.transpose(), 0, N.get(k).toInt() + 1, null))
                                pdistr.get(k)!!.absEq()
                                var sum = 0.0
                                for (i in 0..<pdistr.get(k)!!.length() - 1) {
                                    sum = sum + pdistr.get(k)!!.get(i)
                                }
                                pdistr.get(k)!!.set(pdistr.get(k)!!.length() - 1, abs(1 - sum))
                                pdistr.get(k)!!.scaleEq(1 / pdistr.get(k)!!.elementSum())
                                val a = Matrix(1, N.get(k).toInt() + 1, N.get(k).toInt() + 1)
                                for (i in 0..<a.length()) {
                                    a.set(i, i.toDouble())
                                }
                                val b = Matrix(1, N.get(k).toInt() + 1, N.get(k).toInt() + 1)
                                for (i in 0..<a.length()) {
                                    b.set(i, pdistr.get(k)!!.get(i))
                                }
                                Qret.put(k, Matrix.singleton(max(0.0, min(N.get(k), a.mult(b.transpose()).get(0)))))
                            }
                        }
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
        QN = Q.sumCols()
        QN_diff = QN.add(-1.0, QNc)
        QN_diff.absEq()
    }

    // Handle self-looping classes final pass
    for (k in 0..<K) {
        if (sn.isslc.get(k) == 1.0) {
            for (m in 0..<M) {
                Q.set(m, k, 0.0)
            }
            val ist = sn.refstat.get(k).toInt()
            Q.set(ist, k, sn.njobs.get(k))
            T.set(ist, k, sn.njobs.get(k) * sn.rates.get(ist, k))
            R.set(ist, k, Q.get(ist, k) / T.get(ist, k))
            U.set(ist, k, S.get(ist, k) * T.get(ist, k))
        }
    }
    
    // Final chain normalization
    for (c in 0..<C) {
        val inchain = sn.inchain.get(c)!!
        if (sn.njobs.get(c).isFinite()) {
            var chainSum = 0.0
            for (k in 0..<inchain.length()) {
                chainSum += Q.sumCols().get(inchain.get(k).toInt())
            }
            if (chainSum > 0) {
                for (k in 0..<inchain.length()) {
                    val classIdx = inchain.get(k).toInt()
                    for (m in 0..<M) {
                        Q.set(m, classIdx, sn.njobs.get(c) * Q.get(m, classIdx) / chainSum)
                    }
                }
            }
        }
    }
    
    // Final adjustment for INF stations
    for (ist in 0..<M) {
        if (sn.sched.get(sn.stations.get(ist)) == SchedStrategy.INF) {
            for (k in 0..<K) {
                U.set(ist, k, Q.get(ist, k))
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
    result.iter = it_out
    result.lG = 0.0

    return result
}