/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.mva.handlers

import jline.api.sn.snIsOpenModel
import jline.lang.NetworkStruct
import jline.GlobalConstants
import jline.GlobalConstants.Inf
import jline.lang.constant.NodeType
import jline.lang.constant.SchedStrategy
import jline.solvers.SolverOptions
import jline.solvers.mva.MVAResult
import jline.util.Utils
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath
import java.lang.Double
import kotlin.Int
import kotlin.RuntimeException
import kotlin.math.pow
import kotlin.math.sqrt

fun solver_qna(sn: NetworkStruct, options: SolverOptions): MVAResult {
    val config = options.config
    config.space_max = 1

    val K = sn.nclasses
    val rt = sn.rt.copy()
    val S = sn.rates.elementPow(-1.0)
    val scv = sn.scv.copy()
    scv.removeNaN()

    val I = sn.nnodes
    val M = sn.nstations
    val C = sn.nchains
    val V = Matrix.cellsum(sn.visits)
    val Q = Matrix(M, K, M * K)
    var QN_1 = Q.copy()
    for (i in 0..<QN_1.numRows) {
        for (j in 0..<QN_1.numCols) {
            QN_1.set(i, j, Inf)
        }
    }

    val U = Matrix(M, K, M * K)
    val R = Matrix(M, K, M * K)
    val T = Matrix(M, K, M * K)
    val X = Matrix(1, K, K)

    val lambda = Matrix(1, C, C)

    var it = 0

    for (i in 0..<sn.njobs.length()) {
        if (Double.isFinite(sn.njobs.get(i))) {
            throw RuntimeException("QNA does not support closed classes.")
        }
    }

    val a1 = Matrix(M, K, M * K)
    val a2 = Matrix(M, K, M * K)
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
    val lambdas_inchain: MutableMap<Int?, Matrix?> = HashMap<Int?, Matrix?>()
    val scvs_inchain: MutableMap<Int?, Matrix?> = HashMap<Int?, Matrix?>()
    val d2c = Matrix(1, C, C)
    var last_source_idx = 0
    if (snIsOpenModel(sn)) {
        for (c in 0..<C) {
            val inchain: Matrix = sn.inchain.get(c)!!
            val sourceIdx = sn.refstat.get(inchain.get(0).toInt()).toInt()
            last_source_idx = sourceIdx
            lambdas_inchain.put(c, Matrix(1, inchain.length(), inchain.length()))
            for (i in 0..<inchain.length()) {
                lambdas_inchain.get(c)!!.set(0, i, sn.rates.get(sourceIdx, inchain.get(i).toInt()))
            }
            scvs_inchain.put(c, Matrix(1, inchain.length(), inchain.length()))
            for (i in 0..<inchain.length()) {
                scvs_inchain.get(c)!!.set(0, i, scv.get(sourceIdx, inchain.get(i).toInt()))
            }
            val lambdas_inchain_C = lambdas_inchain.get(c)!!.copy()
            lambdas_inchain_C.removeInfinite()
            lambda.set(c, lambdas_inchain_C.elementSum())
            d2c.set(c, qna_superpos(lambdas_inchain.get(c)!!, scvs_inchain.get(c)!!))
            var openChain = false
            for (i in 0..<inchain.length()) {
                if (Utils.isInf(sn.njobs.get(inchain.get(i).toInt()))) {
                    openChain = true
                }
            }
            if (openChain) {
                for (i in 0..<inchain.length()) {
                    T.set(sourceIdx, inchain.get(i).toInt(), lambdas_inchain.get(c)!!.get(i))
                }
            }
        }
        d2.set(last_source_idx,
            Matrix.extractRows(d2c, last_source_idx, last_source_idx + 1, null).mult(lambda.transpose())
                .get(0) / lambda.elementSum())
    }
    var Q_diff = Q.add(-1.0, QN_1)
    Q_diff.absEq()
    while (Q_diff.elementMax() > options.iter_tol && it <= options.iter_max) {
        it = it + 1
        QN_1 = Q.copy()

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
                            a2.get(i, r) + 1 / lambda_i * f2.get(j * K + s, i * K + r) * T.get(j, s) * rt.get(j * K + s,
                                i * K + r))
                    }
                }
            }
        }

        for (ind in 0..<I) {
            if (sn.isstation.get(ind) == 1.0) {
                val ist = sn.nodeToStation.get(ind).toInt()
                if (sn.nodetype.get(ind) != NodeType.Join) {
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
                                val alpha_mi: kotlin.Double
                                if (rho_ist > 0.7) {
                                    alpha_mi = (rho_ist.pow(mi.toDouble()) + rho_ist) / 2
                                } else {
                                    alpha_mi = FastMath.pow(rho_ist, (mi + 1) / 2.0)
                                }
                                val mubar = lambda_ist / rho_ist
                                c2 = -1.0
                                for (r in 0..<K) {
                                    if (mu_ist.get(r) > 0) {
                                        c2 = c2 + a1.get(ist,
                                            r) / lambda_ist * (mubar / mi / mu_ist.get(r)).pow(2.0) * (scv.get(ist,
                                            r) + 1)
                                    }
                                }
                                val Wiq = (alpha_mi / mubar) * 1 / (1 - rho_ist) * (a2.sumRows(ist) + c2) / 2
                                Q.set(ist, k, a1.get(ist, k) / mu_ist.get(k) + a1.get(ist, k) * Wiq)
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
                            R.set(ist, k, Q.get(ist, k) / T.get(ist, k))
                        }
                    }
                }
            } else {
                if (sn.nodetype.get(ind) == NodeType.Fork) {
                    throw RuntimeException("Fork nodes not supported yet by QNA solver.")
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
        Q_diff = Q.add(-1.0, QN_1)
        Q_diff.absEq()
    }
    val result = MVAResult()
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
    result.logNormConstAggr = 0.0
    return result
}

fun qna_superpos(lambda: Matrix, a2: Matrix): kotlin.Double {
    val lambda_finite_idx: MutableList<Int?> = ArrayList<Int?>()
    for (i in 0..<lambda.length()) {
        if (Double.isFinite(lambda.get(i))) {
            lambda_finite_idx.add(i)
        }
    }
    val a2_new = Matrix(1, lambda_finite_idx.size, lambda_finite_idx.size)
    for (i in lambda_finite_idx.indices) {
        a2_new.set(i, a2.get(lambda_finite_idx.get(i)!!))
    }
    val lambda_new = Matrix(1, lambda_finite_idx.size, lambda_finite_idx.size)
    for (i in lambda_finite_idx.indices) {
        lambda_new.set(i, lambda.get(lambda_finite_idx.get(i)!!))
    }
    return a2_new.mult(lambda_new.transpose()).get(0) / lambda_new.elementSum()
}


