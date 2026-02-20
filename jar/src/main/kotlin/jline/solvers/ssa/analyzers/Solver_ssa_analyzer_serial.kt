package jline.solvers.ssa.analyzers

import jline.api.mam.map_mean
import jline.lang.NetworkStruct
import jline.lang.constant.NodeType
import jline.lang.constant.SchedStrategy
import jline.lang.nodeparam.CacheNodeParam
import jline.lang.nodes.Cache
import jline.lang.nodes.StatefulNode
import jline.lang.state.EventCache
import jline.solvers.SolverOptions
import jline.solvers.ssa.SSAResult
import jline.solvers.ssa.SolverSSA
import jline.solvers.ssa.handlers.solver_ssa
import jline.util.matrix.Matrix

fun solver_ssa_analyzer_serial(sn: NetworkStruct,
                               hash: Boolean,
                               init_state: MutableMap<StatefulNode?, Matrix?>,
                               options: SolverOptions,
                               solverSSA: SolverSSA): SSAResult {
    val M = sn.nstations
    val K = sn.nclasses

    val S = sn.nservers
    val NK = sn.njobs.transpose()
    val schedid = sn.sched

    val PH = sn.proc
    var tranSysState: MutableMap<Int?, Matrix?>? = HashMap<Int?, Matrix?>()
    var tranSync: Matrix? = Matrix(0, 0)

    val XN = Matrix(1, K)
    XN.fill(Double.NaN)
    val UN = Matrix(M, K)
    UN.fill(Double.NaN)
    val QN = Matrix(M, K)
    QN.fill(Double.NaN)
    val RN = Matrix(M, K)
    RN.fill(Double.NaN)
    val TN = Matrix(M, K)
    TN.fill(Double.NaN)
    val CN = Matrix(1, K)
    CN.fill(Double.NaN)

    options.samples++
    solverSSA.eventCache = EventCache(false, options.config.eventcache)
    val result = solver_ssa(sn, solverSSA.eventCache, init_state, options, solverSSA)
    val probSysState = result.pi
    val StateSpaceAggr = result.SSq
    val arvRates = result.arvRates
    val depRates = result.depRates
    tranSysState = result.tranSysState
    tranSync = result.tranSync


    for (k in 0..<K) {
        val refsf = sn.stationToStateful.get(sn.refstat.get(k).toInt())
        val departure: Matrix = depRates.get(k)!!
        val dep_wset_refsf = Matrix(StateSpaceAggr.numRows, 1)
        for (i in 0..<StateSpaceAggr.numRows) {
            dep_wset_refsf.set(i, 0, departure.get(i, refsf.toInt()))
        }
        // toDouble call since 1xn mult nx1
        XN.set(k, probSysState.mult(dep_wset_refsf).toDouble())
    }

    for (i in 0..<M) {
        val isf = sn.stationToStateful.get(i).toInt()
        for (k in 0..<K) {
            val departure: Matrix = depRates.get(k)!!
            val dep_wset_isf = Matrix(StateSpaceAggr.numRows, 1)
            for (j in 0..<StateSpaceAggr.numRows) {
                dep_wset_isf.set(j, 0, departure.get(j, isf))
            }
            TN.set(i, k, probSysState.mult(dep_wset_isf).toDouble())

            val ssaggr_wset_isf =
                Matrix.extract(StateSpaceAggr, 0, StateSpaceAggr.numRows, (i) * K + k, (i) * K + k + 1)
            QN.set(i, k, probSysState.mult(ssaggr_wset_isf).toDouble())
        }
        when (schedid.get(sn.stations.get(i))) {
            SchedStrategy.INF -> {
                var k = 0
                while (k < K) {
                    UN.set(i, k, QN.get(i, k))
                    k++
                }
            }

            SchedStrategy.PS, SchedStrategy.DPS, SchedStrategy.GPS,
            SchedStrategy.PSPRIO, SchedStrategy.DPSPRIO, SchedStrategy.GPSPRIO -> {
                if ((sn.lldscaling == null || sn.lldscaling.isEmpty) && (sn.cdscaling == null || sn.cdscaling.isEmpty())) {
                    var k = 0
                    while (k < K) {
                        if (!PH.get(sn.stations.get(i))!!.get(sn.jobclasses.get(k))!!.isEmpty) {
                            val arrival: Matrix = arvRates.get(k)!!
                            val arv_wset_isf = Matrix(StateSpaceAggr.numRows, 1)
                            var c = 0
                            while (c < StateSpaceAggr.numRows) {
                                arv_wset_isf.set(c, 0, arrival.get(c, isf))
                                c++
                            }
                            // For PS, DPS, GPS: UN = probSysState * arvRates / sn.rates / S
                            UN.set(i, k, probSysState.mult(arv_wset_isf).toDouble() / sn.rates.get(i, k) / S.get(i))
                        }
                        k++
                    }
                } else {
                    // lld/cd cases
                    var col = 0
                    while (col < K) {
                        UN.set(i, col, Double.NaN)
                        col++
                    }
                }
            }

            else -> if ((sn.lldscaling == null || sn.lldscaling.isEmpty) && (sn.cdscaling == null || sn.cdscaling.isEmpty())) {
                var k = 0
                while (k < K) {
                    if (!PH.get(sn.stations.get(i))!!.get(sn.jobclasses.get(k))!!.isEmpty) {
                        val arrival: Matrix = arvRates.get(k)!!
                        val arv_wset_isf = Matrix(StateSpaceAggr.numRows, 1)
                        var c = 0
                        while (c < StateSpaceAggr.numRows) {
                            arv_wset_isf.set(c, 0, arrival.get(c, isf))
                            c++
                        }
                        val map_mean = map_mean(PH.get(sn.stations.get(i))!!.get(sn.jobclasses.get(k))!!.get(0),
                            PH.get(sn.stations.get(i))!!.get(sn.jobclasses.get(k))!!.get(1)) / S.get(i)
                        UN.set(i, k, probSysState.mult(arv_wset_isf).toDouble() * map_mean)
                    }
                    k++
                }
            } else {
                // lld/cd cases
                // int ind = (int) sn.stationToNode.get(i);
                var col = 0
                while (col < K) {
                    UN.set(i, col, Double.NaN)
                    col++
                }
            }

        }
    }

    for (k in 0..<K) {
        for (i in 0..<M) {
            if (TN.get(i, k) > 0) {
                RN.set(i, k, QN.get(i, k) / TN.get(i, k))
            } else {
                RN.set(i, k, 0)
            }
        }
        CN.set(k, NK.get(k) / XN.get(k))
    }

    // update routing probabilities in nodes with state-dependent routing
    val TNcache = Matrix(sn.nstateful, K)
    for (k in 0..<K) {
        for (isf in 0..<sn.nstateful) {
            if (sn.nodetype.get(isf) == NodeType.Cache) {
                val TNcacheValue = probSysState.mult(depRates.get(k)!!.getColumn(isf)).get(0)
                TNcache.set(isf, k, TNcacheValue)
            }
        }
    }

    for (k in 0..<K) {
        for (isf in 0..<sn.nstateful) {
            val statefulNode = solverSSA.model.getStatefulNodes().get(isf)
            if (statefulNode is Cache) {
                if ((sn.nodeparam.get(statefulNode) as CacheNodeParam).hitclass.numCols > k) {
                    val h = (sn.nodeparam.get(statefulNode) as CacheNodeParam).hitclass.get(k).toInt()
                    val m = (sn.nodeparam.get(statefulNode) as CacheNodeParam).missclass.get(k).toInt()
                    if (h == -1 || m == -1) {
                        if ((sn.nodeparam.get(statefulNode) as CacheNodeParam).actualhitprob == null) {
                            val actualhitprobMatrix = Matrix(1, K)
                            actualhitprobMatrix.set(k, Double.NaN)
                            (sn.nodeparam.get(statefulNode) as CacheNodeParam).actualhitprob =
                                actualhitprobMatrix
                        } else {
                            (sn.nodeparam.get(statefulNode) as CacheNodeParam).actualhitprob.set(k, Double.NaN)
                        }
                        // updates cache actual hit and miss data
                        if ((sn.nodeparam.get(statefulNode) as CacheNodeParam).actualmissprob == null) {
                            val actualhitprobMatrix = Matrix(1, K)
                            actualhitprobMatrix.set(k, Double.NaN)
                            (sn.nodeparam.get(statefulNode) as CacheNodeParam).actualmissprob =
                                actualhitprobMatrix
                        } else {
                            (sn.nodeparam.get(statefulNode) as CacheNodeParam).actualmissprob.set(k, Double.NaN)
                        }
                    } else {
                        val actualhitprobvalue = TNcache.get(isf, h) / (TNcache.get(isf, h) + TNcache.get(isf, m))
                        if ((sn.nodeparam.get(statefulNode) as CacheNodeParam).actualhitprob == null) {
                            val actualhitprobMatrix = Matrix(1, K)
                            actualhitprobMatrix.set(k, actualhitprobvalue)
                            (sn.nodeparam.get(statefulNode) as CacheNodeParam).actualhitprob =
                                actualhitprobMatrix
                        } else {
                            (sn.nodeparam.get(statefulNode) as CacheNodeParam).actualhitprob.set(k,
                                actualhitprobvalue)
                        }

                        val actualmissprobvalue = TNcache.get(isf, m) / (TNcache.get(isf, h) + TNcache.get(isf, m))
                        if ((sn.nodeparam.get(statefulNode) as CacheNodeParam).actualmissprob == null) {
                            val actualhitprobMatrix = Matrix(1, K)
                            actualhitprobMatrix.set(k, actualmissprobvalue)
                            (sn.nodeparam.get(statefulNode) as CacheNodeParam).actualmissprob =
                                actualhitprobMatrix
                        } else {
                            (sn.nodeparam.get(statefulNode) as CacheNodeParam).actualmissprob.set(k,
                                actualmissprobvalue)
                        }
                    }
                }
            }
        }
    }

    // matrices QN, CN, RN, UN, XN, TN, where they are Double.isNan set to 0
    QN.apply(Double.NaN, 0.0, "equal")
    CN.apply(Double.NaN, 0.0, "equal")
    RN.apply(Double.NaN, 0.0, "equal")
    UN.apply(Double.NaN, 0.0, "equal")
    XN.apply(Double.NaN, 0.0, "equal")
    TN.apply(Double.NaN, 0.0, "equal")

    return SSAResult(QN, UN, RN, TN, CN, XN, tranSysState, tranSync, sn)
}
