/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.mva.handlers

import jline.GlobalConstants
import jline.GlobalConstants.Inf
import jline.api.pfqn.ld.pfqn_cdfun
import jline.api.pfqn.ld.pfqn_lldfun
import jline.api.pfqn.ld.ljd_linearize
import jline.api.sn.snHasJointDependence
import jline.api.sn.snDeaggregateChainResults
import jline.api.sn.snGetDemandsChain
import jline.io.Ret.snDeaggregateChainResults
import jline.lang.NetworkStruct
import jline.lang.constant.SchedStrategy
import jline.lang.constant.SolverType
import jline.solvers.SolverOptions
import jline.solvers.mva.MVAResult
import jline.lang.nodes.Station
import jline.util.SerializableFunction
import jline.util.Pair
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath
import java.lang.Double
import java.util.*
import kotlin.Int
import kotlin.RuntimeException
import kotlin.math.abs
import kotlin.math.max
import kotlin.math.min


/**
 * Handler for the solver_amvald function
 *//*
 * Note: there are some discrepancies between this and solver_amvald in LINE. There is an error around 0.0001 for
 * ld_multiserver_fcfs and around 0.001 for ld_multiserver_ps_twoclasses. TODO: fix this discrepancy
 * */
fun solver_amvald(sn: NetworkStruct, options: SolverOptions?): MVAResult {
    var options = options ?: SolverOptions(SolverType.MVA)
    val startTime = System.nanoTime()

    val res = snGetDemandsChain(sn)
    val Lchain = res.Dchain
    val STchain = res.STchain
    val Vchain = res.Vchain
    val alpha = res.alpha
    val Nchain = res.Nchain
    val SCVchain = res.SCVchain
    val refstatchain = res.refstatchain

    val M = sn.nstations
    val K = sn.nchains
    var Nt = 0.0
    for (col in 0..<Nchain.numCols) {
        if (Double.isFinite(Nchain.get(0, col))) Nt += Nchain.get(0, col)
    }
    val tol = options.iter_tol
    val nservers = sn.nservers
    val schedparam = sn.schedparam
    val lldscaling = sn.lldscaling
    val cdscaling = sn.cdscaling
    val ljdscaling = sn.ljdscaling
    val ljdcutoffs = sn.ljdcutoffs
    val stations = sn.stations
    val jobClasses = sn.jobclasses

    // Convert ljcdscaling from Map<Station, Map<JobClass, Matrix>> to indexed List<List<Matrix>?>
    // where first index is station index, second is class index
    val ljcdscalingIndexed: List<List<Matrix>?>? = if (sn.ljcdscaling != null) {
        (0 until M).map { i ->
            val station = stations[i]
            val stationMap = sn.ljcdscaling!![station]
            if (stationMap != null) {
                jobClasses.map { jc -> stationMap[jc] ?: Matrix(0, 0) }
            } else {
                null
            }
        }
    } else {
        null
    }

    // Convert ljcdcutoffs from Map<Station, Matrix> to indexed List<Matrix?>
    val ljcdcutoffsIndexed: List<Matrix?>? = if (sn.ljcdcutoffs != null) {
        (0 until M).map { i -> sn.ljcdcutoffs!![stations[i]] }
    } else {
        null
    }

    val sched: List<SchedStrategy> = (0 until M).map { i ->
        sn.sched[sn.stations[i]]!!
    }

    val Uchain = Matrix(M, K)
    val Tchain = Matrix(M, K)
    val Rchain = Matrix(M, K)
    val Cchain_s = Matrix(1, K)

    var Qchain = if (options.init_sol != null) options.init_sol.copy() else null
    if ((Qchain == null) || (Qchain.isEmpty)) {
        // balanced initialization to match MATLAB implementation
        Qchain = Matrix.ones(M, K)

        // Normalize by column sums and multiply by Nchain values
        for (col in 0..<K) {
            val colSum = Qchain.sumCols(col)
            for (row in 0..<M) {
                Qchain.set(row, col, (Qchain.get(row, col) / colSum) * Nchain.get(col))
            }
        }

        // Set infinite values to 0 (open classes)
        Qchain.apply(Inf, 0.0, "equal")

        // For open classes, set reference station queue length to 0
        for (r in 0..<Nchain.numCols) {
            if (Double.isInfinite(Nchain.get(r))) {
                Qchain.set(refstatchain.get(r,).toInt(), r, 0.0)
            }
        }
    }

    val nnzclasses: MutableList<Int> = ArrayList<Int>()
    val Xchain = Matrix(1, STchain.numCols)
    for (r in 0..<Nchain.numCols) {
        if (Double.isInfinite(Nchain.get(0, r))) Xchain.set(0, r, 1.0 / STchain.get(refstatchain.get(r, 0).toInt(), r))
        else Xchain.set(0, r, 1.0 / STchain.sumCols(r))

        if (Nchain.get(0, r) > 0) nnzclasses.add(r)
    }

    for (k in 0..<M) {
        for (r in nnzclasses) {
            if (Double.isInfinite(nservers.get(k, 0))) Uchain.set(k,
                r,
                Vchain.get(k, r) * STchain.get(k, r) * Xchain.get(0, r))
            else Uchain.set(k, r, (Vchain.get(k, r) * STchain.get(k, r) * Xchain.get(0, r)) / nservers.get(k, 0))
        }
    }

    if (options.config.np_priority == null) options.config.np_priority = "default"
    if (options.config.multiserver == null) options.config.multiserver = "default"
    if (options.config.highvar == null) options.config.highvar = "default"

    if (options.method == "default") {
        if (Nt <= 2) options.method = "bs"
        else options.method = "lin"
    }

    //Use list JLineMatrix to represent 3-D Matrix
    val gamma: MutableList<Matrix> = ArrayList<Matrix>()
    val tau = Matrix(K, K)
    when (options.method) {
        "default", "amva_lin", "lin", "amva_qdlin", "qdlin" -> {
            var i = 0
            while (i < K) {
                gamma.add(Matrix(K, M))
                i++
            }
        }

        else -> gamma.add(Matrix(K, M))
    }

    /* Main Loop */
    val omicron = 0.5 // under-relaxation parameter
    var outer_iter = 0.0
    var totiter = 0
    var QchainOuter_1 = Qchain.copy()
    var XchainOuter_1 = Xchain.copy()
    var UchainOuter_1 = Uchain.copy()
    var STeff: Matrix? = null
    while ((outer_iter < 2 || Qchain.sub(QchainOuter_1)
            .elementMaxAbs() > tol) && (outer_iter < FastMath.sqrt(options.iter_max.toDouble()))) {
        outer_iter++

        QchainOuter_1 = Qchain.copy()
        XchainOuter_1 = Xchain.copy()
        UchainOuter_1 = Uchain.copy()

        if (Double.isFinite(Nt) && Nt > 0) {
            // options.method.equals("aql") ||
            //                        options.method.equals("qdaql")
            if (options.method == "default" || options.method == "lin" || options.method == "qdlin") {/* Iteration at population N-1_s */
                for (s in 0..<K) {
                    if (Double.isFinite(Nchain.get(0, s))) {
                        var iter_s = 0.0
                        val Nchain_s = Matrix.oner(Nchain, ArrayList<Int?>(mutableListOf<Int?>(s)))
                        val Qchain_s = Qchain.copy()
                        Qchain_s.scaleEq((Nt - 1) / Nt, Qchain_s)
                        val Xchain_s = Xchain.copy()
                        Xchain_s.scaleEq((Nt - 1) / Nt, Xchain_s)
                        val Uchain_s = Uchain.copy()
                        Uchain_s.scaleEq((Nt - 1) / Nt, Uchain_s)
                        var Qchain_s_1 = Qchain_s.copy()
                        var Xchain_s_1 = Xchain_s.copy()
                        var Uchain_s_1 = Uchain_s.copy()

                        var count_inf = 0
                        var Nt_s = 0.0
                        val deltaclass = Matrix(Nchain_s.numRows, Nchain_s.numCols)
                        val ocl: MutableList<Int?> = ArrayList<Int?>()
                        val ccl: MutableList<Int> = ArrayList<Int>()
                        val nnzclasses_s: MutableList<Int> = ArrayList<Int>()
                        for (col in 0..<Nchain_s.numCols) {
                            val nchainValue = Nchain_s.get(0, col)
                            if (Double.isInfinite(nchainValue)) {
                                count_inf++
                                deltaclass.set(0, col, 1.0)
                                ocl.add(col)
                            } else {
                                Nt_s += nchainValue
                                deltaclass.set(0, col, (nchainValue - 1) / nchainValue)
                                if (nchainValue > 0) ccl.add(col)
                            }

                            if (nchainValue > 0) nnzclasses_s.add(col)
                        }
                        val delta = if (count_inf == Nchain_s.numCols) 1.0 else (Nt_s - 1) / Nt_s

                        //Use List<Integer> instead of JLineMatrix to store the index to avoid type transformation
                        val nnzclasses_eprio: MutableMap<Int?, MutableList<Int?>?> =
                            HashMap<Int?, MutableList<Int?>?>(nnzclasses_s.size)
                        val nnzclasses_hprio: MutableMap<Int?, MutableList<Int>?> =
                            HashMap<Int?, MutableList<Int>?>(nnzclasses_s.size)
                        val nnzclasses_ehprio: MutableMap<Int?, MutableList<Int>?> =
                            HashMap<Int?, MutableList<Int>?>(nnzclasses_s.size)
                        for (r in nnzclasses_s) {
                            val prio = sn.classprio.get(0, r)
                            val eprio_list: MutableList<Int?> = ArrayList<Int?>()
                            val hprio_list: MutableList<Int?> = ArrayList<Int?>()
                            for (i in 0..<sn.classprio.numCols) {
                                if (Double.compare(prio, sn.classprio.get(0, i)) == 0) eprio_list.add(i)
                                else if (Double.compare(prio, sn.classprio.get(0, i)) > 0) hprio_list.add(i) // higher prio (lower value = higher priority)
                            }
                            val eprio_common: MutableList<Int?> = ArrayList<Int?>(nnzclasses_s)
                            val hprio_common: MutableList<Int?> = ArrayList<Int?>(nnzclasses_s)
                            eprio_common.retainAll(eprio_list)
                            hprio_common.retainAll(hprio_list)
                            nnzclasses_eprio.put(r, eprio_common)
                            @Suppress("UNCHECKED_CAST")
                            nnzclasses_hprio.put(r, hprio_common as MutableList<Int>?)

                            val ehprio_common: MutableSet<Int> = LinkedHashSet<Int>(eprio_common)
                            ehprio_common.addAll(hprio_common)
                            @Suppress("UNCHECKED_CAST")
                            nnzclasses_ehprio.put(r, ArrayList<Int?>(ehprio_common) as MutableList<Int>?)
                        }
                        while ((iter_s < 2 || Qchain_s.sub(Qchain_s_1)
                                .elementMaxAbs() > tol) && (iter_s <= FastMath.sqrt(options.iter_max.toDouble()))) {
                            iter_s++

                            Qchain_s_1 = Qchain_s.copy()
                            Xchain_s_1 = Xchain_s.copy()
                            Uchain_s_1 = Uchain_s.copy()

                            val ret = solver_amvald_forward(gamma,
                                tau,
                                Qchain_s_1,
                                Xchain_s_1,
                                Uchain_s_1,
                                STchain,
                                Vchain,
                                Nchain_s,
                                SCVchain,
                                Nt_s,
                                delta,
                                deltaclass,
                                ocl,
                                ccl,
                                nnzclasses_s,
                                nnzclasses_eprio,
                                nnzclasses_hprio,
                                nnzclasses_ehprio,
                                M,
                                K,
                                nservers,
                                schedparam,
                                lldscaling,
                                cdscaling,
                                ljdscaling,
                                ljdcutoffs,
                                ljcdscalingIndexed,
                                ljcdcutoffsIndexed,
                                sched,
                                stations,
                                options)
                            totiter++
                            if (totiter > options.iter_max) {
                                break
                            }
                            val Wchain_s = ret.getLeft()
                            val STeff_s = ret.getRight()

                            for (r in nnzclasses) {
                                if (Wchain_s.sumCols(r) == 0.0) {
                                    Xchain_s.remove(0, r)
                                } else {
                                    if (Double.isInfinite(Nchain_s.get(0, r))) {
                                        var sumValue = 0.0
                                        for (row in 0..<Vchain.numRows) sumValue += Vchain.get(row, r) * Wchain_s.get(
                                            row,
                                            r)
                                        Cchain_s.set(0, r, sumValue)
                                    } else if (Nchain.get(0, r) == 0.0) {
                                        Xchain_s.remove(0, r)
                                        Cchain_s.remove(0, r)
                                    } else {
                                        var sumValue = 0.0
                                        for (row in 0..<Vchain.numRows) sumValue += Vchain.get(row, r) * Wchain_s.get(
                                            row,
                                            r)
                                        Cchain_s.set(0, r, sumValue)
                                        Xchain_s.set(0,
                                            r,
                                            omicron * Nchain_s.get(0, r) / Cchain_s.get(0,
                                                r) + (1 - omicron) * Xchain_s_1.get(0, r))
                                    }
                                }
                                for (k in 0..<M) {
                                    //Rchain_s(k,r) = Vchain(k,r) * Wchain_s(k,r); NOT USED
                                    //Tchain_s(k,r) = Xchain_s(r) * Vchain(k,r); NOT USED
                                    Qchain_s.set(k,
                                        r,
                                        omicron * Xchain_s.get(0, r) * Vchain.get(k, r) * Wchain_s.get(k,
                                            r) + (1 - omicron) * Qchain_s_1.get(k, r))
                                    Uchain_s.set(k,
                                        r,
                                        omicron * Vchain.get(k, r) * STeff_s.get(k, r) * Xchain_s.get(0,
                                            r) + (1 - omicron) * Uchain_s_1.get(k, r))
                                }
                            }
                        }

                        when (options.method) {
                            "default", "lin" -> {
                                var k = 0
                                while (k < M) {
                                    for (r in nnzclasses) {
                                        if (Double.isFinite(Nchain.get(0, r)) && Nchain_s.get(0, r) > 0) gamma.get(r)
                                            .set(s,
                                                k,
                                                Qchain_s_1.get(k, r) / Nchain_s.get(0, r) - QchainOuter_1.get(k,
                                                    r) / Nchain.get(0, r))
                                    }
                                    k++
                                }
                            }

                            else -> {
                                var k = 0
                                while (k < M) {
                                    gamma.get(0)
                                        .set(s, k, Qchain_s_1.sumRows(k) / (Nt - 1) - QchainOuter_1.sumRows(k) / Nt)
                                    k++
                                }
                            }
                        }

                        for (r in nnzclasses) {
                            tau.set(s, r, Xchain_s_1.get(0, r) - XchainOuter_1.get(0, r))
                        }
                    }
                }
            }
        }

        var inner_iter = 0.0
        var Qchain_inner = Qchain.copy()
        var Xchain_inner = Xchain.copy()
        var Uchain_inner = Uchain.copy()

        var count_inf = 0
        var Nt_inner = 0.0
        val deltaclass = Matrix(Nchain.numRows, Nchain.numCols)
        val ocl: MutableList<Int?> = ArrayList<Int?>()
        val ccl: MutableList<Int> = ArrayList<Int>()
        val nnzclasses_inner: MutableList<Int> = ArrayList<Int>()
        for (col in 0..<Nchain.numCols) {
            val nchainValue = Nchain.get(0, col)
            if (Double.isInfinite(nchainValue)) {
                count_inf++
                deltaclass.set(0, col, 1.0)
                ocl.add(col)
            } else {
                Nt_inner += nchainValue
                deltaclass.set(0, col, (nchainValue - 1) / nchainValue)
                if (nchainValue > 0) ccl.add(col)
            }

            if (nchainValue > 0) nnzclasses_inner.add(col)
        }
        val delta = if (count_inf == Nchain.numCols) 1.0 else (Nt_inner - 1) / Nt_inner

        //Use List<Integer> instead of JLineMatrix to store the index to avoid type transformation
        val nnzclasses_eprio: MutableMap<Int?, MutableList<Int?>?> =
            HashMap<Int?, MutableList<Int?>?>(nnzclasses_inner.size)
        val nnzclasses_hprio: MutableMap<Int?, MutableList<Int>?> =
            HashMap<Int?, MutableList<Int>?>(nnzclasses_inner.size)
        val nnzclasses_ehprio: MutableMap<Int?, MutableList<Int>?> =
            HashMap<Int?, MutableList<Int>?>(nnzclasses_inner.size)
        for (r in nnzclasses_inner) {
            val prio = sn.classprio.get(0, r)
            val eprio_list: MutableList<Int?> = ArrayList<Int?>()
            val hprio_list: MutableList<Int?> = ArrayList<Int?>()
            for (i in 0..<sn.classprio.numCols) {
                if (Double.compare(prio, sn.classprio.get(0, i)) == 0) eprio_list.add(i)
                else if (Double.compare(prio, sn.classprio.get(0, i)) < 0) hprio_list.add(i)
            }
            val eprio_common: MutableList<Int?> = ArrayList<Int?>(nnzclasses_inner)
            val hprio_common: MutableList<Int?> = ArrayList<Int?>(nnzclasses_inner)
            eprio_common.retainAll(eprio_list)
            hprio_common.retainAll(hprio_list)
            nnzclasses_eprio.put(r, eprio_common)
            @Suppress("UNCHECKED_CAST")
            nnzclasses_hprio.put(r, hprio_common as MutableList<Int>?)

            val ehprio_common: MutableSet<Int> = LinkedHashSet<Int>(eprio_common)
            ehprio_common.addAll(hprio_common)
            @Suppress("UNCHECKED_CAST")
            nnzclasses_ehprio.put(r, ArrayList<Int?>(ehprio_common) as MutableList<Int>?)
        }
        while ((inner_iter < 2 || Qchain_inner.sub(Qchain)
                .elementMaxAbs() > tol) && (inner_iter <= FastMath.sqrt(options.iter_max.toDouble()))) {
            inner_iter++

            Qchain_inner = Qchain.copy()
            Xchain_inner = Xchain.copy()
            Uchain_inner = Uchain.copy()

            val ret = solver_amvald_forward(gamma,
                tau,
                Qchain_inner,
                Xchain_inner,
                Uchain_inner,
                STchain,
                Vchain,
                Nchain,
                SCVchain,
                Nt_inner,
                delta,
                deltaclass,
                ocl,
                ccl,
                nnzclasses_inner,
                nnzclasses_eprio,
                nnzclasses_hprio,
                nnzclasses_ehprio,
                M,
                K,
                nservers,
                schedparam,
                lldscaling,
                cdscaling,
                ljdscaling,
                ljdcutoffs,
                ljcdscalingIndexed,
                ljcdcutoffsIndexed,
                sched,
                stations,
                options)
            totiter++
            if (totiter > options.iter_max) {
                break
            }
            val Wchain = ret.getLeft()
            STeff = ret.getRight()

            for (r in nnzclasses) {
                if (Wchain.sumCols(r) == 0.0) {
                    Xchain.remove(0, r)
                } else {
                    if (Double.isInfinite(Nchain.get(0, r))) {
                        var sumValue = 0.0
                        for (i in 0..<Vchain.numRows) sumValue += Vchain.get(i, r) * Wchain.get(i, r)
                        Cchain_s.set(0, r, sumValue)
                    } else if (Nchain.get(0, r) == 0.0) {
                        Xchain.remove(0, r)
                        Cchain_s.remove(0, r)
                    } else {
                        var sumValue = 0.0
                        for (i in 0..<Vchain.numRows) sumValue += Vchain.get(i, r) * Wchain.get(i, r)
                        Cchain_s.set(0, r, sumValue)
                        Xchain.set(0,
                            r,
                            omicron * Nchain.get(0, r) / Cchain_s.get(0, r) + (1 - omicron) * Xchain_inner.get(0, r))
                    }
                }
                for (k in 0..<M) {
                    //Rchain(k,r) = Vchain(k,r) * Wchain(k,r); NOT USED
                    Qchain.set(k,
                        r,
                        omicron * Xchain.get(0, r) * Vchain.get(k, r) * Wchain.get(k,
                            r) + (1 - omicron) * Qchain_inner.get(k, r))
                    Tchain.set(k, r, Xchain.get(0, r) * Vchain.get(k, r))
                    Uchain.set(k,
                        r,
                        omicron * Vchain.get(k, r) * STeff.get(k, r) * Xchain.get(0,
                            r) + (1 - omicron) * Uchain_inner.get(k, r))
                }
            }
        }
    }

    for (k in 0..<M) {
        for (r in 0..<K) {
            if (Vchain.get(k, r) * STeff!!.get(k, r) > 0) {
                when (sn.sched.get(sn.stations.get(k))) {
                    SchedStrategy.FCFS, SchedStrategy.SIRO, SchedStrategy.PS, SchedStrategy.LCFSPR, SchedStrategy.DPS, SchedStrategy.HOL -> if (Uchain.sumRows(
                            k) > 1) {
                        var sum_vchain_steff_xchain_k = 0.0
                        var i = 0
                        while (i < Vchain.numCols) {
                            sum_vchain_steff_xchain_k += Vchain.get(k, i) * STeff.get(k, i) * Xchain.get(0, i)
                            i++
                        }
                        Uchain.set(k,
                            r,
                            (min(Uchain.sumRows(k), 1.0) * Vchain.get(k, r) * STeff.get(k, r) * Xchain.get(0,
                                r)) / sum_vchain_steff_xchain_k)
                    }

                    else -> {}
                }
            }
        }
    }

    for (k in 0..<M) {
        for (r in 0..<K)
            if (Qchain.get(k, r) < GlobalConstants.Zero)  {
                Rchain.set(k, r, 0.0)
            } else {
                Rchain.set(k, r, Qchain.get(k, r) / Tchain.get(k, r))
            }
    }
    Xchain.apply(kotlin.Double.Companion.POSITIVE_INFINITY, 0.0, "equal")
    Xchain.apply(kotlin.Double.Companion.NaN, 0.0, "equal")
    Uchain.apply(kotlin.Double.Companion.POSITIVE_INFINITY, 0.0, "equal")
    Uchain.apply(kotlin.Double.Companion.NaN, 0.0, "equal")
    Rchain.apply(kotlin.Double.Companion.POSITIVE_INFINITY, 0.0, "equal")
    Rchain.apply(kotlin.Double.Companion.NaN, 0.0, "equal")

    for (col in 0..<K) {
        if (Nchain.get(0, col) == 0.0) {
            Xchain.remove(0, col)
            for (row in 0..<M) {
                Uchain.remove(row, col)
                Rchain.remove(row, col)
                Tchain.remove(row, col)
            }
        }
    }

    var ret: snDeaggregateChainResults? = null
    if ((sn.lldscaling == null || sn.lldscaling.isEmpty) && (sn.cdscaling == null || sn.cdscaling.size == 0) && !snHasJointDependence(sn)) ret =
        snDeaggregateChainResults(sn, Lchain, null, STchain, Vchain, alpha, null, null, Rchain, Tchain, null, Xchain)
    else ret =
        snDeaggregateChainResults(sn, Lchain, null, STchain, Vchain, alpha, null, Uchain, Rchain, Tchain, null, Xchain)

    val ccl: MutableList<Int?> = ArrayList<Int?>()
    for (i in 0..<Nchain.numCols) {
        if (Double.isFinite(Nchain.get(0, i))) ccl.add(i)
    }
    val Nclosed = Matrix(1, ccl.size)
    val Xclosed = Matrix(1, ccl.size)
    for (i in ccl.indices) {
        Nclosed.set(0, i, Nchain.get(0, ccl.get(i)!!))
        Xclosed.set(0, i, Xchain.get(0, ccl.get(i)!!))
    }
    var lG = 0.0
    for (i in ccl.indices) {
        if (Xclosed.get(0, i) > options.tol) lG += -Nclosed.get(0, i) * FastMath.log(Xclosed.get(0, i))
    }

    val endTime = System.nanoTime()
    val runTime = endTime - startTime

    val result = MVAResult()
    result.method = options.method
    result.QN = ret.Q
    result.RN = ret.R
    result.XN = ret.X
    result.UN = ret.U
    result.TN = ret.T
    result.CN = ret.C
    result.runtime = runTime / 1000000000.0
    result.logNormConstAggr = lG
    return result
}

fun solver_amvald_forward(gamma: MutableList<Matrix>?,
                          tau: Matrix,
                          Qchain_in: Matrix,
                          Xchain_in: Matrix,
                          Uchain_in: Matrix,
                          STchain_in: Matrix,
                          Vchain_in: Matrix,
                          Nchain_in: Matrix,
                          SCVchain_in: Matrix,
                          Nt: kotlin.Double,
                          delta: kotlin.Double,
                          deltaclass: Matrix,
                          ocl: MutableList<Int?>,
                          ccl: MutableList<Int>,
                          nnzclasses: MutableList<Int>,
                          nnzclasses_eprio: MutableMap<Int?, MutableList<Int?>?>?,
                          nnzclasses_hprio: MutableMap<Int?, MutableList<Int>?>,
                          nnzclasses_ehprio: MutableMap<Int?, MutableList<Int>?>,
                          M: Int,
                          K: Int,
                          nservers: Matrix,
                          schedparam: Matrix,
                          lldscaling_in: Matrix?,
                          cdscaling: Map<Station?, SerializableFunction<Matrix?, kotlin.Double>>?,
                          ljdscaling: Map<Station, Matrix>?,
                          ljdcutoffs: Map<Station, Matrix>?,
                          ljcdscaling: List<List<Matrix>?>?,
                          ljcdcutoffs: List<Matrix?>?,
                          sched: List<SchedStrategy>,
                          stations: List<Station?>,
                          options: SolverOptions): Pair<Matrix, Matrix> {
    var gamma = gamma
    var lldscaling: Matrix? = null

    if (gamma!!.size == 0) {
        @Suppress("UNCHECKED_CAST")
        gamma = ArrayList<Matrix?>() as MutableList<Matrix>?
        gamma!!.add(Matrix(K, M))
    }

    /* Evaluate lld and cd correction factors */
    val totArvlQlenSeenByOpen = Matrix(K, M)
    val interpTotArvlQlen = Matrix(M, 1)
    val totArvlQlenSeenByClosed = Matrix(M, K)
    val stationaryQlen = Matrix(M, K)
    val selfArvlQlenSeenByClosed = Matrix(M, K)
    for (k in 0..<M) {
        var sum_qchain_k_nnz = 0.0
        for (r in nnzclasses) sum_qchain_k_nnz += Qchain_in.get(k, r)

        interpTotArvlQlen.set(k, 0, delta * sum_qchain_k_nnz)
        for (r in nnzclasses) {
            selfArvlQlenSeenByClosed.set(k, r, deltaclass.get(r) * Qchain_in.get(k, r))
            if (sched[k] == SchedStrategy.HOL) {
                var sum_qchain_k_ehprio = 0.0
                for (i in nnzclasses_ehprio.get(r)!!) sum_qchain_k_ehprio += Qchain_in.get(k, i)
                totArvlQlenSeenByOpen.set(r, k, sum_qchain_k_ehprio)
                totArvlQlenSeenByClosed.set(k,
                    r,
                    deltaclass.get(r) * Qchain_in.get(k, r) + sum_qchain_k_ehprio - Qchain_in.get(k, r))
            } else {
                totArvlQlenSeenByOpen.set(r, k, sum_qchain_k_nnz)
                totArvlQlenSeenByClosed.set(k,
                    r,
                    deltaclass.get(r) * Qchain_in.get(k, r) + sum_qchain_k_nnz - Qchain_in.get(k, r))
            }
            stationaryQlen.set(k, r, Qchain_in.get(k, r))
        }
    }

    if (lldscaling_in == null || lldscaling_in.isEmpty) {
        //lldscaling = new Matrix(M, (int) FastMath.ceil(Nt)); // fails on LayeredExamples.ex3
        lldscaling = Matrix(M, 1)
        lldscaling.fill(1.0)
    } else {
        lldscaling = lldscaling_in.copy()
    }
    var lldterm = pfqn_lldfun(interpTotArvlQlen.elementIncrease(1.0), lldscaling, null)

    val cdterm = Matrix(M, K)
    cdterm.fill(1.0)
    for (r in nnzclasses) {
        if (!(cdscaling == null || cdscaling.size == 0)) {
            val cdscalingList = (0 until M).map { i -> cdscaling?.get(stations[i]) }
            if (Double.isFinite(Nchain_in.get(0,
                    r))) Matrix.extract(pfqn_cdfun(selfArvlQlenSeenByClosed.elementIncrease(1.0),
                cdscalingList,
                M), 0, M, 0, 1, cdterm, 0, r)
            else Matrix.extract(pfqn_cdfun(stationaryQlen.elementIncrease(1.0), cdscalingList, M),
                0,
                M,
                0,
                1,
                cdterm,
                0,
                r)
        }
    }

    // Joint dependence correction terms
    // Both LJD and LJCD are indexed by population vector n = (n1, n2, ..., nK)
    // - LJD: scaling = table[idx(n)] - same scaling applied to all classes
    // - LJCD: scaling_r = table_r[idx(n)] - each class r has its own table
    // If LJCD is set for a station, it overrides LJD for that station.
    val ljdterm = Matrix(M, K)
    ljdterm.fill(1.0)
    val hasLjcd = ljcdscaling != null && ljcdscaling.isNotEmpty()
    val hasLjd = ljdscaling != null && ljdscaling.isNotEmpty()

    for (k in 0 until M) {
        val station = stations[k]
        if (station == null) continue

        // Get population vector from queue lengths at station k
        val nvec = Matrix(1, K)
        for (r in 0 until K) {
            nvec.set(0, r, Qchain_in.get(k, r))
        }

        // Check LJCD first (indexed by station index k, then class index r)
        val ljcdScalingK = if (hasLjcd) ljcdscaling!![k] else null
        val ljcdCutoffsK = if (hasLjcd && ljcdcutoffs != null) ljcdcutoffs[k] else null

        if (ljcdScalingK != null && ljcdCutoffsK != null) {
            // LJCD: per-class scaling - each class has its own table
            // Clamp to cutoffs
            val nClamped = Matrix(1, K)
            for (r in 0 until K) {
                nClamped.set(0, r, FastMath.min(nvec.get(r), ljcdCutoffsK.get(r)))
            }
            val idx = ljd_linearize(nClamped, ljcdCutoffsK)
            for (r in 0 until K) {
                val classScaling = ljcdScalingK[r]
                if (!classScaling.isEmpty && idx < classScaling.length()) {
                    ljdterm.set(k, r, classScaling.get(idx))
                }
            }
        } else if (hasLjd && ljdscaling!!.containsKey(station)) {
            // LJD: single scaling applied to all classes
            val scaling = ljdscaling[station]
            val cutoffs = ljdcutoffs?.get(station)
            if (scaling != null && !scaling.isEmpty && cutoffs != null) {
                // Clamp to cutoffs
                val nClamped = Matrix(1, K)
                for (r in 0 until K) {
                    nClamped.set(0, r, FastMath.min(nvec.get(r), cutoffs.get(r)))
                }
                val idx = ljd_linearize(nClamped, cutoffs)
                if (idx < scaling.length()) {
                    val scalingValue = scaling.get(idx)
                    for (r in 0 until K) {
                        ljdterm.set(k, r, scalingValue)  // Same value for all classes
                    }
                }
            }
        }
    }

    var msterm = Matrix(0, 0)
    when (options.config.multiserver) {
        "softmin" -> when (options.method) {
            "default", "amva_lin", "lin", "amva_qdlin", "qdlin" -> {
                val g = Matrix(ccl.size, M)
                for (r in ccl) {
                    val param = ((Nt - 1.0) / Nt) * Nchain_in.get(0, r)
                    val gamma_r = gamma.get(r)
                    var i = 0
                    while (i < ccl.size) {
                        val idx = ccl.get(i)
                        var j = 0
                        while (j < M) {
                            g.set(i, j, g.get(i, j) + param * gamma_r.get(idx, j))
                            j++
                        }
                        i++
                    }
                }
                val input = interpTotArvlQlen.add(1.0, g.meanCol().transpose())
                msterm = pfqn_lldfun(input.elementIncrease(1.0), Matrix(0, 0), nservers)
            }

            else -> {
                var g = Matrix(ccl.size, M)
                var i = 0
                while (i < ccl.size) {
                    val idx = ccl.get(i)
                    var j = 0
                    while (j < M) {
                        g.set(i, j, g.get(i, j) + (Nt - 1.0) * gamma.get(0).get(idx, j))
                        j++
                    }
                    i++
                }
                interpTotArvlQlen.add(1.0, g.meanCol().transpose())
                //CommonOps_DSCC.add(1, interpTotArvlQlen, 1, CommonOps_DSCC.transpose(g.meanCol(), null, null), input, null, null);
                msterm = pfqn_lldfun(interpTotArvlQlen.elementIncrease(1.0), Matrix(0, 0), nservers)
            }
        }

        "seidmann" -> {
            nservers.divide(1.0, msterm, false)
            msterm.apply(0.0, 1.0, "equal")
        }

        "default" -> {
            when (options.method) {
                "default", "amva_lin", "lin", "amva_qdlin", "qdlin" -> {
                    val g = Matrix(ccl.size, M)
                    for (r in ccl) {
                        val param = ((Nt - 1.0) / Nt) * Nchain_in.get(0, r)
                        val gamma_r = gamma.get(r)
                        var i = 0
                        while (i < ccl.size) {
                            val idx = ccl.get(i)
                            var j = 0
                            while (j < M) {
                                g.set(i, j, g.get(i, j) + param * gamma_r.get(idx, j))
                                j++
                            }
                            i++
                        }
                    }
                    val input = interpTotArvlQlen.add(1.0, g.meanCol().transpose())
                    //CommonOps_DSCC.add(1, interpTotArvlQlen, 1, CommonOps_DSCC.transpose(g.meanCol(), null, null), input, null, null);
                    msterm = pfqn_lldfun(input.elementIncrease(1.0), Matrix(0, 0), nservers)
                }

                else -> {
                    var g = Matrix(1, M)
                    for (r in ccl) {
                        var j = 0
                        while (j < M) {
                            g.set(0, j, g.get(0, j) + (Nt - 1.0) * gamma.get(0).get(r, j))
                            j++
                        }
                    }
                    msterm =
                        pfqn_lldfun(interpTotArvlQlen.elementIncrease(1 + g.meanRow().value()), Matrix(0, 0), nservers)
                }
            }
            var i = 0
            while (i < M) {
                val schedStrategy: SchedStrategy = sched[i]
                if (schedStrategy == SchedStrategy.FCFS || schedStrategy == SchedStrategy.SIRO || schedStrategy == SchedStrategy.LCFSPR) msterm.set(
                    i,
                    0,
                    1.0 / nservers.get(i, 0))
                i++
            }
        }

        else -> throw RuntimeException("nrecognize multiserver approximation method")
    }

    val Wchain = Matrix(M, K)
    val STeff = Matrix(STchain_in.numRows, STchain_in.numCols)
    lldterm = lldterm.repmat(1, K)
    // ljdterm is already (M x K) - contains either LJD or LJCD values
    for (r in nnzclasses) {
        for (k in 0..<M) {
            STeff.set(k, r, STchain_in.get(k, r) * lldterm.get(k, r) * msterm.get(k, 0) * cdterm.get(k, r) * ljdterm.get(k, r))
        }
    }

    /* if amva.qli or amva.fli, update now totArvlQlenSeenByClosed with STeff */
    when (options.method) {
        "amva_qli", "qli" -> {
            val infset: MutableList<Int> = ArrayList<Int>()
            var i = 0
            while (i < M) {
                if (sched[i] == SchedStrategy.INF) infset.add(i)
                i++
            }

            var k = 0
            while (k < M) {
                if (sched[k] == SchedStrategy.HOL) {
                    for (r in nnzclasses) {
                        //Calculate sum_Qchain_k_r_ephrio;
                        var sum_Qchain_k_r_ephrio = 0.0
                        for (i in nnzclasses_ehprio.get(r)!!) sum_Qchain_k_r_ephrio += Qchain_in.get(k, i)

                        if (abs(Nchain_in.get(0, r) - 1) < 1e-20) {
                            totArvlQlenSeenByClosed.set(k, r, sum_Qchain_k_r_ephrio - Qchain_in.get(k, r))
                        } else {
                            val qlinum = STeff.get(k, r) * (1 + sum_Qchain_k_r_ephrio - Qchain_in.get(k, r))
                            var qliden = 0.0
                            for (i in infset) qliden += STeff.get(i, r)
                            var m = 0
                            while (m < M) {
                                var sum_Qchain_m_r_ephrio = 0.0
                                for (i in nnzclasses_ehprio.get(r)!!) sum_Qchain_m_r_ephrio += Qchain_in.get(m, i)
                                qliden += STeff.get(m, r) * (1 + sum_Qchain_m_r_ephrio - Qchain_in.get(m, r))
                                m++
                            }
                            totArvlQlenSeenByClosed.set(k,
                                r,
                                sum_Qchain_k_r_ephrio - (1 / (Nchain_in.get(0, r) - 1)) * (Qchain_in.get(k,
                                    r) - qlinum / qliden))
                        }
                    }
                } else {
                    for (r in nnzclasses) {
                        //Calculate sum_Qchain_k_r_ephrio;
                        var sum_Qchain_k_r_nnzclasses = 0.0
                        for (i in nnzclasses) sum_Qchain_k_r_nnzclasses += Qchain_in.get(k, i)

                        if (abs(Nchain_in.get(0, r) - 1) < 1e-20) {
                            totArvlQlenSeenByClosed.set(k, r, sum_Qchain_k_r_nnzclasses - Qchain_in.get(k, r))
                        } else {
                            val qlinum = STeff.get(k, r) * (1 + sum_Qchain_k_r_nnzclasses - Qchain_in.get(k, r))
                            var qliden = 0.0
                            for (i in infset) qliden += STeff.get(i, r)
                            var m = 0
                            while (m < M) {
                                var sum_Qchain_m_r_nnzclasses = 0.0
                                for (i in nnzclasses) sum_Qchain_m_r_nnzclasses += Qchain_in.get(m, i)
                                qliden += STeff.get(m, r) * (1 + sum_Qchain_m_r_nnzclasses - Qchain_in.get(m, r))
                                m++
                            }
                            totArvlQlenSeenByClosed.set(k,
                                r,
                                sum_Qchain_k_r_nnzclasses - (1 / (Nchain_in.get(0, r) - 1)) * (Qchain_in.get(k,
                                    r) - qlinum / qliden))
                        }
                    }
                }
                k++
            }
        }

        "amva_fli", "fli" -> {
            var infset = ArrayList<Int>()
            var i = 0
            while (i < M) {
                if (sched[i] == SchedStrategy.INF) infset.add(i)
                i++
            }

            var k = 0
            while (k < M) {
                if (sched[k] == SchedStrategy.HOL) {
                    for (r in nnzclasses) {
                        //Calculate sum_Qchain_k_r_ephrio;
                        var sum_Qchain_k_r_ephrio = 0.0
                        for (i in nnzclasses_ehprio.get(r)!!) sum_Qchain_k_r_ephrio += Qchain_in.get(k, i)

                        if (abs(Nchain_in.get(0, r) - 1) < 1e-20) {
                            totArvlQlenSeenByClosed.set(k, r, sum_Qchain_k_r_ephrio - Qchain_in.get(k, r))
                        } else {
                            val qlinum = STeff.get(k, r) * (1 + sum_Qchain_k_r_ephrio - Qchain_in.get(k, r))
                            var qliden = 0.0
                            for (i in infset) qliden += STeff.get(i, r)
                            var m = 0
                            while (m < M) {
                                var sum_Qchain_m_r_ephrio = 0.0
                                for (i in nnzclasses_ehprio.get(r)!!) sum_Qchain_m_r_ephrio += Qchain_in.get(m, i)
                                qliden += STeff.get(m, r) * (1 + sum_Qchain_m_r_ephrio - Qchain_in.get(m, r))
                                m++
                            }
                            totArvlQlenSeenByClosed.set(k,
                                r,
                                sum_Qchain_k_r_ephrio - (2 / Nchain_in.get(0, r)) * (Qchain_in.get(k,
                                    r) + qlinum / qliden))
                        }
                    }
                } else {
                    for (r in nnzclasses) {
                        //Calculate sum_Qchain_k_r_ephrio;
                        var sum_Qchain_k_r_nnzclasses = 0.0
                        for (i in nnzclasses) sum_Qchain_k_r_nnzclasses += Qchain_in.get(k, i)

                        if (abs(Nchain_in.get(0, r) - 1) < 1e-20) {
                            totArvlQlenSeenByClosed.set(k, r, sum_Qchain_k_r_nnzclasses - Qchain_in.get(k, r))
                        } else {
                            val qlinum = STeff.get(k, r) * (1 + sum_Qchain_k_r_nnzclasses - Qchain_in.get(k, r))
                            var qliden = 0.0
                            for (i in infset) qliden += STeff.get(i, r)
                            var m = 0
                            while (m < M) {
                                var sum_Qchain_m_r_nnzclasses = 0.0
                                for (i in nnzclasses) sum_Qchain_m_r_nnzclasses += Qchain_in.get(m, i)
                                qliden += STeff.get(m, r) * (1 + sum_Qchain_m_r_nnzclasses - Qchain_in.get(m, r))
                                m++
                            }
                            totArvlQlenSeenByClosed.set(k,
                                r,
                                sum_Qchain_k_r_nnzclasses - (2 / Nchain_in.get(0, r)) * (Qchain_in.get(k,
                                    r) + qlinum / qliden))
                        }
                    }
                }
                k++
            }
        }

        else -> {}
    }

    /* Compute response time */
    for (r in nnzclasses) {
        val sd: MutableList<Int> = ArrayList<Int>(nnzclasses)
        val sdprio: MutableList<Int?> = ArrayList<Int?>(nnzclasses_ehprio.get(r))
        sd.remove(r)
        sdprio.remove(r)

        for (k in 0..<M) {
            when (sched[k]) {
                SchedStrategy.INF -> Wchain.set(k, r, STeff.get(k, r))
                SchedStrategy.PS -> when (options.method) {
                    "def", "amva", "amva_qd", "amva_qdamva", "qd", "qdamva", "lin", "qdlin" -> if (options.config.multiserver == "seidmann") {
                        val multiServerTerm = STeff.get(k, r) * (nservers.get(k, 0) - 1)
                        if (ocl.contains(r)) {
                            Wchain.set(k, r, multiServerTerm + STeff.get(k, r) * (1 + totArvlQlenSeenByOpen.get(r, k)))
                        } else {
                            when (options.method) {
                                "default", "amva_lin", "lin", "amva_qdlin", "qdlin" -> {
                                    //Nchain_in(ccl)*permute(gamma(r,k,ccl),3:-1:1)
                                    var tmp = 0.0
                                    for (c in ccl) tmp += Nchain_in.get(0, c) * gamma.get(c).get(r, k)
                                    Wchain.set(k,
                                        r,
                                        multiServerTerm + STeff.get(k, r) * (1 + interpTotArvlQlen.get(k, 0) + tmp - gamma.get(r)
                                            .get(r, k)))
                                }

                                else -> Wchain.set(k,
                                    r,
                                    multiServerTerm + STeff.get(k, r) * (1 + interpTotArvlQlen.get(k, 0) + (Nt - 1) * gamma.get(0)
                                        .get(r, k)))
                            }
                        }
                    } else {
                        if (ocl.contains(r)) {
                            Wchain.set(k, r, STeff.get(k, r) * (1 + totArvlQlenSeenByOpen.get(r, k)))
                        } else {
                            when (options.method) {
                                "default", "amva_lin", "lin", "amva_qdlin", "qdlin" -> {
                                    var tmp = 0.0
                                    for (c in ccl) tmp += Nchain_in.get(0, c) * gamma.get(c).get(r, k)
                                    Wchain.set(k,
                                        r,
                                        Wchain.get(k, r) + STeff.get(k, r) * (1 + interpTotArvlQlen.get(k,
                                            0) + tmp - gamma.get(r).get(r, k)))
                                }

                                else -> Wchain.set(k,
                                    r,
                                    STeff.get(k, r) * (1 + interpTotArvlQlen.get(k, 0) + (Nt - 1) * gamma.get(0)
                                        .get(r, k)))
                            }
                        }
                    }

                    else -> if (options.config.multiserver == "seidmann") {
                        val multiServerTerm = STeff.get(k, r) * (nservers.get(k, 0) - 1)
                        if (ocl.contains(r)) {
                            Wchain.set(k, r, multiServerTerm + STeff.get(k, r) * (1 + totArvlQlenSeenByOpen.get(r, k)))
                        } else {
                            when (options.method) {
                                "default", "amva_lin", "lin", "amva_qdlin", "qdlin" -> {
                                    var tmp = 0.0
                                    for (c in ccl) tmp += Nchain_in.get(0, c) * gamma.get(c).get(r, k)
                                    Wchain.set(k,
                                        r,
                                        multiServerTerm + STeff.get(k, r) * (1 + totArvlQlenSeenByClosed.get(k,
                                            0) + tmp - gamma.get(r).get(r, k)))
                                }

                                else -> Wchain.set(k,
                                    r,
                                    multiServerTerm + STeff.get(k, r) * (1 + totArvlQlenSeenByClosed.get(k,
                                        0) + (Nt - 1) * gamma.get(0).get(r, k)))
                            }
                        }
                    } else {
                        val currentWchain = Wchain.get(k, r)
                        if (ocl.contains(r)) {
                            Wchain.set(k, r, currentWchain + STeff.get(k, r) * (1 + totArvlQlenSeenByOpen.get(r, k)))
                        } else {
                            when (options.method) {
                                "default", "amva_lin", "lin", "amva_qdlin", "qdlin" -> {
                                    var tmp = 0.0
                                    for (c in ccl) tmp += Nchain_in.get(0, c) * gamma.get(c).get(r, k)
                                    Wchain.set(k,
                                        r,
                                        currentWchain + STeff.get(k, r) * (1 + totArvlQlenSeenByClosed.get(k,
                                            0) + tmp - gamma.get(r).get(r, k)))
                                }

                                else -> Wchain.set(k,
                                    r,
                                    currentWchain + STeff.get(k, r) * (1 + totArvlQlenSeenByClosed.get(k,
                                        0) + (Nt - 1) * gamma.get(0).get(r, k)))
                            }
                        }
                    }
                }

                SchedStrategy.DPS -> {
                    if (nservers.get(k,
                            0) > 1) throw RuntimeException("Multi-server DPS not supported yet in AMVA solver")

                    Wchain.set(k, r, STeff.get(k, r) * (1 + selfArvlQlenSeenByClosed.get(k, r)))
                    for (s in sd) {
                        if (schedparam.get(k, s) == schedparam.get(k, r)) Wchain.set(k,
                            r,
                            Wchain.get(k, r) + STeff.get(k, r) * stationaryQlen.get(k, s))
                        else if (schedparam.get(k, s) / schedparam.get(k,
                                r) <= kotlin.Double.Companion.POSITIVE_INFINITY) Wchain.set(k,
                            r,
                            Wchain.get(k, r) + STeff.get(k, r) * stationaryQlen.get(k, s) * schedparam.get(k,
                                s) / schedparam.get(k, r))
                    }
                }

                SchedStrategy.FCFS, SchedStrategy.SIRO, SchedStrategy.LCFSPR -> {
                    if (STeff.get(k, r) <= 0) break

                    //Uchain_r = Uchain_in ./ repmat(Xchain_in,M,1) .* (repmat(Xchain_in,M,1) + repmat(tau(r,:),M,1));
                    val Uchain_r = Matrix(Uchain_in.numRows, Uchain_in.numCols)
                    var i = 0
                    while (i < Uchain_in.numRows) {
                        var j = 0
                        while (j < Uchain_in.numCols) {
                            Uchain_r.set(i,
                                j,
                                (Uchain_in.get(i, j) / Xchain_in.get(0, j)) * (Xchain_in.get(0, j) + tau.get(r, j)))
                            j++
                        }
                        i++
                    }

                    var Bk = Matrix(1, K)
                    if (nservers.get(k, 0) > 1) {
                        val deltaclass_r = Matrix(Xchain_in.numRows, Xchain_in.numCols)
                        deltaclass_r.fill(1.0)
                        deltaclass_r.set(0, r, deltaclass.get(0, r))
                        //Compute load: deltaclass_r .* Xchain_in .* Vchain_in(k,:) .* STeff(k,:)
                        val BK_tmp = Matrix(1, K)
                        var i = 0
                        while (i < K) {
                            BK_tmp.set(0,
                                i,
                                deltaclass_r.get(0, i) * Xchain_in.get(0, i) * Vchain_in.get(k, i) * STeff.get(k, i))
                            i++
                        }
                        if (BK_tmp.elementSum() < 0.75) {
                            Bk = BK_tmp.copy()
                        } else {
                            Bk = BK_tmp.elementPower(nservers.get(k, 0) - 1)
                        }
                    } else {
                        Bk.fill(1.0)
                    }


                    val hasLjd = ljdscaling != null && ljdscaling.isNotEmpty()
                    if (nservers.get(k,
                            0) == 1.0 && (((lldscaling != null) && !lldscaling.isEmpty) || ((cdscaling != null) && (cdscaling.size != 0)) || hasLjd)) {
                        if (options.config.highvar == "hvmva") {
                            var sum_uchain_r = 0.0
                            var weightedSum = 0.0
                            for (s in ccl) {
                                sum_uchain_r += Uchain_r.get(k, s)
                                weightedSum += STeff.get(k, s) * Uchain_r.get(k, s) * (1.0 + SCVchain_in.get(k, s)) / 2.0
                            }
                            Wchain.set(k, r, weightedSum + STeff.get(k, r) * (1 - sum_uchain_r))
                        } else {
                            Wchain.set(k, r, STeff.get(k, r))
                        }

                        var steff_mult_stationaryQlen = 0.0
                        for (s in sd) steff_mult_stationaryQlen += STeff.get(k, s) * stationaryQlen.get(k, s)

                        if (ocl.contains(r)) {
                            Wchain.set(k,
                                r,
                                Wchain.get(k, r) + (STeff.get(k, r) * stationaryQlen.get(k,
                                    r) + steff_mult_stationaryQlen))
                        } else {
                            //case {'default', 'amva.lin', 'lin', 'amva.qdlin','qdlin'} % Linearizer
                            //    Wchain(k,r) = Wchain(k,r) + (STeff(k,r) * selfArvlQlenSeenByClosed(k,r) + STeff(k,sd)*stationaryQlen(k,sd)') + (STeff(k,ccl).*Nchain(ccl)*permute(gamma(r,k,ccl),3:-1:1) - STeff(k,r)*gamma(r,k,r));
                            Wchain.set(k,
                                r,
                                Wchain.get(k, r) + (STeff.get(k, r) * selfArvlQlenSeenByClosed.get(k,
                                    r) + steff_mult_stationaryQlen))
                        }
                    } else {
                        var steff_mult_stationarQlen_mult_Bk = 0.0
                        for (s in sd) steff_mult_stationarQlen_mult_Bk += STeff.get(k, s) * stationaryQlen.get(k,
                            s) * Bk.get(0, s)

                        if (options.config.multiserver == "softmin") {
                            if (ocl.contains(r)) {
                                Wchain.set(k,
                                    r,
                                    STeff.get(k, r) + STeff.get(k, r) * stationaryQlen.get(k, r) * Bk.get(0,
                                        r) + steff_mult_stationarQlen_mult_Bk)
                            } else {
                                //case {'default', 'amva.lin', 'lin', 'amva.qdlin','qdlin'} % Linearizer
                                //    Wchain(k,r) = Wchain(k,r) + STeff(k,r) * selfArvlQlenSeenByClosed(k,r) * Bk(r) + STeff(k,sd) * (stationaryQlen(k,sd) .* Bk(sd))' + (STeff(k,ccl).*Nchain(ccl)*permute(gamma(r,k,ccl),3:-1:1) - STeff(k,r)*gamma(r,k,r));
                                Wchain.set(k,
                                    r,
                                    STeff.get(k, r) + STeff.get(k, r) * selfArvlQlenSeenByClosed.get(k, r) * Bk.get(0,
                                        r) + steff_mult_stationarQlen_mult_Bk)
                            }
                        } else {
                            if (ocl.contains(r)) {
                                Wchain.set(k,
                                    r,
                                    STeff.get(k, r) * (nservers.get(k, 0) - 1) + STeff.get(k, r) + (STeff.get(k,
                                        r) * deltaclass.get(0, r) * stationaryQlen.get(k, r) * Bk.get(0,
                                        r) + steff_mult_stationarQlen_mult_Bk))
                            } else {
                                //case {'default', 'amva.lin', 'lin', 'amva.qdlin','qdlin'} % Linearizer
                                //    Wchain(k,r) = Wchain(k,r) + (STeff(k,r) * selfArvlQlenSeenByClosed(k,r)*Bk(r) + STeff(k,sd).*Bk(sd)*stationaryQlen(k,sd)') + (STeff(k,ccl).*Nchain(ccl)*permute(gamma(r,k,ccl),3:-1:1) - STeff(k,r)*gamma(r,k,r));
                                Wchain.set(k,
                                    r,
                                    STeff.get(k, r) * (nservers.get(k, 0) - 1) + STeff.get(k, r) + (STeff.get(k,
                                        r) * selfArvlQlenSeenByClosed.get(k, r) * Bk.get(0,
                                        r) + steff_mult_stationarQlen_mult_Bk))
                            }
                        }
                    }
                }

                SchedStrategy.HOL -> {
                    if (STeff.get(k, r) <= 0) break

                    var Uchain_r = Matrix(Uchain_in.numRows, Uchain_in.numCols)
                    var i = 0
                    while (i < Uchain_in.numRows) {
                        var j = 0
                        while (j < Uchain_in.numCols) {
                            Uchain_r.set(i,
                                j,
                                (Uchain_in.get(i, j) / Xchain_in.get(0, j)) * (Xchain_in.get(0, j) + tau.get(r, j)))
                            j++
                        }
                        i++
                    }

                    var prioScaling = 0.0
                    when (options.config.np_priority) {
                        "default", "cl" -> {
                            var UHigherPrio = 0.0
                            for (h in nnzclasses_hprio.get(r)!!) UHigherPrio += Vchain_in.get(k, h) * STeff.get(k,
                                h) * (Xchain_in.get(0, h) - Qchain_in.get(k, h) * tau.get(h))
                            prioScaling = FastMath.min(max(options.tol, 1 - UHigherPrio), 1 - options.tol)
                        }

                        "shadow" -> {
                            var UHigherPrio = 0.0
                            for (h in nnzclasses_hprio.get(r)!!) UHigherPrio += Vchain_in.get(k, h) * STeff.get(k,
                                h) * Xchain_in.get(0, h)
                            prioScaling = FastMath.min(max(options.tol, 1 - UHigherPrio), 1 - options.tol)
                        }
                    }

                    var Bk = Matrix(1, K)
                    if (nservers.get(k, 0) > 1) {
                        //Compute load: deltaclass .* Xchain_in .* Vchain_in(k,:) .* STeff(k,:)
                        val BK_tmp = Matrix(1, K)
                        var i = 0
                        while (i < K) {
                            BK_tmp.set(0,
                                i,
                                deltaclass.get(0, i) * Xchain_in.get(0, i) * Vchain_in.get(k, i) * STeff.get(k, i))
                            i++
                        }

                        if (BK_tmp.elementSum() < 0.75) {
                            when (options.config.multiserver) {
                                "softmin" -> Bk = BK_tmp.copy()
                                "default", "seidmann" -> BK_tmp.divide(nservers.get(k, 0), Bk, true)
                            }
                        } else {
                            when (options.config.multiserver) {
                                "softmin" -> Bk = BK_tmp.elementPower(nservers.get(k, 0))
                                "default", "seidmann" -> {
                                    BK_tmp.divide(nservers.get(k, 0), BK_tmp, true)
                                    Bk = BK_tmp.elementPower(nservers.get(k, 0))
                                }
                            }
                        }
                    } else {
                        Bk.fill(1.0)
                    }

                    val hasLjd2 = ljdscaling != null && ljdscaling.isNotEmpty()
                    if (nservers.get(k,
                            0) == 1.0 && (((lldscaling != null) && !lldscaling.isEmpty) || ((cdscaling != null) && (cdscaling.size != 0)) || hasLjd2)) {
                        if (options.config.highvar == "hvmva") {
                            var sum_uchain_r = 0.0
                            for (s in ccl) sum_uchain_r += Uchain_r.get(k, s)
                            Wchain.set(k, r, (STeff.get(k, r) / prioScaling) * (1 - sum_uchain_r))
                            for (s in ccl) {
                                var UHigherPrio_s = 0.0
                                for (h in nnzclasses_hprio.get(s)!!) UHigherPrio_s += Vchain_in.get(k, h) * STeff.get(k,
                                    h) * (Xchain_in.get(0, h) - Qchain_in.get(k, h) * tau.get(h))
                                val prioScaling_s = FastMath.min(max(options.tol, 1 - UHigherPrio_s), 1 - options.tol)
                                Wchain.set(k,
                                    r,
                                    Wchain.get(k, r) + (STeff.get(k, s) / prioScaling_s) * Uchain_r.get(k,
                                        s) * (1 + SCVchain_in.get(k, s)) / 2)
                            }
                        } else {
                            Wchain.set(k, r, STeff.get(k, r) / prioScaling)
                        }

                        if (ocl.contains(r)) {
                            Wchain.set(k,
                                r,
                                Wchain.get(k, r) + (STeff.get(k, r) * stationaryQlen.get(k, r)) / prioScaling)
                        } else {
                            //case {'default', 'amva.lin', 'lin', 'amva.qdlin','qdlin'} % Linearizer
                            //    %Wchain(k,r) = Wchain(k,r) + (STeff(k,r) * selfArvlQlenSeenByClosed(k,r) + STeff(k,sdprio)*stationaryQlen(k,sdprio)') + (STeff(k,[r,sdprio]).*Nchain([r,sdprio])*permute(gamma(r,k,[r,sdprio]),3:-1:1) - STeff(k,r)*gamma(r,k,r));
                            //    Wchain(k,r) = Wchain(k,r) + (STeff(k,r) * selfArvlQlenSeenByClosed(k,r) - STeff(k,r)*gamma(r,k,r)) / prioScaling;
                            Wchain.set(k,
                                r,
                                Wchain.get(k, r) + (STeff.get(k, r) * selfArvlQlenSeenByClosed.get(k, r)) / prioScaling)
                        }
                    } else {
                        when (options.config.multiserver) {
                            "softmin" -> if (ocl.contains(r)) {
                                Wchain.set(k,
                                    r,
                                    (STeff.get(k, r) / prioScaling) + STeff.get(k, r) * stationaryQlen.get(k,
                                        r) * Bk.get(0, r) / prioScaling)
                            } else {
                                //case {'default', 'amva.lin', 'lin', 'amva.qdlin','qdlin'} % Linearizer
                                //    %Wchain(k,r) = Wchain(k,r) + STeff(k,r) * selfArvlQlenSeenByClosed(k,r) * Bk(r) + STeff(k,sdprio) * (stationaryQlen(k,sdprio) .* Bk(sdprio))' + (STeff(k,[r,sdprio]).*Nchain([r,sdprio])*permute(gamma(r,k,[r,sdprio]),3:-1:1) - STeff(k,r)*gamma(r,k,r));
                                //    Wchain(k,r) = Wchain(k,r) + STeff(k,r) * selfArvlQlenSeenByClosed(k,r) * Bk(r) / prioScaling + (STeff(k,[r]).*Nchain([r])*permute(gamma(r,k,[r]),3:-1:1) - STeff(k,r)*gamma(r,k,r)) / prioScaling;
                                Wchain.set(k,
                                    r,
                                    (STeff.get(k, r) / prioScaling) + STeff.get(k, r) * selfArvlQlenSeenByClosed.get(k,
                                        r) * Bk.get(0, r) / prioScaling)
                            }

                            "seidmann", "default" -> if (ocl.contains(r)) {
                                Wchain.set(k,
                                    r,
                                    (STeff.get(k, r) * (nservers.get(k, 0) - 1) / prioScaling) + (STeff.get(k,
                                        r) / prioScaling) + (STeff.get(k, r) * stationaryQlen.get(k, r) * Bk.get(0,
                                        r)) / prioScaling)
                            } else {
                                //case {'default', 'amva.lin', 'lin', 'amva.qdlin','qdlin'} % Linearizer
                                //    %Wchain(k,r) = Wchain(k,r) + (STeff(k,r) * selfArvlQlenSeenByClosed(k,r)*Bk(r) + STeff(k,sdprio).*Bk(sdprio)*stationaryQlen(k,sdprio)') + (STeff(k,[r,sdprio]).*Nchain([r,sdprio])*permute(gamma(r,k,[r,sdprio]),3:-1:1) - STeff(k,r)*gamma(r,k,r));
                                //    Wchain(k,r) = Wchain(k,r) + STeff(k,r) * selfArvlQlenSeenByClosed(k,r)*Bk(r)/prioScaling + (STeff(k,r).*Nchain(r)*permute(gamma(r,k,r),3:-1:1) - STeff(k,r)*gamma(r,k,r))/prioScaling;
                                Wchain.set(k,
                                    r,
                                    (STeff.get(k, r) * (nservers.get(k, 0) - 1) / prioScaling) + (STeff.get(k,
                                        r) / prioScaling) + (STeff.get(k, r) * selfArvlQlenSeenByClosed.get(k,
                                        r) * Bk.get(0, r)) / prioScaling)
                            }
                        }
                    }
                }

                else -> {}
            }
        }
    }
    return Pair<Matrix, Matrix>(Wchain, STeff)
}

