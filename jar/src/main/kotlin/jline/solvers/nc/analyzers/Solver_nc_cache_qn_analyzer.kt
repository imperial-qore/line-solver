/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.nc.analyzers

import jline.api.cache.*
import jline.api.mc.dtmc_stochcomp
import jline.api.sn.snRefreshVisits
import jline.io.line_error
import jline.io.mfilename
import jline.lang.NetworkStruct
import jline.lang.constant.NodeType
import jline.lang.nodeparam.CacheNodeParam
import jline.solvers.SolverOptions
import jline.solvers.nc.NCResult
import jline.solvers.nc.analyzers.solver_nc_analyzer
import jline.solvers.nc.analyzers.solver_ncld_analyzer
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import java.io.*
import java.lang.Double
import java.util.*
import kotlin.Any
import kotlin.Array
import kotlin.Int


/**
 * NC Analyzer class for solver_nc_cacheqn_analyzer
 */
fun solver_nc_cache_qn_analyzer(sn: NetworkStruct, options: SolverOptions): NCResult {
    var sn = sn
    var snorig = sn
    var res = NCResult()
    try {
        val bos = ByteArrayOutputStream()
        val out = ObjectOutputStream(bos)
        out.writeObject(sn)
        val bis = ByteArrayInputStream(bos.toByteArray())
        val `in` = ObjectInputStream(bis)
        snorig = `in`.readObject() as NetworkStruct
    } catch (e: IOException) {
        line_error(mfilename(object : Any() {}),
            "Could not create a copy of the NetworkStruct in SolverNCCacheQNAnalyzer")
    } catch (e: ClassNotFoundException) {
        line_error(mfilename(object : Any() {}),
            "Could not create a copy of the NetworkStruct in SolverNCCacheQNAnalyzer")
    }
    sn = snorig
    val random = Random(options.seed.toLong())
    val I = sn.nnodes
    val K = sn.nclasses
    val statefulNodes = ArrayList<Int?>()
    for (i in 0..<sn.isstateful.length()) {
        if (sn.isstateful.get(i) == 1.0) {
            statefulNodes.add(i)
        }
    }
    val statefulNodeClasses = Matrix(1, statefulNodes.size * K)
    var idx = 0
    for (ind in statefulNodes) {
        for (i in 0..<K) {
            statefulNodeClasses.set(idx, (ind!! * K + i).toDouble())
            idx++
        }
    }
    var lambda = Matrix(1, K)
    var lambda_1 = Matrix(1, K)
    val caches = ArrayList<Int?>()
    for (i in sn.nodetype.indices) {
        if (sn.nodetype.get(i) == NodeType.Cache) {
            caches.add(i)
        }
    }
    val hitprob = Matrix(sn.nodetype.size, K)
    val missprob = Matrix(sn.nodetype.size, K)

    for (it in 1..options.iter_max) {
        val inputClass = ArrayList<Int?>()
        for (ind in caches) {
            val ch = sn.nodeparam.get(sn.nodes.get(ind!!)) as CacheNodeParam
            val hitClass = ch.hitclass
            val missClass = ch.missclass
            for (i in 0..<hitClass.length()) {
                if (hitClass.get(i) != -1.0) {
                    inputClass.add(i)
                }
            }

            // Solution of the isolated cache
            val m = ch.itemcap
            val n = ch.nitems
            val h = m.length()
            val u = lambda.length()

            if (it == 1) {
                if (n < m.elementSum() + 2) {
                    line_error(mfilename(object : Any() {}),
                        "NC requires the number of items to exceed the cache capacity at least by 2.")
                }
                for (i in inputClass.indices) {
                    lambda_1.set(inputClass.get(i)!!, random.nextDouble())
                }
                lambda = Matrix(lambda_1)
                sn.nodetype.set(ind, NodeType.ClassSwitch)
            }

            val missrate = Matrix(sn.nodetype.size, u)

            val lambda_cache = Array(u) { Matrix(n, h + 1) }
            for (i in 0..<u) {
                lambda_cache[i] = Matrix(n, h + 1)
            }

            for (v in 0..<u) {
                for (k in 0..<n) {
                    for (l in 0..<h + 1) {
                        if (ch.pread.getOrDefault(v, null) != null) {
                            lambda_cache[v]!!.set(k, l, lambda.get(v) * ch.pread.get(v)!!.get(k))
                        }
                    }
                }
            }

            val R = ch.accost
            val gamma = cache_gamma_lp(lambda_cache, R).gamma
            
            if (options.method == "exact") {
                val pij = cache_prob_erec(gamma, m)
                res.method = "exact"
                for (i in 0..<missprob.numCols) {
                    missprob.set(ind, i, 0)
                }
                for (v in 0..<u) {
                    val A = Matrix.extractColumn(lambda_cache[v], 0, null).transpose()
                    val B = Matrix.extractColumn(pij, 0, null)
                    missrate.set(ind, v, A.mult(B).get(0))
                }
            } else {
                // Use cache_miss_spm for non-exact method (consistent with MATLAB)
                // Convert Array<Matrix> to MatrixCell for function call
                val lambdaCell = MatrixCell(lambda_cache)
                val cacheMissResult = cache_miss_spm(gamma, m, lambdaCell)
                res.method = "spm"
                for (i in 0..<missprob.numCols) {
                    missprob.set(ind, i, 0)
                }
                // Extract miss rates directly from MU field
                for (v in 0..<u) {
                    missrate.set(ind, v, cacheMissResult.getMU()[v])
                }
            }
            for (i in 0..<missprob.numCols) {
                val missValue = missrate.get(ind, i) / lambda.get(i)
                if (Double.isNaN(missValue)) {
                    missprob.set(ind, i, 0)
                } else {
                    missprob.set(ind, i, missValue)
                }
            }
            for (i in 0..<hitprob.numCols) {
                val hitValue = 1 - missrate.get(ind, i) / lambda.get(i)
                if (Double.isNaN(hitValue)) {
                    hitprob.set(ind, i, 0)
                } else {
                    hitprob.set(ind, i, hitValue)
                }
            }
            // Bring back isolated model results into the queueing model
            for (i in inputClass.indices) {
                val r: Int = inputClass.get(i)!!
                for (j in 0..<sn.rtnodes.numCols) {
                    sn.rtnodes.set(ind * K + r, j, 0)
                }
                for (jnd in 0..<I) {
                    if (sn.connmatrix.get(ind, jnd) == 1.0) {
                        sn.rtnodes.set(ind * K + r, (jnd * K + hitClass.get(r)).toInt(), hitprob.get(ind, r))
                        sn.rtnodes.set(ind * K + r, (jnd * K + missClass.get(r)).toInt(), missprob.get(ind, r))
                    }
                }
            }
            val statefulNodeClassesList: MutableList<Int?> = ArrayList<Int?>(statefulNodeClasses.length())
            for (i in 0..<statefulNodeClasses.length()) {
                statefulNodeClassesList.add(statefulNodeClasses.get(i).toInt())
            }
            sn.rt = dtmc_stochcomp(sn.rtnodes, statefulNodeClassesList)
        }
        sn = snRefreshVisits(sn, sn.chains, sn.rt, sn.rtnodes)

        when (options.method) {
            else -> if (!(sn.lldscaling == null || sn.lldscaling.isEmpty) || !(sn.cdscaling == null || sn.cdscaling.isEmpty())) {
                res = solver_ncld_analyzer(sn, options)
            } else {
                res = solver_nc_analyzer(sn, options)
            }
        }

        var nodevisits: Matrix? = null
        for (key in sn.nodevisits.keys) {
            if (nodevisits == null) {
                nodevisits = sn.nodevisits.get(key)
            } else {
                nodevisits = nodevisits.add(1.0, sn.nodevisits.get(key))
            }
        }
        for (ind in caches) {
            for (i in inputClass.indices) {
                val r: Int = inputClass.get(i)!!
                var c = 0
                while (c < sn.chains.numRows && sn.chains.get(c, r) == 0.0) {
                    c++
                }
                val inchain = ArrayList<Int?>()
                for (j in 0..<sn.chains.numCols) {
                    if (sn.chains.get(c, j) != 0.0) {
                        inchain.add(j)
                    }
                }
                var sumXN = 0.0
                for (ix in inchain) {
                    sumXN += res.XN.get(ix!!)
                }
                if (sn.refclass.get(c) > -1) {
                    lambda.set(r,
                        sumXN * nodevisits!!.get(ind!!, r) / nodevisits.get(sn.stationToNode.get(sn.refstat.get(r)
                            .toInt()).toInt(), sn.refclass.get(c).toInt()))
                } else {
                    lambda.set(r,
                        sumXN * nodevisits!!.get(ind!!, r) / nodevisits.get(sn.stationToNode.get(sn.refstat.get(r)
                            .toInt()).toInt(), r))
                }
            }
        }
        if (lambda.sub(lambda_1).elementSum() < options.iter_tol) {
            res.it = it
            break
        }
        lambda_1 = lambda
    }
    res.hitProb = hitprob
    res.missProb = missprob
    return res
}
