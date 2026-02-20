/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.nc.analyzers

import jline.api.cache.cache_gamma_lp
import jline.api.cache.cache_miss_is
import jline.api.cache.cache_miss_spm
import jline.api.cache.cache_prob_erec
import jline.api.cache.cache_prob_is
import jline.io.line_error
import jline.io.mfilename
import jline.lang.NetworkStruct
import jline.lang.constant.NodeType
import jline.lang.constant.ReplacementStrategy
import jline.lang.nodeparam.CacheNodeParam
import jline.lang.nodes.Cache
import jline.solvers.SolverOptions
import jline.solvers.nc.NCResult
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import java.lang.Double
import kotlin.Any
import kotlin.Array
import kotlin.math.abs

fun solver_nc_cache_analyzer(sn: NetworkStruct, options: SolverOptions): NCResult {
    val startTime = System.nanoTime()
    val QN = Matrix(sn.nnodes, sn.nclasses)
    val UN = Matrix(sn.nnodes, sn.nclasses)
    val RN = Matrix(sn.nnodes, sn.nclasses)
    var TN = Matrix(sn.nnodes, sn.nclasses)
    val CN = Matrix(1, sn.nclasses)
    val AN = Matrix(sn.nnodes, sn.nclasses)
    val WN = Matrix(sn.nnodes, sn.nclasses)
    val XN = Matrix(1, sn.nclasses)
    val lG = Double.NaN

    var source_ist = -1
    for (i in sn.nodetype.indices) {
        if (sn.nodetype.get(i) == NodeType.Source) {
            source_ist = sn.nodeToStation.get(i).toInt()
            break
        }
    }
    val sourceRate = Matrix(1, sn.rates.numCols)
    for (i in 0..<sn.rates.numCols) {
        if (!Double.isNaN(sn.rates.get(source_ist, i))) {
            sourceRate.set(i, sn.rates.get(source_ist, i))
        }
    }
    TN = Matrix(sn.nodetype.size, sn.nclasses)
    for (i in 0..<sourceRate.numCols) {
        TN.set(source_ist, i, sourceRate.get(i))
    }

    var cache: Cache? = null
    for (i in sn.nodetype.indices) {
        if (sn.nodetype.get(i) == NodeType.Cache) {
            cache = sn.nodes.get(i) as Cache?
            break
        }
    }
    val ch = sn.nodeparam.get(cache) as CacheNodeParam

    val m = ch.itemcap
    val n = ch.nitems // Number of items

    /*
 * if n<m+2
    line_error(mfilename,'NC requires the number of items to exceed the cache capacity at least by 2.');
   end
 */
    val h = m.length() // Number of lists
    val u = sn.nclasses // Number of users
    val lambda = Array(u) { Matrix(n, h + 1) }
    for (i in 0..<u) {
        lambda[i] = Matrix(n, h + 1)
    }

    for (v in 0..<u) {
        for (k in 0..<n) {
            for (l in 0..<h + 1) {
                if (ch.pread.getOrDefault(v, null) != null) {
                    lambda[v]!!.set(k, l, sourceRate.get(v) * ch.pread.get(v)!!.get(k))
                }
            }
        }
    }

    val lambdaCell = MatrixCell(lambda)

    val R = ch.accost
    val gamma = cache_gamma_lp(lambda, R).gamma

    var pij: Matrix? = null
    val missRate = Matrix(1, u)

    val res = NCResult()
    when (options.method) {
        "exact" -> {
            res.method = "exact"
            when (cache!!.replacementStrategy) {
                ReplacementStrategy.RR, ReplacementStrategy.FIFO -> {
                    pij = cache_prob_erec(gamma, m)
                    val newPij = Matrix(pij.numRows, 1 + pij.numCols)
                    var i = 0
                    while (i < pij.numRows) {
                        var j = 0
                        while (j < newPij.numCols) {
                            if (j < 1) {
                                newPij.set(i, j, abs(1 - pij.sumRows(i)))
                            } else {
                                newPij.set(i, j, pij.get(i, j - 1))
                            }
                            j++
                        }
                        i++
                    }
                    pij = newPij
                    var v = 0
                    while (v < u) {
                        missRate.set(v,
                            Matrix.extractColumn(lambda[v], 0, null).transpose()
                                .mult(Matrix.extractColumn(pij, 0, null)).get(0))
                        v++
                    }
                    line_error(mfilename(object : Any() {}),
                        "NC does not support exact solution of the specified cache replacement policy.")
                }

                else -> line_error(mfilename(object : Any() {}),
                    "NC does not support exact solution of the specified cache replacement policy.")
            }
        }

        "sampling" -> {
            res.method = "sampling"
            when (cache!!.replacementStrategy) {
                ReplacementStrategy.RR, ReplacementStrategy.FIFO -> {
                    val samples = options.samples.toInt()
                    val cacheMissIs = cache_miss_is(gamma, m, lambdaCell, samples)
                    val missRateList = cacheMissIs.mu
                    var v = 0
                    while (v < u) {
                        missRate.set(v, missRateList[v])
                        v++
                    }
                }

                else -> line_error(mfilename(object : Any() {}),
                    "NC does not support sampling solution of the specified cache replacement policy.")
            }
        }

        else ->         /*refactor this code
            [~,missRate,~,~,lE] = cache_miss_spm(gamma, m, lambda);
            pij = cache_prob_spm(gamma, m, lE);
              */
            when (cache!!.replacementStrategy) {
                ReplacementStrategy.RR, ReplacementStrategy.FIFO -> {
                    val cacheMissSpm = cache_miss_spm(gamma, m, lambdaCell)
                    val missRateList = cacheMissSpm.mu
                    var v = 0
                    while (v < u) {
                        missRate.set(v, missRateList[v])
                        v++
                    }
                    cacheMissSpm.getlE() //never used
                    res.method = "rayint"
                }

                else -> line_error(mfilename(object : Any() {}),
                    "NC does not support approximate solution of the specified cache replacement policy.")
            }
    }

    for (r in 0..<sn.nclasses) {
        if (ch.hitclass.length() > r && ch.missclass.get(r) > -1 && ch.hitclass.get(r) > -1) {
            XN.set(ch.missclass.get(r).toInt(), XN.get(ch.missclass.get(r).toInt()) + missRate.get(r))
            XN.set(ch.hitclass.get(r).toInt(), XN.get(ch.hitclass.get(r).toInt()) + sourceRate.get(r) - missRate.get(r))
        }
    }

    val endTime = System.nanoTime()
    res.QN = QN
    res.UN = UN
    res.RN = RN
    res.TN = TN
    res.CN = CN
    res.XN = XN
    res.AN = AN
    res.WN = WN
    res.lG = lG //???
    res.runtime = (endTime - startTime) / 1000000000.0
    //res.iter = iter;
    return res
}

