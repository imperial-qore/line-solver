/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.mva.analyzers

import jline.api.cache.cache_gamma_lp
import jline.api.cache.cache_mva
import jline.api.cache.cache_prob_fpi
import jline.api.cache.cache_ttl_lrum
import jline.io.line_error
import jline.io.mfilename
import jline.lang.NetworkStruct
import jline.lang.constant.NodeType
import jline.lang.constant.ReplacementStrategy
import jline.lang.nodeparam.CacheNodeParam
import jline.lang.nodes.Cache
import jline.solvers.SolverOptions
import jline.solvers.mva.MVAResult
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath
import java.lang.Double
import kotlin.Any
import kotlin.Array

/**
 * MVA Analyzer class for non-rentrant caches
 */
fun solver_mva_cache_analyzer(sn: NetworkStruct, options: SolverOptions): MVAResult {
    val res = MVAResult()
    val startTime = System.nanoTime()
    val method = options.method
    val QN = Matrix(sn.nnodes, sn.nclasses)
    val UN = Matrix(sn.nnodes, sn.nclasses)
    val RN = Matrix(sn.nnodes, sn.nclasses)
    var TN = Matrix(sn.nnodes, sn.nclasses)
    val CN = Matrix(1, sn.nclasses)
    val AN = Matrix(sn.nnodes, sn.nclasses)
    val WN = Matrix(sn.nnodes, sn.nclasses)
    val XN = Matrix(1, sn.nclasses)
    val lG = Double.NaN
    val iter = 1

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
    val h = m.length() // Number of lists
    val u = sn.nclasses // Number of users
    val lambda = Array(u) { Matrix(n, h + 1) }
    for (v in 0..<u) {
        for (k in 0..<n) {
            for (l in 0..<h + 1) {
                if (ch.pread.getOrDefault(v, null) != null) {
                    lambda[v].set(k, l, sourceRate.get(v) * ch.pread.get(v)!!.get(k))
                }
            }
        }
    }

    val R = ch.accost
    val gamma = cache_gamma_lp(lambda, R).gamma
    var pij: Matrix? = null
    var missRate: Matrix?
    missRate = Matrix(1, u)
    when (method) {
        "exact" -> when (cache!!.replacementStrategy) {
            ReplacementStrategy.RR, ReplacementStrategy.FIFO -> {
                pij = cache_mva(gamma, m).pij
                val newPij = Matrix(pij.numRows, 1 + pij.numCols)
                var i = 0
                while (i < pij.numRows) {
                    var j = 0
                    while (j < newPij.numCols) {
                        if (j < 1) {
                            newPij.set(i, j, FastMath.abs(1 - pij.sumRows(i)))
                        } else {
                            newPij.set(i, j, pij.get(i, j - 1))
                        }
                        j++
                    }
                    i++
                }
                pij = newPij
                missRate = Matrix(1, u)
                var v = 0
                while (v < u) {
                    missRate.set(v,
                        Matrix.extractColumn(lambda[v], 0, null).transpose().mult(Matrix.extractColumn(pij, 0, null))
                            .get(0))
                    v++
                }
                line_error(mfilename(object : Any() {}),
                    "MVA does not support exact solution of the specified cache replacement policy.")
            }

            else -> line_error(mfilename(object : Any() {}),
                "MVA does not support exact solution of the specified cache replacement policy.")
        }

        else -> {
            when (cache!!.replacementStrategy) {
                ReplacementStrategy.RR, ReplacementStrategy.FIFO -> pij = cache_prob_fpi(gamma, m) // FPI method
                ReplacementStrategy.LRU ->                         // doesn't work for large caches e.g.n = 10000 m = 2000
                    @Suppress("UNCHECKED_CAST")
                    pij = cache_ttl_lrum(lambda as Array<Matrix?>,
                        m) // without considering different graph of different items linear
                else -> line_error(mfilename(object : Any() {}),
                    "MVA does not support approximate solution of the specified cache replacement policy.")
            }
            missRate = Matrix(1, u)
            var v = 0
            while (v < u) {
                missRate.set(v,
                    Matrix.extractColumn(lambda[v], 0, null).transpose().mult(Matrix.extractColumn(pij, 0, null))
                        .get(0))
                v++
            }
        }
    }
    if (method !== "ttl.map") {
        for (r in 0..<sn.nclasses) {
            if (ch.hitclass.length() > r && ch.missclass.get(r) > -1 && ch.hitclass.get(r) > -1) {
                XN.set(ch.missclass.get(r).toInt(), XN.get(ch.missclass.get(r).toInt()) + missRate.get(r))
                XN.set(ch.hitclass.get(r).toInt(),
                    XN.get(ch.hitclass.get(r).toInt()) + sourceRate.get(r) - missRate.get(r))
            }
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
    res.logNormConstAggr = lG
    res.runtime = (endTime - startTime) / 1000000000.0
    res.iter = iter
    
    // Set the actual method used
    res.method = when {
        method == "exact" -> "exact"
        else -> {
            when (cache!!.replacementStrategy) {
                ReplacementStrategy.RR, ReplacementStrategy.FIFO -> "fpi"
                ReplacementStrategy.LRU -> "ttl"
                else -> method
            }
        }
    }
    
    return res
}

