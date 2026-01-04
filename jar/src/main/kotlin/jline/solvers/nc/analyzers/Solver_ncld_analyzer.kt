/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.nc.analyzers

import jline.lang.NetworkStruct
import jline.GlobalConstants
import jline.solvers.SolverOptions
import jline.solvers.nc.NCResult
import jline.solvers.nc.handlers.solver_ncld
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath
import java.lang.Double
import kotlin.RuntimeException

fun solver_ncld_analyzer(sn: NetworkStruct, options: SolverOptions): NCResult {
    val Tstart = System.nanoTime()
    val nservers = sn.nservers
    val nserversFinite = nservers.copy()
    nserversFinite.removeInfinity()
    if (nserversFinite.elementMax() > 1 && sn.njobs.hasInfinite() && options.method == "exact") {
        throw RuntimeException("NC solver cannot provide exact solutions for open or mixed queueing networks. Remove the 'exact' option.")
    }
    val snfloor: NetworkStruct
    val snceil: NetworkStruct
    snfloor = sn.copy()
    snceil = sn.copy()

    snfloor.njobs = sn.njobs.copy()
    snceil.njobs = sn.njobs.copy()

    val eta = Matrix(sn.njobs.numRows, sn.njobs.numCols)
    var nonIntegerJob = false
    for (i in 0..<eta.numRows) {
        for (j in 0..<eta.numCols) {
            snfloor.njobs.set(i, j, FastMath.floor(sn.njobs.get(i, j)))
            snceil.njobs.set(i, j, FastMath.ceil(sn.njobs.get(i, j)))
            eta.set(i, j, FastMath.abs(sn.njobs.get(i, j) - snfloor.njobs.get(i, j)))
            if (eta.get(i, j) > GlobalConstants.FineTol) {
                nonIntegerJob = true
            }
        }
    }
    var res = NCResult()
    if (nonIntegerJob) {
        if (options.method == "exact") {
            throw RuntimeException("NC load-dependent solver cannot provide exact solutions for fractional populations.")
        }
        val retfloor = solver_ncld(snfloor, options)
        val retceil = solver_ncld(snceil, options)
        res.QN = retfloor.Q.add(1.0, eta.elementMult(retceil.Q.sub(retfloor.Q), null))
        res.UN = retfloor.U.add(1.0, eta.elementMult(retceil.U.sub(retfloor.U), null))
        res.RN = retfloor.R.add(1.0, eta.elementMult(retceil.R.sub(retfloor.R), null))
        res.TN = retfloor.T.add(1.0, eta.elementMult(retceil.T.sub(retfloor.T), null))
        res.XN = retfloor.X.add(1.0, eta.elementMult(retceil.X.sub(retfloor.X), null))
        
        // Calculate CN using Little's Law: C(k) = N(k) / X(k) for each class
        res.CN = Matrix(1, sn.nclasses)
        for (k in 0..<sn.nclasses) {
            res.CN.set(0, k, sn.njobs.get(0, k) / res.XN.get(0, k))
        }
        
        // Handle lG interpolation: lG = lGf + sum(eta .* (lGf - lGc))
        res.lG = retfloor.lG + eta.elementSum() * (retfloor.lG - retceil.lG)
        res.it = retfloor.it + retceil.it
        res.method = retceil.method
    } else {
        val ret = solver_ncld(sn, options)
        res.QN = ret.Q
        res.UN = ret.U
        res.RN = ret.R
        res.TN = ret.T
        res.XN = ret.X
        
        // Calculate CN using Little's Law: C(k) = N(k) / X(k) for each class
        res.CN = Matrix(1, sn.nclasses)
        for (k in 0..<sn.nclasses) {
            res.CN.set(0, k, sn.njobs.get(0, k) / res.XN.get(0, k))
        }
        
        res.lG = ret.lG
        res.it = ret.it
        res.method = ret.method
    }
    res.runtime = (System.nanoTime() - Tstart).toDouble() / 1000000000.0
    return res
}
