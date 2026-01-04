/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.nc.analyzers

import jline.lang.NetworkStruct
import jline.GlobalConstants
import jline.solvers.SolverOptions
import jline.solvers.nc.NCResult
import jline.solvers.nc.handlers.solver_nc
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath

fun solver_nc_analyzer(sn: NetworkStruct, options: SolverOptions): NCResult {
    val Tstart = System.nanoTime()
    val nservers = sn.nservers
    val nserversFinite = nservers.copy()
    nserversFinite.removeInfinity()
    if (nserversFinite.elementMax() > 1 && options.method == "exact") {
        throw RuntimeException("NC solver cannot provide exact solutions for open or mixed queueing networks. Remove the 'exact' option.")
    }
    val snfloor: NetworkStruct = sn.copy()
    val snceil: NetworkStruct = sn.copy()

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
        val retfloor = solver_nc(snfloor, options)
        res.runtime = retfloor.runtime
        val retceil = solver_nc(snceil, options)
        res.runtime += retceil.runtime
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
        val ret = solver_nc(sn, options)
        res.QN = ret.Q.copy()
        res.UN = ret.U.copy()
        res.RN = ret.R.copy()
        res.TN = ret.T.copy()
        res.XN = ret.X.copy()
        
        // Calculate CN using Little's Law: C(k) = N(k) / X(k) for each class
        res.CN = Matrix(1, sn.nclasses)
        for (k in 0..<sn.nclasses) {
            res.CN.set(0, k, sn.njobs.get(0, k) / res.XN.get(0, k))
        }
        
        res.lG = ret.lG
        res.it = ret.it
        res.runtime = ret.runtime
        res.method = ret.method
    }
    res.runtime = (System.nanoTime() - Tstart).toDouble() / 1000000000.0

    return res
}
