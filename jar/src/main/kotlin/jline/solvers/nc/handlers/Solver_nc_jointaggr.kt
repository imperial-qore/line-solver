package jline.solvers.nc.handlers

import jline.api.pfqn.ld.pfqn_ncld
import jline.api.sn.snGetDemandsChain
import jline.lang.NetworkStruct
import jline.lang.state.State
import jline.lang.state.ToMarginal.toMarginal
import jline.solvers.SolverOptions
import jline.solvers.nc.SolverNC
import jline.solvers.nc.handlers.solver_nc
import jline.util.Utils
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath

/**
 * Computes aggregated joint probabilities for the NC solver.
 * This is a port of the MATLAB solver_nc_jointaggr.m function.
 *
 * @param sn Network structure
 * @param options Solver options
 * @return SolverNCJointReturn containing joint probability and normalization constant
 */
fun solver_nc_jointaggr(sn: NetworkStruct, options: SolverOptions): SolverNC.SolverNCJointReturn {
    val startTimeMillis = System.nanoTime()

    // Calculate visit matrix (cellsum equivalent)
    var V = Matrix(sn.nstateful, sn.nclasses)
    for (i in 0..<sn.visits.size) {
        V = V.add(1.0, sn.visits.get(i))
    }

    val M = sn.nstations
    val C = sn.nchains
    val Nchain = Matrix(1, C)

    // Build chain populations
    for (c in 0..<C) {
        val inchain = sn.inchain.get(c)
        var chainPop = 0.0
        for (col in 0..<inchain!!.numCols) {
            val k = inchain.get(0, col).toInt()
            chainPop += sn.njobs.get(k)
        }
        Nchain.set(c, chainPop)
    }

    // Transform everything into a load-dependent model
    val nservers = sn.nservers
    val mu = Matrix(M, Nchain.elementSum().toInt())
    for (ist in 0..<M) {
        if (Utils.isInf(nservers.get(ist))) { // Infinite server
            for (j in 0..<Nchain.elementSum().toInt()) {
                mu.set(ist, j, (j + 1).toDouble())
            }
        } else {
            for (j in 0..<Nchain.elementSum().toInt()) {
                mu.set(ist, j, FastMath.min((j + 1).toDouble(), nservers.get(ist)))
            }
        }
    }
    val state = sn.state

    val lG: Double
    val ST: Matrix

    when (options.method) {
        "exact" -> {
            val ret = snGetDemandsChain(sn)
            val Lchain = ret.Dchain
            val Nchain_ret = ret.Nchain
            val Zchain = Nchain_ret.copy()
            Zchain.fill(0.0)
            lG = pfqn_ncld(Lchain, Nchain_ret, Zchain, mu, options).lG as Double
            
            // Build service time matrix
            ST = sn.rates.copy()
            for (i in 0..<ST.numRows) {
                for (j in 0..<ST.numCols) {
                    ST.set(i, j, 1.0 / ST.get(i, j))
                }
            }
            ST.removeNaN()
        }
        else -> {
            // Use general NC solver - note: this may not fully consider LD model transformation
            val ncResult = solver_nc(sn, options)
            lG = ncResult.lG
            ST = ncResult.STeff
        }
    }
    
    val G = FastMath.exp(lG)
    var lPr = 0.0

    // Sum over all stations
    for (ist in 0..<M) {
        val isf = sn.stationToStateful.get(ist).toInt()
        val ind = sn.stationToNode.get(ist).toInt()
        val marginalResult = toMarginal(sn, ind, state.get(sn.stateful.get(isf)), null, null, null, null, null)
        val nivec = marginalResult.nir
        
        // Get unique rows (equivalent to MATLAB unique(nivec,'rows'))
        val uniqueNivec = getUniqueRows(nivec)
        
        for (row in 0..<uniqueNivec.numRows) {
            val currentNivec = Matrix.extractRows(uniqueNivec, row, row + 1, null)
            val nivec_chain = currentNivec.mult(sn.chains.transpose())
            
            // Check if any population is positive
            var hasPositivePop = false
            for (j in 0..<nivec_chain.length()) {
                if (nivec_chain.get(j) > 0) {
                    hasPositivePop = true
                    break
                }
            }
            
            if (hasPositivePop) {
                // Build service time matrix for this station
                val ST_ist = Matrix(1, sn.nclasses)
                val V_ist = Matrix(1, sn.nclasses)
                for (k in 0..<sn.nclasses) {
                    ST_ist.set(0, k, ST.get(ist, k))
                    V_ist.set(0, k, V.get(ist, k))
                }
                val ST_V_ist = ST_ist.elementMult(V_ist)
                val mu_ist = Matrix.extractRows(mu, ist, ist + 1, null)
                val Znivec = currentNivec.copy()
                Znivec.fill(0.0)
                
                val lF_i = pfqn_ncld(ST_V_ist, currentNivec, Znivec, mu_ist, options).lG as Double
                lPr += lF_i
            }
        }
    }
    
    lPr -= lG
    val Pr = FastMath.exp(lPr)

    val endTimeMillis = System.nanoTime()
    val runtime = (endTimeMillis - startTimeMillis) / 1000000000.0

    return SolverNC.SolverNCJointReturn(Pr, G, lG, runtime)
}

/**
 * Helper function to get unique rows from a matrix (equivalent to MATLAB's unique(matrix,'rows'))
 * @param matrix Input matrix
 * @return Matrix containing unique rows
 */
private fun getUniqueRows(matrix: Matrix): Matrix {
    val uniqueRowsList = mutableListOf<Matrix>()
    val seenRows = mutableSetOf<String>()
    
    for (i in 0..<matrix.numRows) {
        val row = Matrix.extractRows(matrix, i, i + 1, null)
        val rowString = (0..<row.numCols).map { row.get(0, it).toString() }.joinToString(",")
        
        if (!seenRows.contains(rowString)) {
            seenRows.add(rowString)
            uniqueRowsList.add(row)
        }
    }
    
    if (uniqueRowsList.isEmpty()) {
        return Matrix(0, matrix.numCols)
    }
    
    var result = uniqueRowsList[0]
    for (i in 1..<uniqueRowsList.size) {
        result = Matrix.concatRows(result, uniqueRowsList[i], null)
    }
    
    return result
}