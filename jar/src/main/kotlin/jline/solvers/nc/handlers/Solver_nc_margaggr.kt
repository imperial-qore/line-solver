package jline.solvers.nc.handlers

import jline.api.mam.map_lambda
import jline.api.pfqn.ld.pfqn_ncld
import jline.api.sn.snGetDemandsChain
import jline.lang.NetworkStruct
import jline.lang.state.State
import jline.lang.state.ToMarginal.toMarginal
import jline.solvers.SolverOptions
import jline.solvers.nc.SolverNC
import jline.util.Utils
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import org.apache.commons.math3.util.FastMath
import kotlin.math.ln

/**
 * Computes aggregated marginal probabilities for the NC solver.
 * This is a port of the MATLAB solver_nc_margaggr.m function.
 *
 * @param sn Network structure
 * @param options Solver options
 * @param lG Log normalization constant (optional)
 * @return SolverNCMargReturn containing marginal probabilities and normalization constant
 */
fun solver_nc_margaggr(sn: NetworkStruct, options: SolverOptions, lG: Double?): SolverNC.SolverNCMargReturn {
    val startTimeMillis = System.nanoTime()
    
    val M = sn.nstations
    val K = sn.nclasses
    val state = sn.state
    val S = sn.nservers
    val NK = sn.njobs.transpose() // Initial population per class
    val C = sn.nchains
    val PH = sn.proc

    // Determine service times
    val ST = Matrix(M, K)
    for (k in 0..<K) {
        for (ist in 0..<M) {
            val processCell: MatrixCell = PH.get(sn.stations.get(ist))!!.get(sn.jobclasses.get(k))!! as MatrixCell
            val lambda = map_lambda(processCell)
            ST.set(ist, k, 1.0 / lambda)
        }
    }
    ST.removeNaN()

    val ret = snGetDemandsChain(sn)
    val Lchain = ret.Dchain
    val STchain = ret.STchain
    val Nchain = ret.Nchain

    // Build visit matrix V
    val V = Matrix(sn.nstations, sn.nclasses)
    for (c in 0..<sn.nchains) {
        val inchain = sn.inchain.get(c)
        for (ist in 0..<sn.nstations) {
            for (col in 0..<inchain!!.numCols) {
                val k = inchain.get(0, col).toInt()
                V.set(ist, k, sn.visits.get(c)?.get(ist, k) ?: 0.0)
            }
        }
    }

    val M_chain = STchain.numRows
    
    // Build mu matrix for load-dependent scaling
    val mu = Matrix(M_chain, Nchain.elementSum().toInt())
    for (ist in 0..<M_chain) {
        if (Utils.isInf(S.get(ist))) { // Infinite server
            for (j in 0..<Nchain.elementSum().toInt()) {
                mu.set(ist, j, (j + 1).toDouble())
            }
        } else {
            for (j in 0..<Nchain.elementSum().toInt()) {
                mu.set(ist, j, FastMath.min((j + 1).toDouble(), S.get(ist)))
            }
        }
    }

    // Compute normalization constant if not provided
    val lG_val = lG ?: run {
        val Zchain_tmp = Nchain.copy()
        Zchain_tmp.fill(0.0)
        pfqn_ncld(Lchain, Nchain, Zchain_tmp, mu, options).lG as Double
    }
    val G = FastMath.exp(lG_val)

    // Compute marginal probabilities
    val lPr = Matrix(sn.nstations, 1)
    for (ist in 0..<sn.nstations) {
        val ind = sn.stationToNode.get(ist).toInt()
        val isf = sn.stationToStateful.get(ist).toInt()
        val marginalResult = toMarginal(sn, ind, state.get(sn.stateful.get(isf)), null, null, null, null, null)
        val nivec = marginalResult.nir
        
        if (nivec.elementMin() < 0) { // User flags that state should be ignored
            lPr.set(ist, Double.NaN)
        } else {
            // Compute reduced system without station ist
            var Lchain_minus_i = Matrix(0, Lchain.numCols)
            var mu_minus_i = Matrix(0, mu.numCols)
            for (j in 0..<sn.nstations) {
                if (j != ist) {
                    val Lchain_row = Matrix.extractRows(Lchain, j, j + 1, null)
                    val mu_row = Matrix.extractRows(mu, j, j + 1, null)
                    Lchain_minus_i = Matrix.concatRows(Lchain_minus_i, Lchain_row, null)
                    mu_minus_i = Matrix.concatRows(mu_minus_i, mu_row, null)
                }
            }
            
            // Convert to chain populations
            val nivec_chain = nivec.mult(sn.chains.transpose())
            val Nchain_minus_i = Nchain.copy()
            for (j in 0..<Nchain_minus_i.length()) {
                Nchain_minus_i.set(j, Nchain_minus_i.get(j) - nivec_chain.get(j))
            }
            val Zchain_minus_i = Nchain.copy()
            Zchain_minus_i.fill(0.0)
            
            // Compute normalization constant for reduced system
            val lG_minus_i = pfqn_ncld(Lchain_minus_i, Nchain_minus_i, Zchain_minus_i, mu_minus_i, options).lG as Double
            
            // Compute local normalization constant for station ist
            val ST_ist = Matrix(1, K)
            val V_ist = Matrix(1, K)
            for (k in 0..<K) {
                ST_ist.set(0, k, ST.get(ist, k))
                V_ist.set(0, k, V.get(ist, k))
            }
            val ST_V_ist = ST_ist.elementMult(V_ist)
            val mu_ist = Matrix.extractRows(mu, ist, ist + 1, null)
            val Znivec = nivec.copy()
            Znivec.fill(0.0)
            
            val lF_i = pfqn_ncld(ST_V_ist, nivec, Znivec, mu_ist, options).lG as Double
            
            // Marginal probability in log space
            lPr.set(ist, lF_i + lG_minus_i - lG_val)
        }
    }
    
    val Pr = lPr.copy()
    for (i in 0..<Pr.numRows) {
        for (j in 0..<Pr.numCols) {
            Pr.set(i, j, FastMath.exp(Pr.get(i, j)))
        }
    }
    Pr.removeNaN()
    
    val endTimeMillis = System.nanoTime()
    val runtime = (endTimeMillis - startTimeMillis) / 1000000000.0
    
    return SolverNC.SolverNCMargReturn(Pr, G, lG?.toDouble() ?: Double.MIN_VALUE, runtime)
}