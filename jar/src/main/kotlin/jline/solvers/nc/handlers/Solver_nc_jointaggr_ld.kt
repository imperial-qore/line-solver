package jline.solvers.nc.handlers

import jline.api.pfqn.ld.pfqn_ncld
import jline.api.sn.snGetDemandsChain
import jline.lang.NetworkStruct
import jline.lang.state.State
import jline.lang.state.ToMarginal.toMarginal
import jline.solvers.SolverOptions
import jline.solvers.nc.SolverNC
import jline.util.Utils
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath

/**
 * Computes load-dependent aggregated joint probabilities for the NC solver.
 * This is a port of the MATLAB solver_nc_jointaggr_ld.m function.
 *
 * @param sn Network structure
 * @param options Solver options
 * @return SolverNCJointReturn containing joint probability, normalization constant, and runtime
 */
fun solver_nc_jointaggr_ld(sn: NetworkStruct, options: SolverOptions): SolverNC.SolverNCJointReturn {
    val startTimeMillis = System.nanoTime()

    // Initialization
    val state = sn.state
    val S = sn.nservers
    val rates = sn.rates
    
    // Determine service times
    val ST = rates.copy()
    for (i in 0..<ST.numRows) {
        for (j in 0..<ST.numCols) {
            ST.set(i, j, 1.0 / ST.get(i, j))
        }
    }
    ST.removeNaN()

    val ret = snGetDemandsChain(sn)
    val Lchain = ret.Dchain
    val STchain = ret.STchain
    val alpha = ret.alpha
    val Nchain = ret.Nchain

    val M = STchain.numRows
    val K = STchain.numCols

    // Build load-dependent scaling matrix mu_chain
    val mu_chain = Matrix(M, Nchain.elementSum().toInt())
    for (ist in 0..<M) {
        if (Utils.isInf(S.get(ist))) { // Infinite server
            for (j in 0..<Nchain.elementSum().toInt()) {
                mu_chain.set(ist, j, (j + 1).toDouble())
            }
        } else {
            for (j in 0..<Nchain.elementSum().toInt()) {
                mu_chain.set(ist, j, FastMath.min((j + 1).toDouble(), S.get(ist)))
            }
        }
    }

    // Compute global normalization constant
    val Zchain = Nchain.copy()
    Zchain.fill(0.0)
    val lG = pfqn_ncld(Lchain, Nchain, Zchain, mu_chain, options).lG as Double

    // Compute joint probability in log space
    var lPr = 0.0
    for (ist in 0..<M) {
        val isf = sn.stationToStateful.get(ist).toInt()
        val marginalResult = toMarginal(sn, ist, state.get(sn.stateful.get(isf)), null, null, null, null, null)
        val nivec = marginalResult.nir
        
        // Convert to chain populations
        val nivec_chain = nivec.mult(sn.chains.transpose())
        
        // Extract parameters for station ist
        val Lchain_ist = Matrix.extractRows(Lchain, ist, ist + 1, null)
        val mu_chain_ist = Matrix.extractRows(mu_chain, ist, ist + 1, null)
        val Znivec_chain = nivec_chain.copy()
        Znivec_chain.fill(0.0)
        
        // Compute lF_i: normalization constant for Lchain at this station
        val lF_i = pfqn_ncld(Lchain_ist, nivec_chain, Znivec_chain, mu_chain_ist, options).lG as Double
        
        // Build alpha-scaled service times for station ist  
        val ST_alpha_ist = Matrix(1, sn.nclasses)
        for (k in 0..<sn.nclasses) {
            ST_alpha_ist.set(0, k, ST.get(ist, k) * alpha.get(ist, k))
        }
        val Znivec = nivec.copy()
        Znivec.fill(0.0)
        
        // Compute lg0_i: normalization constant for alpha-scaled service times
        val lg0_i = pfqn_ncld(ST_alpha_ist, nivec, Znivec, mu_chain_ist, options).lG as Double
        
        // Extract STchain for station ist
        val STchain_ist = Matrix.extractRows(STchain, ist, ist + 1, null)
        
        // Compute lG0_i: normalization constant for STchain
        val lG0_i = pfqn_ncld(STchain_ist, nivec_chain, Znivec_chain, mu_chain_ist, options).lG as Double
        
        // Accumulate log probability: lF_i + (lg0_i - lG0_i)
        lPr += lF_i + (lg0_i - lG0_i)
    }
    
    // Final probability computation
    val Pr = FastMath.exp(lPr - lG)
    val G = FastMath.exp(lG)
    
    val endTimeMillis = System.nanoTime()
    val runtime = (endTimeMillis - startTimeMillis) / 1000000000.0

    return SolverNC.SolverNCJointReturn(Pr, G, lG, runtime)
}