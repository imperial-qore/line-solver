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

fun solver_nc_joint(sn: NetworkStruct, options: SolverOptions): SolverNC.SolverNCJointReturn {
    val state = sn.state
    val S = sn.nservers
    val rates = sn.rates
    val ST = rates.copy()
    for (i in 0..<ST.numRows) {
        for (j in 0..<ST.numCols) {
            ST.set(i, j, 1.0 / ST.get(i, j))
        }
    }
    ST.removeNaN()
    val ret = snGetDemandsChain(sn)
    val STchain = ret.STchain
    val Vchain = ret.Vchain
    val alpha = ret.alpha
    val Nchain = ret.Nchain
    val startTimeMillis = System.nanoTime()
    val M = STchain.numRows
    val K = ST.numCols
    var Lchain = Matrix(0, K)

    var mu_chain = Matrix(0, Nchain.elementSum().toInt())
    for (i in 0..<M) {
        val tmp = Matrix(1, Nchain.elementSum().toInt())
        if (Utils.isInf(S.get(i))) {
            for (j in 0..<tmp.length()) {
                tmp.set(j, (j + 1).toDouble())
            }
        } else {
            for (j in 0..<tmp.length()) {
                tmp.set(j, FastMath.min((j + 1).toDouble(), S.get(i)))
            }
        }
        mu_chain = Matrix.concatRows(mu_chain, tmp, null)

        val tmp1 =
            Matrix.extractRows(STchain, i, i + 1, null).elementMult(Matrix.extractRows(Vchain, i, i + 1, null), null)
        Lchain = Matrix.concatRows(Lchain, tmp1, null)
    }

    val Z_tmp = Nchain.copy()
    Z_tmp.fill(0.0)
    val lG = pfqn_ncld(Lchain, Nchain, Z_tmp, mu_chain, options).lG
    var lPr = 0.0

    for (i in 0..<M) {
        val isf = sn.stationToStateful.get(i).toInt()
        val ret1 = toMarginal(sn, i, state.get(sn.stateful.get(isf)), null, null, null, null, null)
        val nivec = ret1.nir
        val nivec_chain = nivec.mult(sn.chains.transpose())
        val Lchain_i = Matrix.extractRows(Lchain, i, i + 1, null)
        val ST_i = Matrix.extractRows(ST, i, i + 1, null)
        val alpha_i = Matrix.extractRows(alpha, i, i + 1, null)
        val STchain_i = Matrix.extractRows(STchain, i, i + 1, null)
        val mu_chain_i = Matrix.extractRows(mu_chain, i, i + 1, null)
        val Zvec_chain = nivec_chain.copy()
        Zvec_chain.fill(0.0)
        val Zvec = nivec.copy()
        Zvec.fill(0.0)

        val lF_i = pfqn_ncld(Lchain_i, nivec_chain, Zvec_chain, mu_chain_i, options).lG
        val lg0_i = pfqn_ncld(ST_i.elementMult(alpha_i, null), nivec, Zvec, mu_chain_i, options).lG
        val lG0_i = pfqn_ncld(STchain_i, nivec_chain, Zvec_chain, mu_chain_i, options).lG
        lPr = lPr + lF_i + (lg0_i - lG0_i)
    }
    val Pr = FastMath.exp(lPr - lG)
    val endTimeMillis = System.nanoTime()
    val runtime = (endTimeMillis - startTimeMillis) / 1000000000.0
    val G = FastMath.exp(lG)
    return SolverNC.SolverNCJointReturn(Pr, G, lG, runtime)
}