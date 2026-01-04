package jline.api.sn

import jline.io.Ret
import jline.lang.NetworkStruct
import jline.GlobalConstants
import jline.GlobalConstants.Inf
import jline.util.Utils
import jline.util.matrix.Matrix

/**
 * Calculate new queueing network parameters after aggregating classes into chains
 *
 * @param sn - NetworkStruct object for the queueing network model
 * @return chain parameters
 */

fun snGetDemandsChain(sn_in: NetworkStruct): Ret.snGetDemands {
    val sn: NetworkStruct = sn_in.copy();
    val M = sn.nstations
    val K = sn.nclasses
    val C = sn.nchains
    val N = sn.njobs

    val scv = sn.scv
    scv.apply(Double.NaN, 1.0, "equal")

    val ST = Matrix(0, 0)
    sn.rates.divide(1.0, ST, false)
    ST.removeNaN()

    val alpha = Matrix(M, K)
    val Vchain = Matrix(M, C)
    for (c in 0..<C) {
        val inchain = sn.inchain[c]
        if (sn.refclass[0, c] > -1) {
            for (i in 0..<M) {
                //Vchain(i,c) = sum(sn.visits{c}(i,inchain)) / sum(sn.visits{c}(sn.refstat(inchain(1)),sn.refclass(c)));
                val visits = sn.visits[c]
                var res = 0.0
                val iIdx = sn.stationToStateful[i].toInt()
                for (col in 0..<inchain!!.numCols) res += visits!![iIdx, inchain[0, col].toInt()]
                Vchain[i, c] = res / visits!![sn.stationToStateful[sn.refstat[inchain.value()
                    .toInt(), 0].toInt()].toInt(), sn.refclass[0, c].toInt()]
                //alpha(i,k) = alpha(i,k) + sn.visits{c}(i,k) / sum(sn.visits{c}(i,inchain));
                for (col in 0..<inchain.numCols) {
                    val k = inchain[0, col].toInt()
                    alpha[i, k] = alpha[i, k] + visits[i, k] / res
                }
            }
        } else {
            for (i in 0..<M) {
                //Vchain(i,c) = sum(sn.visits{c}(i,inchain)) / sum(sn.visits{c}(sn.refstat(inchain(1)),inchain));
                val visits = sn.visits[c]
                var res1 = 0.0
                var res2 = 0.0
                val refIdx = sn.stationToStateful[sn.refstat[inchain!!.value().toInt(), 0].toInt()].toInt()
                val iIdx = sn.stationToStateful[i].toInt()
                for (col in 0..<inchain.numCols) {
                    val idx = inchain[0, col].toInt()
                    res1 += visits!![iIdx, idx]
                    res2 += visits[refIdx, idx]
                }
                Vchain[i, c] = res1 / res2
                //alpha(i,k) = alpha(i,k) + sn.visits{c}(i,k) / sum(sn.visits{c}(i,inchain));
                for (col in 0..<inchain.numCols) {
                    val k = inchain[0, col].toInt()
                    alpha[i, k] = alpha[i, k] + visits!![iIdx, k] / res1
                }
            }
        }
    }

    Vchain.apply(Inf, 0.0, "equal")
    Vchain.apply(Double.NaN, 0.0, "equal")
    for (c in 0..<C) {
        val vchainRef = Vchain[sn.refstat[sn.inchain[c]!!.value().toInt(), 0].toInt(), c]
        for (i in Vchain.colIndexes[c]..<Vchain.colIndexes[c + 1]) Vchain.nonZeroValues[i] /= vchainRef
    }
    alpha.apply(Inf, 0.0, "equal")
    alpha.apply(Double.NaN, 0.0, "equal")
    alpha.apply(GlobalConstants.Zero, 0.0, "less")

    val Lchain = Matrix(M, C)
    val STchain = Matrix(M, C)
    val SCVchain = Matrix(M, C)
    val Nchain = Matrix(1, C)
    val refstatchain = Matrix(C, 1)
    for (c in 0..<C) {
        val inchain = sn.inchain[c]
        //Nchain(c) = sum(N(inchain)); isOpenChain = any(isinf(N(inchain)));
        var isOpenChain = false
        var sum = 0.0
        for (col in 0..<inchain!!.numCols) {
            sum += N[inchain[0, col].toInt()]
            if (Utils.isInf(sum)) {
                isOpenChain = true
                break
            }
        }
        Nchain[0, c] = sum

        for (i in 0..<M) {
            sum = 0.0
            if (isOpenChain && i.toDouble() == sn.refstat[inchain.value().toInt(), 0]) {
                //STchain(i,c) = 1.0 / sumfinite(sn.rates(i,inchain));
                for (col in 0..<inchain.numCols) {
                    val rateValue = sn.rates[i, inchain[0, col].toInt()]
                    if (java.lang.Double.isFinite(rateValue)) sum += rateValue
                }
                STchain[i, c] = 1.0 / sum
            } else {
                //STchain(i,c) = ST(i,inchain) * alpha(i,inchain)';
                for (col in 0..<inchain.numCols) {
                    val idx = inchain[0, col].toInt()
                    sum += ST[i, idx] * alpha[i, idx]
                }
                STchain[i, c] = sum
            }
            Lchain[i, c] = Vchain[i, c] * STchain[i, c]
            //alphachain = sum(alpha(i,inchain(isfinite(SCV(i,inchain))))');
            var alphachain = 0.0
            for (col in 0..<inchain.numCols) {
                val idx = inchain[0, col].toInt()
                val scvValue = scv[i, idx]
                if (java.lang.Double.isFinite(scvValue)) alphachain += alpha[i, idx]
            }
            if (alphachain > GlobalConstants.Zero) {
                sum = 0.0
                for (col in 0..<inchain.numCols) {
                    val idx = inchain[0, col].toInt()
                    sum += scv[i, idx] * alpha[i, idx]
                }
                SCVchain[i, c] = sum / alphachain
            }
        }
        refstatchain[c, 0] = sn.refstat[inchain.value().toInt(), 0]
        for (col in 1..<inchain.numCols) {
            val classIdx = inchain[0, col].toInt()
            if (sn.refstat[classIdx, 0] != refstatchain[c, 0]) throw RuntimeException("Class have different reference station")
        }
    }
    Lchain.apply(Inf, 0.0, "equal")
    Lchain.apply(Double.NaN, 0.0, "equal")
    STchain.apply(Inf, 0.0, "equal")
    STchain.apply(Double.NaN, 0.0, "equal")
    return Ret.snGetDemands(Lchain, STchain, Vchain, alpha, Nchain, SCVchain, refstatchain)
}
/**
 * Stochastic network GetDemandsChain algorithms
 */
@Suppress("unused")
class SngetdemandschainAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}