package jline.solvers.mam.handlers

import jline.lang.NetworkStruct
import jline.lang.constant.NodeType
import jline.util.matrix.Matrix

data class MetricsResult(
    val QN: Matrix,
    val UN: Matrix,
    val RN: Matrix,
    val TN: Matrix,
    val CN: Matrix,
    val XN: Matrix
)

fun solver_mam_metrics(
    sn: NetworkStruct,
    inapResult: INAPResult,
    rcat: RCATModel
): MetricsResult {
    val M = sn.nstations
    val K = sn.nclasses
    val processMap = rcat.processMap
    val pi = inapResult.pi
    val N = rcat.N

    val QN = Matrix(M, K)
    val UN = Matrix(M, K)
    val RN = Matrix(M, K)
    val TN = Matrix(M, K)

    for (ist in 0 until M) {
        for (r in 0 until K) {
            val p = processMap.get(ist, r).toInt()
            if (p >= 0 && p < pi.size && pi[p].getNumCols() > 0) {
                val Np = N[p]
                val piP = pi[p]

                var queueLength = 0.0
                for (n in 0 until Np) {
                    queueLength += n * piP.get(0, n)
                }
                QN.set(ist, r, queueLength)

                val utilization = 1.0 - piP.get(0, 0)
                UN.set(ist, r, utilization)

                val muIr = sn.rates.get(ist, r)
                if (!muIr.isNaN() && muIr > 0) {
                    TN.set(ist, r, muIr * utilization)
                }
            }
        }
    }

    for (ist in 0 until M) {
        for (r in 0 until K) {
            val tput = TN.get(ist, r)
            if (tput > 0) {
                RN.set(ist, r, QN.get(ist, r) / tput)
            } else {
                RN.set(ist, r, 0.0)
            }
        }
    }

    val CN = Matrix(1, K)
    val XN = Matrix(1, K)

    for (r in 0 until K) {
        val njobs = sn.njobs.get(r)

        if (njobs.isInfinite()) {
            for (ist in 0 until M) {
                val nodeIdx = sn.stationToNode.get(ist).toInt()
                if (sn.nodetype[nodeIdx] == NodeType.Source) {
                    val rate = sn.rates.get(ist, r)
                    if (!rate.isNaN()) {
                        XN.set(0, r, rate)
                    }
                    break
                }
            }
            var totalRespTime = 0.0
            for (ist in 0 until M) {
                totalRespTime += RN.get(ist, r)
            }
            CN.set(0, r, totalRespTime)
        } else {
            val refst = sn.refstat.get(r).toInt()
            if (refst >= 0 && refst < M) {
                val xr = TN.get(refst, r)
                XN.set(0, r, xr)
                if (xr > 0) {
                    CN.set(0, r, njobs / xr)
                }
            }
        }
    }

    return MetricsResult(QN, UN, RN, TN, CN, XN)
}
